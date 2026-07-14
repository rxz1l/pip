library(tidyverse)
library(data.table)
library(clusterProfiler)
library(enrichplot)
library(ggridges)
library(GseaVis)
library(ggplot2)
library(grDevices)

# ==============================================================================
# 1. 保持原样：格式化通路名称函数 [cite: 12, 13]
# ==============================================================================
format_pathway_name <- function(x) {
  db <- ifelse(grepl("^KEGG_", x), "KEGG",
        ifelse(grepl("^HALLMARK_", x), "HALLMARK", "OTHER"))
  term <- gsub("^(KEGG_|HALLMARK_)", "", x)
  term <- gsub("_", " ", term)
  term <- stringr::str_to_lower(term)
  term <- stringr::str_to_title(term)
  term <- gsub("\\bTca\\b", "TCA", term)
  term <- gsub("\\bAtp\\b", "ATP", term)
  term <- gsub("\\bDna\\b", "DNA", term)
  term <- gsub("\\bRna\\b", "RNA", term)
  term <- paste0("[", db, "] ", term)
  return(term)
}

# ==============================================================================
# 2. 修改的部分：重写数据准备函数 (不再读取外部 GSEA 结果文件)
# ==============================================================================
process_gsea_data <- function(deg_path, target_ids, species = "hs") {
  
  # ==============================================================================
  # 1. 数据准备：差异表处理 (保持稳定逻辑)
  # ==============================================================================
  gene_data <- read.table(deg_path, header = TRUE, sep = "\t", check.names = FALSE)
  
  # 基础过滤：去空、去NA
  gene_data <- gene_data[!is.na(gene_data$symbol) & gene_data$symbol != "" & !is.na(gene_data$log2FC), ]
  
  # 排序并去重 (保留绝对值最大的 log2FC)
  gene_data <- gene_data[order(abs(gene_data$log2FC), decreasing = TRUE), ]
  gene_data <- gene_data[!duplicated(gene_data$symbol), ]
  
  # 构建 genelist
  genelist <- gene_data$log2FC
  names(genelist) <- as.character(gene_data$symbol)
  
  # 增加微小扰动防止 GSEA 计算逻辑崩溃 (针对小规模数据集)
  set.seed(123)
  genelist <- genelist + rnorm(length(genelist), sd = 1e-9)
  genelist <- sort(genelist, decreasing = TRUE)

  # ==============================================================================
  # 2. 优化：单例加载逻辑 (使用标志变量 .current_sp_in_mem)
  # ==============================================================================
  base_dir <- "/share/nas1/wangsu/project/BMKSCTools_test/wangmy/GSEA_GSVA/bin/"
  rdata_file <- ifelse(species == "hs", "hs_msigdbr_data.RData", "mm_msigdbr_data.RData")
  rdata_path <- file.path(base_dir, rdata_file)
  
  # 检查标志变量是否匹配当前请求的物种
  need_load <- TRUE
  if (exists(".current_sp_in_mem", envir = .GlobalEnv)) {
    if (.GlobalEnv$.current_sp_in_mem == species) {
      need_load <- FALSE
    }
  }

  if (need_load) {
    message(">>> 正在加载本地数据库: ", rdata_file)
    if (!file.exists(rdata_path)) stop(paste0("Database not found: ", rdata_path))
    
    # 加载到全局环境
    load(rdata_path, envir = .GlobalEnv)
    
    # 加载后设置标志变量
    assign(".current_sp_in_mem", species, envir = .GlobalEnv)
  } else {
    message(">>> 检测到 ", species, " 数据库已在内存中，直接使用。")
  }

  # ==============================================================================
  # 3. 构建 TERM2GENE 映射表
  # ==============================================================================
  # 使用 .GlobalEnv$ 确保调用的是全局环境中的 msigdbr_data
  term2gene <- .GlobalEnv$msigdbr_data %>% 
    dplyr::filter(gs_name %in% target_ids) %>% 
    dplyr::select(term = gs_name, gene = gene_symbol) %>%
    dplyr::distinct()
  
  # 格式化通路名称以适配绘图
  term2gene$term <- sapply(term2gene$term, format_pathway_name)
  
  return(list(genelist = genelist, term2gene = term2gene))
}

# ==============================================================================
# 3. 绘图
# ==============================================================================
#' 动态配色并绘制 GSEA 图
#' @param processed_data 包含 genelist 和 term2gene 的列表
#' @param target_ids 目标通路 ID 向量
#' @param output_prefix 输出文件路径前缀
generate_dynamic_gsea_plots <- function(processed_data, target_ids, output_prefix) {
  
  # 1. 运行 GSEA 分析
  gsea_plot_obj <- GSEA(processed_data$genelist, 
                        TERM2GENE = processed_data$term2gene, 
                        pvalueCutoff = 1, 
                        minGSSize = 1, 
                        maxGSSize = 10000)
  # 2. 根据通路数量取出对应数量的配色
  # 使用 colorRampPalette 生成平滑渐变色，或使用内置调色盘
  n_ids <- length(target_ids)
  # 这里使用经典的 Set1 配色并扩展到 n_ids 个
  dynamic_colors <- colorRampPalette(c("#D62728", "#1F77B4", "#2CA02C", "#FF7F0E", "#9467BD"))(n_ids)
  
  # 3. 绘制 gseaNb 组合图
  p_nb <- gseaNb(object = gsea_plot_obj,
                 geneSetID = target_ids,
                 subPlot = 3,
                 curveCol = dynamic_colors, # 传入动态生成的颜色
                 htHeight = 0.5,
                 addPval = TRUE,
                 termWidth = 40
                 )
  
  ggsave(p_nb, filename = paste0(output_prefix, "_gsea.pdf"), width = 12, height = 10)
  ggsave(p_nb, filename = paste0(output_prefix, "_gsea.png"), width = 12, height = 10)
  
  # 4. 绘制山峦图 
  # 提取 GSEA 结果并筛选目标通路
res <- gsea_plot_obj@result %>% 
  dplyr::filter(ID %in% target_ids)

gene_list <- gsea_plot_obj@geneList

# ======================
# ridge 数据
# ======================
ridge_df <- res %>%
  dplyr::select(ID, core_enrichment, pvalue) %>%
  tidyr::separate_rows(core_enrichment, sep = "/") %>%
  dplyr::rename(gene = core_enrichment) %>%
  dplyr::mutate(
    log2FC = gene_list[gene]
  ) %>%
  dplyr::filter(!is.na(log2FC))

ridge_df$ID <- factor(ridge_df$ID, levels = rev(target_ids))

# ======================
# NES 数据
# ======================
nes_df <- res %>%
  dplyr::select(ID, NES) %>%
  dplyr::mutate(
    ID = factor(ID, levels = rev(target_ids))
  )

# ======================
# 坐标范围
# ======================
x_min <- min(ridge_df$log2FC, na.rm = TRUE)
x_max <- max(ridge_df$log2FC, na.rm = TRUE)
x_range <- x_max - x_min

# 👉 左边界
x_left <- x_min - 0.4 * x_range

# 👉 NES 位置（靠左边框）
x_nes <- x_left + 0.05 * x_range

# ======================
# 绘图
# ======================
p_ridge <- ggplot(ridge_df, aes(x = log2FC, y = ID, fill = pvalue)) +

  ggridges::geom_density_ridges(
    alpha = 0.9,
    color = "white",
    scale = 0.9
  ) +

  scale_fill_gradient(
    low = "#D73027",
    high = "#4575B4",
    name = "p-value"
  ) +

  # NES 点
  geom_point(
    data = nes_df,
    aes(
      x = x_nes,
      y = ID,
      color = NES,
      size = abs(NES)
    ),
    inherit.aes = FALSE
  ) +

  scale_color_gradient2(
    low = "#4575B4",
    mid = "white",
    high = "#D73027",
    midpoint = 0,
    name = "NES"
  ) +

  scale_size(range = c(2, 5), guide = "none") +

  # ✅ 关键：给边界留一点空间（防止被外框切掉）
  scale_x_continuous(
    limits = c(x_left, x_max),
    expand = expansion(mult = c(0.02, 0.02))  # ⭐ 关键修复
  ) +

  labs(
    x = "log2FoldChange",
    y = NULL
  ) +

  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10, color = "black"),
    legend.position = "right",

    # ❌ 全部去掉网格（关键）
    panel.grid = element_blank(),

    # ✅ 外框
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    ),

    plot.margin = margin(5.5, 20, 5.5, 40)
  )

# ======================
# 显示
# ======================
print(p_ridge)

# ======================
# 保存
# ======================
ggsave(
  filename = paste0(output_prefix, "_ridge.pdf"),
  plot = p_ridge,
  width = 10,
  height = 8
)

ggsave(
  filename = paste0(output_prefix, "_ridge.png"),
  plot = p_ridge,
  width = 10,
  height = 8,
  dpi = 300
)

  return(list(gsea_obj = gsea_plot_obj, nb_plot = p_nb, ridge_plot = p_ridge))
}


