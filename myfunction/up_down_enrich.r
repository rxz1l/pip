library(clusterProfiler)
library(org.Mm.eg.db) 
library(ggplot2)
library(dplyr)
library(stringr)
library(readr)

plot_enrichment_v3 <- function(df, 
                               method = "GO", 
                               ont = "BP", 
                               org_db = org.Mm.eg.db, 
                               topN = 5,
                               title = "Enrichment Analysis",
                               out_dir = NULL) {
  
  # 设置 KEGG 下载方式，缓解连接问题
  options(clusterProfiler.download.method = "auto") 
  
  # 1. 筛选基因
  up_symbols <- df %>% filter(log2FC > 0 & Qvalue < 0.05) %>% pull(symbol)
  down_symbols <- df %>% filter(log2FC < 0 & Qvalue < 0.05) %>% pull(symbol)
  
  # 2. 富集逻辑
  run_analysis <- function(genes, is_up = TRUE) {
    if (length(genes) == 0) return(NULL)
    
    # 尝试运行富集，捕获 KEGG 网络错误
    res <- tryCatch({
      if (method == "GO") {
        enrichGO(gene = genes, OrgDb = org_db, keyType = "SYMBOL", 
                 ont = ont, pAdjustMethod = "BH", pvalueCutoff = 0.05)
      } else {
        species <- ifelse(grepl("Mm", class(org_db)[1]), "mmu", "hsa")
        ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org_db)
        enrichKEGG(gene = ids$ENTREZID, organism = species, pvalueCutoff = 0.05)
      }
    }, error = function(e) {
      message(paste0("Error in ", method, " analysis: ", e$message))
      return(NULL)
    })
    
    if (is.null(res) || nrow(res) == 0) return(NULL)
    
    # --- 核心改进：为 Description 增加唯一后缀防止重叠 ---
    status_tag <- ifelse(is_up, " (Up)", " (Down)")
    
    as.data.frame(res) %>%
      slice_head(n = topN) %>%
      mutate(Status = ifelse(is_up, "Upregulated", "Downregulated"),
             logP = -log10(p.adjust),
             # 关键：创建唯一 ID 用于 y 轴
             UniqueID = paste0(Description, status_tag))
  }
  
  up_res <- run_analysis(up_symbols, is_up = TRUE)
  down_res <- run_analysis(down_symbols, is_up = FALSE)
  plot_df <- rbind(up_res, down_res)
  
  if (is.null(plot_df) || nrow(plot_df) == 0) {
    message("No results to plot.")
    return(NULL)
  }
  
  # --- 3. 绘图 (使用 UniqueID 排序，但显示 Description) ---
  p <- ggplot(plot_df, aes(x = logP, y = reorder(UniqueID, logP), fill = Status)) +
    geom_col(color = "black", width = 0.7, size = 0.3) +
    scale_fill_manual(values = c("Upregulated" = "#F8766D", "Downregulated" = "#71CAFF")) +
    # 核心改进：映射回原始名称，并自动换行
    scale_y_discrete(labels = function(x) {
      original_names <- gsub(" \\(Up\\)| \\(Down\\)", "", x)
      str_wrap(original_names, width = 50)
    }) +
    theme_classic() +
    labs(title = title, x = expression(-log[10](Q~value)), y = NULL) +
    theme(axis.text.y = element_text(size = 9, color = "black"),
          legend.position = "bottom", 
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # 4. 保存
  if (!is.null(out_dir)) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    file_prefix <- paste0(method, if(method=="GO") paste0("_", ont) else "", "_", gsub(" ", "_", title))
    ggsave(filename = file.path(out_dir, paste0(file_prefix, ".pdf")), plot = p, width = 8, height = 6)
    ggsave(filename = file.path(out_dir, paste0(file_prefix, ".png")), plot = p, width = 8, height = 6)
    write_csv(plot_df, file.path(out_dir, paste0(file_prefix, "_data.csv")))
  }
  
  return(list(plot = p, data = plot_df))
}
