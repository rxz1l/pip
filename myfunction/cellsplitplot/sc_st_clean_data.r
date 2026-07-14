#!/usr/bin/env Rscript

# 检查并自动安装缺少的依赖包
required_packages <- c("argparse", "Seurat", "tidyverse")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, repos = "https://cloud.r-project.org")

library(argparse)
library(Seurat)
library(tidyverse)

# 1. 设置命令行参数解析器 ----
parser <- ArgumentParser(description = "处理空间转录组(ST)多方法细胞类型注释结果的集成脚本")

parser$add_argument("-i", "--rds_path", type = "character", required = TRUE,
                    help = "输入的 Seurat RDS 文件路径 (例如: .../GL261.single_seruat.Rds)")

parser$add_argument("-c", "--csv_dir", type = "character", required = TRUE,
                    help = "包含 MIA, Anchor, Spotlight, RCTD 结果原始文件的根目录 (例如: .../result/sc_st/)")

parser$add_argument("-o", "--out_dir", type = "character", required = TRUE,
                    help = "输出的 CSV 结果保存目录 (例如: .../script/cellsplit/barcode2cluster/)")

args <- parser$parse_args()

# 创建输出目录（如果不存在的话）
if (!dir.exists(args$out_dir)) {
    dir.create(args$out_dir, recursive = TRUE)
    message("创建输出目录: ", args$out_dir)
}

# 2. 读取 ST 数据 ----
message("正在读取 Seurat RDS 文件...")
st <- readRDS(args$rds_path)


# 3. 处理 MIA 数据 ----
mia_file <- file.path(args$csv_dir, "BMK2_MIA/MIA_cluster_cellType.xls")
if (file.exists(mia_file)) {
    message("正在处理 MIA 数据...")
    MIA <- read.table(file = mia_file, header = TRUE)
    
    MIA_result <- st@meta.data %>% 
        select(seurat_clusters) %>%                                
        rownames_to_column(var = "Barcode") %>%                                              
        inner_join(MIA %>% mutate(cluster = as.factor(cluster)), by = c("seurat_clusters" = "cluster")) %>%                                                                           
        select(Barcode, cellType) %>% 
        rename(Cluster = cellType)                             
    
    write.csv(MIA_result, file = file.path(args$out_dir, "MIA_result.csv"), row.names = FALSE)
} else {
    warning("未找到 MIA 配置文件: ", mia_file)
}


# 4. 处理 Anchor 数据 ----
anchor_file <- file.path(args$csv_dir, "BMK3_Anchor/Anchor_barcode_cellType.xls")
if (file.exists(anchor_file)) {
    message("正在处理 Anchor 数据...")
    Anchor <- read.table(anchor_file, header = TRUE)
    
    Anchor_result <- Anchor %>% 
        rename(Barcode = barcode, Cluster = cellType)                        
    
    write.csv(Anchor_result, file = file.path(args$out_dir, "Anchor.csv"), row.names = FALSE)
} else {
    warning("未找到 Anchor 配置文件: ", anchor_file)
}


# 5. 处理 Spotlight 数据 ----
spotlight_file <- file.path(args$csv_dir, "BMK4_Spotlight/Spotlight_barcode_cellType.xls")
if (file.exists(spotlight_file)) {
    message("正在处理 Spotlight 数据...")
    Spotlight <- read.table(spotlight_file, header = TRUE)
    
    Spotlight_result <- Spotlight %>% 
        rename(Barcode = barcode, Cluster = cellType)                  
    
    write.csv(Spotlight_result, file = file.path(args$out_dir, "Spotlight.csv"), row.names = FALSE)
} else {
    warning("未找到 Spotlight 配置文件: ", spotlight_file)
}


# 6. 处理 RCTD 数据 ----
rctd_file <- file.path(args$csv_dir, "BMK5_RCTD/RCTD.xls")
if (file.exists(rctd_file)) {
    message("正在处理 RCTD 数据...")
    RCTD <- read.table(rctd_file, header = TRUE)                                                     
    
    RCTD_result <- RCTD %>% 
        select(barcode, RCTD_celltype) %>%                                       
        rename(Barcode = barcode, Cluster = RCTD_celltype)                                   
    
    write.csv(RCTD_result, file = file.path(args$out_dir, "RCTD.csv"), row.names = FALSE) 
} else {
    warning("未找到 RCTD 配置文件: ", rctd_file)
}

message("所有任务运行完成！")
