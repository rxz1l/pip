# load packages -----
library(Seurat)
library(optparse)
library(tidyverse)
# library(permute)
# library(data.table)
# library(SCopeLoomR)
# library(dataNIC)

# 传参 -----
option_list = list(
  make_option(c("-i", "--input_path"), type = "character", default = FALSE,
              action = "store", help = "input path to data"
  ),
  make_option(c("-p", "--percentage"), type = "double", default = FALSE,
              action = "store", help = "The proportion of data extracted for each cell type"
  ),
  make_option(c("-o", "--out_dir"), type = "character", default = FALSE,
              action = "store", help = "The output directory"
  )
)
opt = parse_args(OptionParser(option_list = option_list))

print(opt)

input_path = opt$input_path
percentage = opt$percentage
out_dir = opt$out_dir

# load data -----
data = readRDS(input_path)

# 过滤 -----
# data1 = subset(data, nFeature_RNA > 200 & nFeature_RNA < 4000)
# data1 = PercentageFeatureSet(data1, "^MT-", col.name = "percent_mito")
# data1 = subset(data1, percent_mito < 15)
# data2 = data1
# data2 = data2[rowSums(data2@assays$RNA@counts>0)>=1000,]
# data3 = data1
# data3 = data3[rowSums(data3@assays$RNA@counts>0)>=10,]

# sce_test <- subset(x = uterus, downsample = 2000)

cell_ids <- colnames(data)
cell_source <- factor(data$celltype)
cell_ids_list <- lapply(split(cell_ids, cell_source), function(x){sample(x, length(x)*percentage)})
cell_id <- unlist(cell_ids_list)
data4= data[,cell_id]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(t(as.matrix(data4@assays$RNA@counts)),file = str_glue("{out_dir}/scenic.data.csv"))
saveRDS(data4, file = str_glue("{out_dir}/scenic.RDS"))
