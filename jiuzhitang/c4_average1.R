library(Seurat)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(ConfigParser)
library(reticulate)
library(KernSmooth)
library(fields)


library(argparse)
parser = ArgumentParser()
parser$add_argument("--RDSor10X", help="输入文件格式")
parser$add_argument("--cellranger_path", help="path/to/outs/*_bc_matrix/:tsv.gz or tsv files")
parser$add_argument("--sample", help="the sample name of project")
parser$add_argument("--outdir", help="outdir of project")
parser$add_argument("--species", help="the species name of project")
args <- parser$parse_args()
str(args)

RDSor10X=args$RDSor10X
cellranger_path=args$cellranger_path
sample=args$sample
outdir=args$outdir
species=args$species

if (!dir.exists(paste0(outdir))){
  dir.create(paste0(outdir))
}


if(RDSor10X == "RDS"){
pbmc = readRDS(cellranger_path)
}else if(RDSor10X == "10X"){
pbmc.data <- Read10X(data.dir = cellranger_path ,gene.column = 1)
pbmc <- CreateSeuratObject(counts = pbmc.data, min.cells  = 0, min.features= 0,project = sample)
}

pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize",scale.factor = 10000)

a <- pbmc@assays$RNA@data

###############################################

row_ave <- rowMeans(a, na.rm=TRUE)

################################################

# # 获取 pbmc@assays$RNA@counts 矩阵
# counts_matrix <- pbmc@assays$RNA@counts

# # 计算每行中的非零元素个数
# non_zero_counts <- rowSums(counts_matrix != 0)

# # 计算每行中非零元素占总列数的比例
# rows <- non_zero_counts / ncol(counts_matrix)
##############################################
rows = diff(t(pbmc@assays$RNA@counts)@p)/length(pbmc@meta.data$orig.ident)

tj <- data.frame(rownames(pbmc@assays$RNA@counts),rows,row_ave)
colnames(tj) <- c("gene","cell_ratio","cell_average")


write.table(tj,file = paste0(outdir,'/',sample,'样本基因占比统计表1','.xls'),col.names =T,row.names= F ,sep='\t',quote = F)


