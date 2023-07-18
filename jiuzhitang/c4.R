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
parser$add_argument("--species", help="the species name of project")
parser$add_argument("--mark", help="mark gene + gene type")
parser$add_argument("--outdir", help="outdir of project")
args <- parser$parse_args()
str(args)

RDSor10X=args$RDSor10X
cellranger_path=args$cellranger_path
sample=args$sample
species=args$species
mark=args$mark
outdir=args$outdir

if (!dir.exists(paste0(outdir))){
  dir.create(paste0(outdir))
}


# #pbmc.data<-Read10X(data.dir = cellranger_path)
# pbmc.data <- Read10X(data.dir = cellranger_path ,gene.column = 1)
# pbmc <- CreateSeuratObject(counts = pbmc.data, min.cells  = 0, min.features= 0,project = sample)
if(RDSor10X == "RDS"){
pbmc = readRDS(cellranger_path)
}else if(RDSor10X == "10X"){
pbmc.data <- Read10X(data.dir = cellranger_path ,gene.column = 1)
pbmc <- CreateSeuratObject(counts = pbmc.data, min.cells  = 0, min.features= 0,project = sample)
}
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize",scale.factor = 10000)

a <- pbmc@assays$RNA@data
#write.table(a,file = paste0('/public/work/Project/Single_cell/Project_c4/c4-4_jiuzhitang_20220624/analyse/seurat/b.csv'),col.names =F,sep=',')
b <- colnames(a)
c <- c('gene_name')
d <- c(c,b)
f <- as.data.frame(d)
g <- t(f)
#write.table(g,file = paste0('/public/work/Project/Single_cell/Project_c4/c4-4_jiuzhitang_20220624/analyse/seurat/a.csv'),col.names =F,row.names= F ,sep=',')

row_n <- data.frame(rownames(pbmc@assays$RNA@counts))
row_add <- cbind(row_n ,pbmc@assays$RNA@counts)
g <- as.data.frame(g)
names(row_add) <- names(g)
a_add <- rbind(g ,row_add)

write.table(a_add,file = paste0(outdir,'/',sample,'.csv'),col.names =F,row.names= F ,sep=',',quote = F)

# 获取 pbmc@assays$RNA@counts 矩阵
# counts_matrix <- pbmc@assays$RNA@counts

# # 计算每行中的非零元素个数
# non_zero_counts <- rowSums(counts_matrix != 0)

# # 计算每行中非零元素占总列数的比例
# rows <- non_zero_counts / ncol(counts_matrix)

rows = diff(t(pbmc@assays$RNA@counts)@p)/length(pbmc@meta.data$orig.ident)

# names(rows) = rownames(pbmc)

tj <- data.frame(rownames(pbmc@assays$RNA@counts),rows)
colnames(tj) <- c("gene","cell_ratio")
tj <- arrange(tj, -cell_ratio)
rownames(tj) = NULL
if (species=='Human'){
 alias <- read.csv('/public/work/Pipline/Single_RNA/standard_analysis/C4/tongji_ratio/data/human_alias.xls',sep='\t')
 protein <- read.csv('/public/work/Pipline/Single_RNA/standard_analysis/C4/tongji_ratio/data/human_protein.xls',sep='\t')
}else{
 alias <- read.csv('/public/work/Pipline/Single_RNA/standard_analysis/C4/tongji_ratio/data/mouse_alias.xls',sep='\t')
 protein <- read.csv('/public/work/Pipline/Single_RNA/standard_analysis/C4/tongji_ratio/data/mouse_protein.xls',sep='\t')
}


tj_alias <- left_join(tj, alias, by="gene")
tj_alias_protein <- left_join(tj_alias, protein, by="gene")

write.table(tj_alias_protein,file = paste0(outdir,'/',sample,'样本基因占比统计表','.xls'),col.names =T,row.names= F ,sep='\t',quote = F)


library(dplyr)
#tj <- arrange(tj, -cell_ratio)
tj$ord <- rownames(tj)
tj$ord_ratio <- as.numeric(rownames(tj))/nrow(tj)

type <- read.csv(mark,sep='\t')
tj_type <- inner_join(tj, type, by="gene")

library(scales)
library(ggrepel)

p <- ggplot(data = tj)+geom_point(aes(x=ord_ratio,y=cell_ratio),tj,color='grey',size=0.6)+labs(x="Rank of gene expression percentage(%)",y="Percent of cell",title="Percent of cell with mRNA detected")+scale_y_continuous(labels=percent) +scale_x_continuous(labels=percent,expand = c(0, 0))+geom_point(aes(x=ord_ratio,y=cell_ratio,color=type),tj_type,size=3)+theme_classic()+theme(plot.title = element_text(hjust = 0.5,size=20,margin = margin(0.5,0,0.5,0,'cm')),axis.title.x =element_text(size=18,margin = margin(0.8,1,0,1,'cm')),axis.title.y =element_text(margin = margin(0,0.5,0,0,'cm'),size=18),axis.text=element_text(size=12))+scale_color_manual(values=c("Negative marker" = "#800000","Positive marker" = "#007799"))+geom_label_repel(aes(ord_ratio,cell_ratio , fill=type, label=gene),tj_type,force = 1,nudge_x=0.15,box.padding = unit(0.5, "lines"),show.legend = FALSE)

ggsave(paste0(outdir,'/',sample,'样本基因占比统计图','.pdf'),p,height=9,width =9)
ggsave(paste0(outdir,'/',sample,'样本基因占比统计图','.png'),p,type = 'cairo-png',height=9,width =9)


# save.image(file=paste0(outdir,'/',sample,'.RData'))

