library(tidyverse)
library(Seurat)
library(argparse)

parser = ArgumentParser()
parser$add_argument("--sample1", help="sample1基因占比文件路径")
parser$add_argument("--sample2", help="sample2基因占比文件路径")
parser$add_argument("--sample1_alias", help=" ")
parser$add_argument("--sample2_alias", help=" ")
parser$add_argument("--sample1_name", help=" ")
parser$add_argument("--sample2_name", help=" ")
parser$add_argument("--outdir", help=" ")
parser$add_argument("--compare", help=" ")
parser$add_argument("--RDSor10X", help="输入文件格式")
parser$add_argument("--path", help="输入文件路径")

args <- parser$parse_args()
str(args)

sample1=args$sample1
sample2=args$sample2
sample1_alias=args$sample1_alias
sample2_alias=args$sample2_alias
sample1_name=args$sample1_name
sample2_name=args$sample2_name
outdir=args$outdir
compare=args$compare
RDSor10X=args$RDSor10X
path=args$path

if(!dir.exists(outdir)){
	dir.create(outdir, recursive = T)
}

x <- read.csv(sample1,sep='\t')
y <- read.csv(sample2,sep='\t')


ob.list <- list()
samples<-strsplit(compare,'vs')[[1]]


if(RDSor10X == "RDS"){
numsap=1
for (each in samples){
	ob <- readRDS(str_glue("{path}/{each}/QC/{each}_QC.rds"))
	ob <- RenameCells(ob, new.names = paste0(colnames(ob),'-',numsap))
	ob$stim <-each
	ob <- NormalizeData(ob)
	ob <- FindVariableFeatures(ob,  selection.method = "vst",nfeatures = 2000)
	numsap=numsap+1
	ob.list[[each]] <- ob
}
}else if(RDSor10X == "10X"){
numsap=1
for (each in samples){
	pbmc.data <- Read10X(data.dir = paste0(path,'/',each) ,gene.column = 1)
	colnames(pbmc.data) <- paste0(colnames(pbmc.data),'-',numsap)
	ob <- CreateSeuratObject(counts =pbmc.data ,project =each,min.cells = 0 , min.features = 0 )
	ob$stim <-each
	ob <- NormalizeData(ob)
	ob <- FindVariableFeatures(ob,  selection.method = "vst",nfeatures = 2000)
	numsap=numsap+1
	ob.list[[each]] <- ob
}
}else{
    print("请输入RDS或10X")
}


anchors <- FindIntegrationAnchors(object.list = ob.list, dims = 1:20, anchor.features = 2000)
combined <- IntegrateData(anchorset = anchors, dims = 1:20)


DefaultAssay(combined) <- "RNA"
markers_com <- FindMarkers(combined, ident.1 = samples[1], group.by = 'stim', ident.2 = samples[2],logfc.threshold=0,min.pct=0)

save.image(file=paste0(outdir,'/load.RData'))

aa <- full_join(x, y, by = "gene")
aa$cell_ratio_diff <- (aa$cell_ratio.x-aa$cell_ratio.y)/(aa$cell_ratio.x+aa$cell_ratio.y)*2
# aa$cell_average_diff <- markers_com$avg_log2FC

# aa$p_val_adj <- markers_com$p_val_adj
aa <- aa %>% full_join(markers_com %>% rownames_to_column(var = "gene") %>% select(gene, avg_log2FC, p_val_adj), by = "gene") %>% 
	rename("cell_average_diff" = "avg_log2FC")

aa <- aa[,c(1,2,4,6,3,5,7,8)]

colnames(aa) <- c('Gene_name',paste0('cell_ratio_',sample1_alias),paste0('cell_ratio_',sample2_alias),'cell_ratio_diff',paste0('cell_average_',sample1_alias),paste0('cell_average_',sample2_alias),'cell_average_diff','p_val_adj')

write.table(aa,file = paste0(outdir,paste0('/',sample1_name,'vs',sample2_name,'.xls')),sep='\t',quote = F,row.names =F)

cell_ratio_up <- subset(aa , cell_ratio_diff  > 0.4  )

cell_ratio_down <- subset(aa , cell_ratio_diff  < -0.4  )

cell_average_up <- subset(aa , cell_average_diff  > 0.25 & p_val_adj < 0.05 )

cell_average_down <- subset(aa , cell_average_diff  < -0.25 & p_val_adj < 0.05 )

write.table(cell_ratio_up,file = paste0(outdir,paste0('/',sample1_name,'vs',sample2_name,'_ratio_up.xls')),sep='\t',quote = F,row.names =F) 
write.table(cell_ratio_down,file = paste0(outdir,paste0('/',sample1_name,'vs',sample2_name,'_ratio_down.xls')),sep='\t',quote = F,row.names =F)
write.table(cell_average_up,file = paste0(outdir,paste0('/',sample1_name,'vs',sample2_name,'_average_up.xls')),sep='\t',quote = F,row.names =F)
write.table(cell_average_down,file = paste0(outdir,paste0('/',sample1_name,'vs',sample2_name,'_average_down.xls')),sep='\t',quote = F,row.names =F)
