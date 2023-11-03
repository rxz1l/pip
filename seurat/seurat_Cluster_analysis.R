### Modified on April 6th, 2023
### Author: Mengyao Ma

library(Seurat)
library(dplyr)
library(stringr)
library(argparse)
library(reshape2)
library(tidyr)
library(ggplot2)
library(scales)
library(cowplot)
library(reticulate)
library(tidyverse)
library(psych)
library(pheatmap)
# library(COSG)
# use_python("/public/work/Personal/liujinxu/software/Anaconda/Anaconda3-2021.11/envs/report_python/bin/python")


parser = ArgumentParser()
parser$add_argument("--path", help="/path/to/samples gene bar")
parser$add_argument("--compare", help="samplevssamplevssample,using vs as the split")
parser$add_argument("--species", help="GRCh38 or mm10")
parser$add_argument("--outdir", help="outdir of project")
parser$add_argument("--Nfeatures",help="nfeatures of FindVariableFeatures()", default='2000', type = "integer")
parser$add_argument("--resolution", help="resolution for cluster",default='0.6',type = "double")
parser$add_argument("--min_cells", help="genes expressed in min cells",default='3', type = "integer")
parser$add_argument("--dim_use", help="dim_use for cluster",default='30', type = "integer")
parser$add_argument('--pc',help='pc usage', default = "50", type = "integer")
#parser$add_argument('--only.pos',help='if only return positive markers in FindAllmarkers', default = "F")
parser$add_argument("--dimension_method", help=' Which dimensionality reduction to use (umap tsne) ',default="umap")
parser$add_argument("--prefix", help="prefix of results")
parser$add_argument("--singler", help="performing SingleR:yes or no", default="no")
# parser$add_argument("--ref_gene", help='gene matrix')
# parser$add_argument("--cluster", help=' if you have other cluster method give me')
args <- parser$parse_args()
str(args)

path=args$path
compare=args$compare
species = args$species
outdir=args$outdir
resolution=args$resolution
min_cells=args$min_cells
Nfeatures = args$Nfeatures
dim.use=args$dim_use
pc.usage=args$pc
dimension_method=args$dimension_method
#only.pos=args$only.pos
prefix=args$prefix
singler=args$singler

# ref_gene=args$ref_gene
# cluster_method = args$cluster

#################################################### step1:create directory

if (!dir.exists(outdir)){
  dir.create(outdir)
}
if (!dir.exists(paste0(outdir,'/Anchors'))){
  dir.create(paste0(outdir,'/Anchors'))
}
if (!dir.exists(paste0(outdir,'/Dimension_reduction'))){
  dir.create(paste0(outdir,'/Dimension_reduction'))
}

if (!dir.exists(paste0(outdir,'/Marker'))){
  dir.create(paste0(outdir,'/Marker'))
}

if (!dir.exists(paste0(outdir,'/CellsRatio'))){
  dir.create(paste0(outdir,'/CellsRatio'))
}

#################################################### step2:read-in data
ob.list <- list()
samples <- strsplit(compare,'vs')[[1]] # sample names

# 多个样品整合
if (length(samples) > 1) {
  numsap=1
  for (each in samples){
    ob <- readRDS(paste0(path,'/',each,'/QC/',each,'_QC.rds'))
    # if(grep("-1", colnames(ob@assays$RNA@counts)[1])){ # 10x file
    if(length(grep('-1',colnames(ob@assays$RNA@counts)[1])) == 1){
      colnames(ob@assays$RNA@counts) <- str_replace_all(colnames(ob@assays$RNA@counts), '-1',paste0('-',numsap))
    }else{
      colnames(ob@assays$RNA@counts) <- paste0(colnames(ob@assays$RNA@counts),'-',numsap) # C4 file
    }
    ob <- CreateSeuratObject(counts =ob@assays$RNA@counts,project =each,min.cells = min_cells)
    ob$batch <- each
    ob <- ob %>%
      NormalizeData() %>%
      FindVariableFeatures(selection.method = "vst",nfeatures = Nfeatures)
    numsap=numsap+1
    ob.list[[each]] <- ob
  }
  # Find anchor by Seurat or # SCT harmony
  anchors <- FindIntegrationAnchors(object.list = ob.list, dims = 1:20, anchor.features = Nfeatures)
  object.combined <- IntegrateData(anchorset = anchors, dims = 1:20)

  # store anchor genes
  object.combined[['RNA']] %>% rownames() %>%
    as_tibble() %>% rename("Gene_name" = "value") %>%
    write.table(file= paste0(outdir,'/Anchors/',"anchored_gene.csv"),sep=',',quote = F,row.names =F)
  DefaultAssay(object.combined) <- "integrated"

} else if (length(samples)==1){ # 单一样品
  object.combined <- readRDS(paste0(path,'/',samples,'/QC/',samples,'_QC.rds'))
}

object.combined <- object.combined %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = pc.usage, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:dim.use, n.components = 3L, 
          umap.method = "umap-learn",metric = "correlation") %>%
  RunTSNE(reduction = "pca", dims = 1:dim.use) %>%
  FindNeighbors(reduction = "pca", dims = 1:dim.use) %>%
  FindClusters(resolution = resolution, algorithm = 1)

# store high variable genes
object.combined %>% VariableFeatures() %>% as_tibble() %>%
  rename("Gene_name" = "value") %>%
  write.table(file = paste0(outdir,'/Dimension_reduction/',prefix,"_high_varience_gene.csv"),sep=',',quote = F,row.names =F)

# object.combined %>% saveRDS(paste0(outdir,'/before_UMAPplot.rds'))
# # Visualization: PCA
# top10 <- object.combined %>% VariableFeatures() %>% head(10)
# plot_feature <- VariableFeaturePlot(object.combined) # 基因在所有细胞中的平均表达量，红色为高变基因。
# plot_feature2 <- LabelPoints(plot = plot_feature, points = top10, repel = TRUE)
# # ggsave(paste0(outdir, "/Dimension_reduction/", compare, "_Dispersion.pdf"), plot_all, width = 12, height = 9)
# pdf(paste0(outdir,'/Dimension_reduction/',compare,'_PCA_Dispersion.pdf'),width = 12, height = 9)
# print(CombinePlots(plots = list(plot_feature, plot_feature2),legend="bottom"))
# dev.off()
#
# png(paste0(outdir,'/Dimension_reduction/',compare,'_PCA_Dispersion.png'),type = 'cairo-png', width = 12, height = 9)
# print(CombinePlots(plots = list(plot_feature, plot_feature2),legend="bottom"))
# dev.off()
#
# pdf(paste0(outdir,'/Dimension_reduction/',compare,'_PCA_Elbow.pdf'),width = 12, height = 9)
# print(ElbowPlot(object.combined))
# dev.off()
#
# png(paste0(outdir,'/Dimension_reduction/',compare,'_PCA_Elbow.png'),type = 'cairo-png', width = 12, height = 9)
# print(ElbowPlot(object.combined))
# dev.off()
#
# pdf(paste0(outdir,'/Dimension_reduction/',compare,'_PCA_Heatmap.pdf'), width = 12, height = 9)
# DimHeatmap(object.combined, dims = 1:dim.use, cells = 500, balanced = TRUE)
# dev.off()
#
# png(paste0(outdir,'/Dimension_reduction/',compare,'_PCA_Heatmap.png'),type="cairo-png", width = 12, height = 9)
# DimHeatmap(object.combined, dims = 1:dim.use, cells = 500, balanced = TRUE)
# dev.off()


# Visualization: UMAP
object.combined@reductions$umap@cell.embeddings %>% as.data.frame() %>%
  rownames_to_column(var = 'Barcode') %>%
  write.table(file=paste0(outdir,'/Dimension_reduction/',prefix,'_UMAP.csv'),row.names=F,col.names=T,sep=',',quote=F)

p1 <- object.combined %>% DimPlot(reduction = "umap", group.by = "batch")
p2 <- object.combined %>% DimPlot(reduction = "umap", label = TRUE)
p_split_sample_UMAP <- object.combined %>% DimPlot(reduction = "umap", split.by = "batch",ncol = 2)
ggsave(paste0(outdir,'/Dimension_reduction/',prefix,'_sample_UMAP.pdf'),p1, width = 12, height = 9)
ggsave(paste0(outdir,'/Dimension_reduction/',prefix,'_sample_UMAP.png'),p1,type="cairo-png", width = 12, height = 9)
ggsave(paste0(outdir,'/Dimension_reduction/',prefix,'_clusters_UMAP.pdf'),p2, width = 12, height = 9)
ggsave(paste0(outdir,'/Dimension_reduction/',prefix,'_clusters_UMAP.png'),p2,type="cairo-png", width = 12, height = 9)
ggsave(paste0(outdir,'/Dimension_reduction/',prefix,'_sample_split_UMAP.pdf'),p_split_sample_UMAP,width = 22, height = 10)
ggsave(paste0(outdir,'/Dimension_reduction/',prefix,'_sample_split_UMAP.png'),p_split_sample_UMAP,type="cairo-png",width = 22, height = 10)

# Visualization: t-SNE
object.combined@reductions$tsne@cell.embeddings %>% as.data.frame() %>%
  rownames_to_column(var = 'Barcode') %>%
  write.table(file=paste0(outdir,'/Dimension_reduction/',prefix,'_tSNE.csv'),row.names=F, col.names=T,sep=',',quote=F)

p3 <- object.combined %>% DimPlot(reduction = "tsne", group.by = "batch")
p4 <- object.combined %>% DimPlot(reduction = "tsne", label = TRUE)
p_split_sample_TSNE <- object.combined %>% DimPlot(reduction = "tsne", split.by = "batch",ncol = 2)

ggsave(paste0(outdir,'/Dimension_reduction/',prefix,'_sample_TSNE.pdf'),p3)
ggsave(paste0(outdir,'/Dimension_reduction/',prefix,'_sample_TSNE.png'),p3,type="cairo-png", width = 12, height = 9)
ggsave(paste0(outdir,'/Dimension_reduction/',prefix,'_clusters_TSNE.pdf'),p4)
ggsave(paste0(outdir,'/Dimension_reduction/',prefix,'_clusters_TSNE.png'),p4,type="cairo-png", width = 12, height = 9)
ggsave(paste0(outdir,'/Dimension_reduction/',prefix,'_sample_split_TSNE.pdf'),p_split_sample_TSNE,width = 22, height = 10)
ggsave(paste0(outdir,'/Dimension_reduction/',prefix,'_sample_split_TSNE.png'),p_split_sample_TSNE,type="cairo-png",width = 22, height = 10)

object.combined %>% Idents() %>% as.data.frame() %>%
  rownames_to_column(var = "Barcode") %>% rename("Cluster" = ".") %>%
  write.table(file=paste0(outdir,'/Dimension_reduction/',prefix,'_cluster.csv'),row.names=T, col.names = T,quote=F,sep=',')


### Calculate cell ratios
### 计算每个sample中有多少个细胞（num 和%）在每个cluster中
### Calculate cell ratios
sample_clus_count <- table(object.combined$batch,object.combined$seurat_clusters)
sample_clus_count <- sample_clus_count %>% as.data.frame() %>% dplyr::rename(sample = Var1, cluster = Var2, num = Freq)
sample_clus_count %>% write.table(file=paste0(outdir,'/CellsRatio/',prefix,'_each_sample_celltype_abundance.csv'),sep=',',row.names = F, quote=F)

sample_clus_perc <- prop.table(table(object.combined$batch,object.combined$seurat_clusters),1)
sample_clus_perc <- sample_clus_perc %>% as.data.frame() %>% dplyr::rename(sample = Var1, cluster = Var2, ratio = Freq)
sample_clus_perc %>% write.table(file=paste0(outdir,'/CellsRatio/',prefix,'_each_sample_celltype_perc.csv'),sep=',',row.names = F, quote=F)

# For each cluster, calculate the ratios of cells in each sample
clus_sample_perc <- prop.table(table(object.combined$seurat_clusters, object.combined$batch),1)
clus_sample_perc <- clus_sample_perc %>% as.data.frame() %>% dplyr::rename(cluster = Var1, sample = Var2, ratio = Freq)
clus_sample_perc %>% write.table(file=paste0(outdir,'/CellsRatio/',prefix,'_each_cluster_sample_perc.csv'),sep=',',row.names = F, quote=F)

# 2 barplots
barplot1 <- clus_sample_perc %>%
  mutate(cluster = as.numeric(cluster)) %>% arrange(cluster) %>%
  ggplot(aes(x = factor(cluster), y = ratio, fill = sample)) +
  geom_bar(position = 'stack', stat = 'identity') +
  labs(x = "Cluster Name", y = "Cell Ratios") +
  theme(panel.background=element_rect(fill='transparent', color='black'),
        legend.key=element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(color="black")) +
  scale_y_continuous(expand=c(0.001,0.001)) +
  guides(fill = guide_legend(keywidth = 1, keyheight = 1,ncol=1,title = 'Sample')) +
  ggtitle("Percentage of cells in each group for different cell types")

ggsave(filename = paste0(outdir,'/CellsRatio/',prefix,'_each_cluster_sample_percent.pdf'), barplot1, dpi = 320,width = 12, height = 9)
ggsave(filename = paste0(outdir,'/CellsRatio/',prefix,'_each_cluster_sample_percent.png'), barplot1, dpi = 320, width = 12, height = 9, type="cairo-png")

#柱状图
barplot2 <- sample_clus_perc %>%
  mutate(cluster = as.numeric(cluster)) %>% arrange(cluster) %>%
  ggplot(aes(x = sample, y = ratio, fill = factor(cluster))) +
  geom_bar(position = 'stack', stat = 'identity') +
  labs(x = "Cluster Name", y = "Cell Ratios") +
  theme(panel.background=element_rect(fill='transparent', color='black'),
        legend.key=element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(color="black")) +
  scale_y_continuous(expand=c(0.001,0.001)) +
  guides(fill = guide_legend(keywidth = 1, keyheight = 1,ncol=1,title = 'Cluster')) +
  ggtitle("Percentage of cells in each cluster for different samples")

ggsave(filename = paste0(outdir,'/CellsRatio/',prefix,'_each_sample_cluster_percent.pdf'), barplot2, dpi = 320, width = 12, height = 9)
ggsave(filename = paste0(outdir,'/CellsRatio/',prefix,'_each_sample_cluster_percent.png'), barplot2, dpi = 320, type="cairo-png",  width = 12, height = 9)

# 导出cluster以及相关信息
if(dim(object.combined)[2] > 50){
  cluster_ID <- object.combined %>% Idents() %>% as.data.frame()
  cluster_cor <- object.combined %>% Embeddings(reduction = "umap") %>% as.data.frame()
  coor <- cluster_ID %>%
    cbind(cluster_cor, object.combined[["nCount_RNA"]], object.combined[['nFeature_RNA']]) %>%
    set_names(c("Cluster", "UMAP_1", "UMAP_2", "UMAP_3", "nUMI", "nGene")) %>%
    arrange(Cluster) %>%
    rownames_to_column(var = "Barcode") %>%
    add_count(Cluster, name = 'cellNum') %>%
    write.table(file=paste0(outdir,'/CellsRatio/',prefix,'_cluster_cell_stat.csv'),sep=',',row.names = F, quote=F)

  ### Find All Markers（cluster之间差异基因分析）
  ### 超过一个cluster
  cluster_num <- object.combined@meta.data$seurat_clusters %>% unique() %>% length()

  DefaultAssay(object.combined) <- "RNA"
  if (cluster_num>1){
    object.combined@misc$markers <- FindAllMarkers(object = object.combined, only.pos= T,min.pct = 0.25)
    combined.markers <- object.combined@misc$markers
    colnames(combined.markers)[7] <- 'Gene_name'
    combined.markers <- combined.markers[,c(7,6,2,1,5,3,4)]
    # combined.markers %>% write.table(file=paste0(outdir,'/DEG/','Cluster_DEG.txt'),quote=F,row.names=F,col.names=T,sep='\t')
    # combined.markers %>% write.table(file=paste0(outdir,'/DEG/','Cluster_DEG.csv'),quote=F,row.names=F,col.names=T,sep=',')

    # calculate the average expression of all the DEGs in each cluster and save the results into txt/csv.
    all_markers <- combined.markers$Gene_name
    gene_cell_exp <- object.combined %>%
      Seurat::AverageExpression(assays = "RNA", features = all_markers, group.by = "seurat_clusters")
    gene_cell_exp <- gene_cell_exp$RNA %>% as.data.frame() %>% rownames_to_column(var = "Gene_name")
    combined.markers.avg.df <- combined.markers %>% left_join(gene_cell_exp, by = "Gene_name")

    combined.markers.avg.df %>% write.table(paste0(outdir, '/Marker/avg_exp_allmarkers.csv'), quote = F, row.names = F, col.names = T, sep = ",")
    combined.markers.avg.df %>% write.table(paste0(outdir, '/Marker/avg_exp_allmarkers.txt'), quote = F, row.names = F, col.names = T, sep = "\t")

    # select top 10 markers based on the avg_log2FC values
    # save all the related statistics of the top 10 markers for each cluster into txt/csv
    combined.markers.avg.df %>% group_by(cluster) %>% top_n(10, avg_log2FC) %>%
      write.table(file=paste0(outdir,'/Marker/','top10_marker_each_cluster.txt'),quote=F,row.names=F,col.names=T,sep='\t')
    combined.markers.avg.df %>% group_by(cluster) %>% top_n(10, avg_log2FC) %>%
      write.table(file=paste0(outdir,'/Marker/','top10_marker_each_cluster.csv'),quote=F,row.names=F,col.names=T,sep=',')

    ######----------- Extract the top 4 markers in each cluster and generate several plots for better
    top_markers <- combined.markers %>% group_by(cluster) %>% top_n(4, avg_log2FC)
    top_gene <- top_markers$Gene_name %>% unique()
    clus<-unique(top_markers$cluster)

    # violin and tsne graph
    for (each in clus){
      top_clus <- top_markers %>% filter(cluster == each)
      top_gene_name <- top_clus$Gene_name
      top2_features<-(top_clus$Gene_name)[1:2]

      plot1_marker <- FeaturePlot(object.combined, features = top2_features, blend = TRUE)
      ggsave(paste0(outdir,'/Marker/',prefix,'_Cluster_',each,'_co-expression.pdf'), plot1_marker, width = 12, height = 9)
      ggsave(paste0(outdir,'/Marker/',prefix,'_Cluster_',each,'_co-expression.png'), plot1_marker, width = 12, height = 9, type = 'cairo-png')

      plot2_marker <- VlnPlot(object.combined,top_gene_name,ncol =2,pt.size = 0) + xlab("Cluster") + ylab("log(UMI)")
      ggsave(paste0(outdir,'/Marker/',prefix,'_Cluster_',each,'_violin.pdf'), plot2_marker, width = 12, height = 9)
      ggsave(paste0(outdir,'/Marker/',prefix,'_Cluster_',each,'_violin.png'), plot2_marker, width = 12, height = 9, type = 'cairo-png')

      plot3_marker <- FeaturePlot(object.combined,top_gene_name,cols = c("grey", "blue"), reduction = "tsne")
      ggsave(paste0(outdir,'/Marker/',prefix,'_Cluster_',each,'_tsne.pdf'), plot3_marker, width = 12, height = 9)
      ggsave(paste0(outdir,'/Marker/',prefix,'_Cluster_',each,'_tsne.png'), plot3_marker, width = 12, height = 9, type="cairo-png")

      plot4_marker <- FeaturePlot(object.combined,top_gene_name,cols = c("grey", "blue"), reduction = "umap")
      ggsave(paste0(outdir,'/Marker/',prefix,'_Cluster_',each,'_umap.pdf'), plot4_marker,width = 12, height = 9)
      ggsave(paste0(outdir,'/Marker/',prefix,'_Cluster_',each,'_umap.png'), plot4_marker,width = 12, height = 9, type="cairo-png")
    }

    # object.combined %>% saveRDS(paste0(outdir, '/intestine_before_coor_test.rds'))
    # calculate the correlation between the markers
    # AverageExp_top4markers <- AverageExpression(object.combined, features=unique(top_markers$Gene_name),group.by = "seurat_clusters")
    AverageExp_top4markers <- AverageExpression(object.combined, assays = "RNA",features=unique(top_markers$Gene_name), group.by = "seurat_clusters")
    coorda <- corr.test(as.matrix(AverageExp_top4markers$RNA),as.matrix(AverageExp_top4markers$RNA),method="spearman")
    plot_heatmap <- pheatmap(coorda$r)
    ggsave(paste0(outdir,"/Marker/", prefix, '_correlation_top4_markers.pdf'), plot_heatmap, width = 22, height = 13)
    ggsave(paste0(outdir,"/Marker/", prefix, '_correlation_top4_markers.png'), plot_heatmap, width = 22, height = 13, type = 'cairo-png')

    object.combined <- ScaleData(object = object.combined, features = top_gene)
    features_heatmap <- DoHeatmap(object = object.combined, features = top_gene, disp.min = -3, disp.max = 3,label = FALSE)
    ggsave(paste0(outdir, '/Marker/', prefix, '_top4markers_Heatmap.pdf'), features_heatmap, width = 22, height = 13)
    ggsave(paste0(outdir, '/Marker/', prefix, '_top4markers_Heatmap.png'), features_heatmap, width = 22, height = 13, type = 'cairo-png')

    plot_dot <- DotPlot(object.combined, features = top_gene) +
      theme(axis.text.x = element_text(size = 9, angle = 90, vjust = 1, hjust=1),
            plot.margin=unit(rep(1.5,4),'lines'))
    ggsave(paste0(outdir, '/Marker/', prefix, '_top4markers_Dotplot.pdf'), plot_dot, width = 22, height = 13)
    ggsave(paste0(outdir, '/Marker/', prefix, '_top4markers_Dotplot.png'), plot_dot, width = 22, height = 13, type = 'cairo-png')
    } else{
    cat(paste(",cluster,p_val_adj,p_val,avg_log2FC,pct.1,pct.2","\n",sep=""),
        file=paste0(outdir,'/Marker/avg_exp_allmarkers.csv'))
    cat(paste("Number of cells used for clustering,0",
              "\n",sep=""),file=paste0(outdir,'/Marker/avg_exp_allmarkers.csv',sep=""))
    }
  } else{
    cat(paste(",cluster,p_val_adj,p_val,avg_log2FC,pct.1,pct.2","\n",sep=""),
        file=paste0(outdir,'/Marker/avg_exp_allmarkers.csv'))
    cat(paste("Number of cells used for clustering,0",
              "\n",sep=""),file=paste0(outdir,'/Marker/avg_exp_allmarkers.csv',sep=""))
    }

##### cell-type annotation SingleR
# choose to use SingleR or not
if(singler == "no"){
  object.combined %>% saveRDS(paste0(outdir, "/", prefix, "_combined.rds"))
  # save.image(file=paste0(outdir,'/', prefix,'.RData'))
} else{
  dir.create(paste0(outdir,'/singleR'))  # create the singleR directory

  # Perform SingleR for cell types annotation
  library(SingleR)
  library(scater)
  # Seven reference databases (5 for human and 2 for mouse)
  if (species=="GRCh38"){
    #     train=HumanPrimaryCellAtlasData()
    load('/public/work/Pipline/Single_RNA/standard_analysis/Transcriptome/05.SingleR/HUMAN.RData')
    # load('./data/HUMAN.RData')
  } else if(species=="mm10"){
    #     train=MouseRNAseqData()
    load('/public/work/Pipline/Single_RNA/standard_analysis/Transcriptome/05.SingleR/Mouse.RData')
  }

  # Annotation for single cells
  object.combined@meta.data$Barcode <- object.combined@meta.data %>% rownames()
  query  <- as.SingleCellExperiment(object.combined)

  common <- intersect(rownames(query), rownames(train)) # find common genes
  query <-query[common,]
  train <-train[common,]

  # Annotation for single cells
  # if(!is.null(gene)){
  #   pred <- SingleR(query, train, labels = train$label)
  # }else{
  pred <- SingleR(query, train, labels = train$label.main, de.method="wilcox")
  # }

  final <- table(Barcode = colnames(object.combined), Celltype = pred$pruned.labels)
  write.table(final, file = paste0(outdir ,"/singleR/single_cell_final_celltype.txt"), quote = FALSE, sep = "\t")

  pred$final_celltype <- str_replace(str_replace(str_replace(pred$pruned.labels,"\\d$",""),"\\.$",""),"\\.$","")
  celltype_df <- data.frame(Barcode = rownames(pred) ,tuning.scores = pred$tuning.scores,
                            labels = pred$labels, pruned.labels=pred$pruned.labels,
                            first.labels = pred$first.labels, final_celltype = pred$final_celltype)

  celltype_df %>% write.table(file = paste0(outdir ,"/singleR/single_cell_singler_prediction_statistics.txt"), quote = FALSE, sep = "\t",row.names = F)
  Barcode <- as.vector(celltype_df$Barcode)
  final_celltype <- (celltype_df$final_celltype)
  names(final_celltype) <- Barcode
  object.combined@active.ident <- factor(final_celltype)
  object.combined@meta.data$sc_pred_celltype <- final_celltype # Add single-cell cell_type predictions into meta.data

  plot_singleR <- object.combined %>% DimPlot(label= T, reduction = dimension_method, label.size = 2, repel = T)
  ggsave(paste0(outdir, '/singleR/Single_cell_DimPlot.pdf'), plot_singleR, width = 12, height = 9)
  ggsave(paste0(outdir, '/singleR/Single_cell_DimPlot.png'), plot_singleR, width = 12, height = 9, type = 'cairo-png')

  # Annotation for clusters
  # if(!is.null(ref_gene)){
  #   pred2 <- SingleR(query, train, labels = train$label, method = "cluster", clusters = colData(query)$seurat_clusters)
  # }else{
  pred2 <- SingleR(query, train, labels = train$label.main, method = "cluster", clusters = colData(query)$seurat_clusters, de.method="wilcox")
  # }
  #pred2 <- SingleR(query, train, labels = train$label, method = "cluster", clusters = colData(query)$seurat_clusters)

  pred2$final_celltype <-str_replace(str_replace(str_replace(pred2$pruned.labels,"\\d+$",""),"\\.$",""),"\\.$","")
  celltype_df2 <- data.frame(Clusters = rownames(pred2) ,tuning.scores = pred2$tuning.scores,
                             labels = pred2$labels, pruned.labels=pred2$pruned.labels,
                             first.labels = pred2$first.labels, final_celltype = pred2$final_celltype)
  celltype_df2 %>% write.table(file =  paste0(outdir,"/singleR/cluster_singler_prediction_statistics.txt"), quote = FALSE, sep = "\t",row.names = F)

  # change column names
  celltype_df2 <- celltype_df2 %>%
    dplyr::rename(seurat_clusters = Clusters, cluster_pred_celltype = final_celltype) %>%
    select(seurat_clusters, cluster_pred_celltype)
  # object.combined@meta.data <- (merge(object.combined@meta.data, celltype_df2, by = "seurat_clusters"))
  object.combined@meta.data <- object.combined@meta.data %>% left_join(celltype_df2, by = "seurat_clusters")

  rownames(object.combined@meta.data) <- object.combined@meta.data$Barcode
  BB <- as.vector(object.combined@meta.data$Barcode)
  TT <- (object.combined@meta.data$cluster_pred_celltype)
  names(TT) <-BB
  object.combined@active.ident  <- factor(TT)

  plot_singleR_cluster <- object.combined %>% DimPlot(label= T, reduction = dimension_method, label.size = 4, repel = T)
  ggsave(paste0(outdir ,'/singleR/seurat_clustersDimPlot.pdf'), plot_singleR_cluster, width = 12, height = 9)
  ggsave(paste0(outdir ,'/singleR/seurat_clustersDimPlot.png'), plot_singleR_cluster, width = 12, height = 9, type="cairo-png")

  object.combined %>% saveRDS(file=paste0(outdir,'/singleR/', prefix,'_combined.rds'))
}



