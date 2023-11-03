library(Seurat)
# library(tidyseurat)
library(tidyverse)
library(dplyr)
library(stringr)
library(argparse)
library(reshape2)
library(tidyr)
library(ggplot2)
library(scales)
library(cowplot)
library(reticulate)
library(psych)
library(pheatmap)
library(dittoSeq)

parser = ArgumentParser()
parser$add_argument("--path", help="/path to singleR directory/rds")
# parser$add_argument("--compare", help="samplevssamplevssample,using vs as the split")
parser$add_argument("--outdir", help="outdir of project")
parser$add_argument("--prefix", help="prefix of results")
parser$add_argument("--celltype", help="column of celltype")
parser$add_argument("--grouplist", help="plot of group")



args <- parser$parse_args()
str(args)

path=args$path
# compare=args$compare
outdir=args$outdir
prefix=args$prefix
celltype=args$celltype
grouplist=args$grouplist
grouplist = str_split(grouplist, ",")[[1]]

###--------- Create directories
if (!dir.exists(outdir)){
  dir.create(outdir)
}

if (!dir.exists(paste0(outdir,'/Final_Marker'))){
  dir.create(paste0(outdir,'/Final_Marker'))
}

if (!dir.exists(paste0(outdir,'/CellsRatio'))){
  dir.create(paste0(outdir,'/CellsRatio'))
}

if (!dir.exists(paste0(outdir,'/UMAP'))){
  dir.create(paste0(outdir,'/UMAP'))
}

###---------- set colors
friendly_cols <- dittoSeq::dittoColors()
custom_theme <-
  list(
    scale_fill_manual(values = friendly_cols),
    scale_color_manual(values = friendly_cols),
    theme_bw() +
      theme(
        panel.border = element_blank(),
        axis.line = element_line(),
        panel.grid.major = element_line(size = 0.2),
        panel.grid.minor = element_line(size = 0.1),
        text = element_text(size = 18),
        legend.position = "right",
        strip.background = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10,r = 10,b = 10,l = 10)),
        axis.title.y = element_text(margin = margin(t = 10,r = 10,b = 10,l = 10)),
        axis.text.x = element_text(angle = 30,hjust = 1,vjust = 1,)))

###---------- Start to read files
object.combined <- readRDS(file = path) # read rds

DefaultAssay(object.combined) <- "RNA"
Idents(object.combined) <- celltype

####---------- 1st. 3 UMAP plots. (celltype; colored by group; split by group)

umap_1 <- object.combined %>% DimPlot(reduction = "umap", label = T,label.size = 6, repel = T) +
  custom_theme # celltype UMAP +
# umap_2 <- object.combined %>% DimPlot(reduction = "umap", group.by = "group") +
#   custom_theme # color by groups
# umap_split <- object.combined %>% DimPlot(reduction = "umap", split.by = "group", ncol = 2) +
#   custom_theme # Split by group
# umap_sample <- object.combined %>% DimPlot(reduction = "umap", group.by = "batch", ncol = 2) +
#   custom_theme 
# umap_split_2 <- object.combined %>% DimPlot(reduction = "umap", split.by = "batch", ncol = 2) +
#   custom_theme # Split by sample
# umap_tissue <- object.combined %>% DimPlot(reduction = "umap", group.by = "tissue", ncol = 2) +
#   custom_theme # color by tissue
# umap_tissue_split <- object.combined %>% DimPlot(reduction = "umap", split.by = "tissue", ncol = 2) +
#   custom_theme # color by tissue
walk(grouplist, \(x) {
  umap <- object.combined %>% DimPlot(reduction = "umap", group.by = x) +
    custom_theme # color by groups
  umap_split <- object.combined %>% DimPlot(reduction = "umap", split.by = x, ncol = 2) +
    custom_theme # Split by group
  ggsave(paste0(outdir,'/UMAP/', prefix, "_", x, '_UMAP.pdf'), umap, width = 12, height = 9)
  ggsave(paste0(outdir,'/UMAP/', prefix, "_", x, '_UMAP.png'), umap, type="cairo-png", width = 12, height = 9)
  ggsave(paste0(outdir,'/UMAP/', prefix, "_", x, '_split_UMAP.pdf'), umap_split, width = 22, height = 10)
  ggsave(paste0(outdir,'/UMAP/', prefix, "_", x, '_split_UMAP.png'), umap_split, type="cairo-png", width = 22, height = 10)
})

ggsave(paste0(outdir,'/UMAP/', prefix, '_celltype_UMAP.pdf'), umap_1, width = 12, height = 9)
ggsave(paste0(outdir,'/UMAP/', prefix, '_celltype_UMAP.png'), umap_1, type="cairo-png", width = 12, height = 9)
# ggsave(paste0(outdir,'/UMAP/', prefix, '_group_UMAP.pdf'), umap_2, width = 12, height = 9)
# ggsave(paste0(outdir,'/UMAP/', prefix, '_group_UMAP.png'), umap_2, type="cairo-png", width = 12, height = 9)
# ggsave(paste0(outdir,'/UMAP/', prefix, '_group_split_UMAP.pdf'), umap_split, width = 22, height = 10)
# ggsave(paste0(outdir,'/UMAP/', prefix, '_group_split_UMAP.png'), umap_split, type="cairo-png", width = 22, height = 10)
# ggsave(paste0(outdir,'/UMAP/', prefix, '_sample_UMAP.pdf'), umap_sample, width = 12, height = 9)
# ggsave(paste0(outdir,'/UMAP/', prefix, '_sample_UMAP.png'), umap_sample, type="cairo-png", width = 12, height = 9)
# ggsave(paste0(outdir,'/UMAP/', prefix, '_sample_split_UMAP.pdf'), umap_split_2, width = 22, height = 10)
# ggsave(paste0(outdir,'/UMAP/', prefix, '_sample_split_UMAP.png'), umap_split_2, type="cairo-png", width = 22, height = 10)
# ggsave(paste0(outdir,'/UMAP/', prefix, '_tissue_UMAP.pdf'), umap_tissue, width = 12, height = 9)
# ggsave(paste0(outdir,'/UMAP/', prefix, '_tissue_UMAP.png'), umap_tissue, type="cairo-png", width = 12, height = 9)
# ggsave(paste0(outdir,'/UMAP/', prefix, '_tissue_split_UMAP.pdf'), umap_tissue_split, width = 22, height = 10)
# ggsave(paste0(outdir,'/UMAP/', prefix, '_tissue_split_UMAP.png'), umap_tissue_split, type="cairo-png", width = 22, height = 10)

###----------- 2nd. Calculate related information of cell ratios
### 计算每个sample中有多少个细胞类型（num 和%）在每个celltype中
### samples celltype percent
### group by each (celltype, sample) combination

###group by celltype
walk(grouplist, \(x) {
  data = prop.table(table(object.combined$celltype,object.combined@meta.data[[x]]),1)
  data <- data %>% as.data.frame() %>% dplyr::rename("celltype" = "Var1",  {{x}} := "Var2", "ratio" = "Freq")
  sample_clus_perc %>% write.table(file=paste0(outdir,'/CellsRatio/',prefix,'_each_celltype_', x, '_perc.csv'),sep=',',row.names = F, quote=F)
  
  barplot4 <- data %>%
    ggplot(aes_string(x = "celltype", y = "ratio", fill = x)) +
    geom_bar(position = 'stack', stat = 'identity') +
    labs(x = "Cell Type", y = "Cell Ratios") +
    theme(panel.background=element_rect(fill='transparent', color='black'),
        legend.key=element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(color="black")) +
    scale_y_continuous(expand=c(0.001,0.001)) +
  # scale_fill_manual(values=colour1) +
    guides(fill = guide_legend(keywidth = 1, keyheight = 1,ncol=1,title = x)) +
    custom_theme +
    ggtitle(str_glue("Percentage of cells in each cell type for different {x}"))

  ggsave(filename = paste0(outdir,'/CellsRatio/',prefix,'_each_celltype_', x, '_percent.pdf'), barplot4, dpi = 320,width = 12, height = 9)
  ggsave(filename = paste0(outdir,'/CellsRatio/',prefix,'_each_celltype_', x, '_percent.png'), barplot4, dpi = 320, width = 12, height = 9, type="cairo-png")

})

walk(grouplist, \(x) {
  data = prop.table(table(object.combined@meta.data[[x]], object.combined$celltype),1)
  data <- data %>% as.data.frame() %>% dplyr::rename({{x}} := "Var1", "celltype" = "Var2", "ratio" = "Freq")
  sample_clus_perc %>% write.table(file=paste0(outdir,'/CellsRatio/',prefix, '_each_', x, '_celltype', '_perc.csv'),sep=',',row.names = F, quote=F)
  
  barplot4 <- data %>%
    ggplot(aes_string(x = x, y = "ratio", fill = "celltype")) +
    geom_bar(position = 'stack', stat = 'identity') +
    labs(x = "Cell Type", y = "Cell Ratios") +
    theme(panel.background=element_rect(fill='transparent', color='black'),
        legend.key=element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(color="black")) +
    scale_y_continuous(expand=c(0.001,0.001)) +
  # scale_fill_manual(values=colour1) +
    guides(fill = guide_legend(keywidth = 1, keyheight = 1,ncol=1,title = x)) +
    custom_theme +
    ggtitle(str_glue("Percentage of cells in each cell type for different {x}"))

  ggsave(filename = paste0(outdir,'/CellsRatio/',prefix, '_each_', x, '_celltype', '_percent.pdf'), barplot4, dpi = 320,width = 12, height = 9)
  ggsave(filename = paste0(outdir,'/CellsRatio/',prefix, '_each_', x, '_celltype', '_percent.png'), barplot4, dpi = 320, width = 12, height = 9, type="cairo-png")

})




# sample_clus_count <- table(object.combined$batch,object.combined$celltype) # for each sample, how many cell types
# sample_clus_perc <- prop.table(table(object.combined$batch,object.combined$celltype),1) # for each sample, the cell type percentages
# group_clus_perc <- prop.table(table(object.combined$group, object.combined$celltype),1)
# clus_sample_perc <- prop.table(table(object.combined$celltype, object.combined$batch),1)
# clus_group_relative_perc <- prop.table(table(object.combined$celltype, object.combined$group),1)
# # clus_group_absolute_perc <- prop.table(table(object.combined$celltype, object.combined$group))
# clus_tissue_perc <- prop.table(table(object.combined$celltype, object.combined$tissue),1)
# tissue_clus_perc <- prop.table(table(object.combined$tissue, object.combined$celltype),1)


# # Save to files
# sample_clus_count <- sample_clus_count %>% as.data.frame() %>% dplyr::rename(sample = Var1, celltype = Var2, num = Freq) 
# sample_clus_count %>% write.table(file=paste0(outdir,'/CellsRatio/',prefix,'_each_sample_celltype_abundance.csv'),sep=',',row.names = F, quote=F)
# sample_clus_perc <- sample_clus_perc %>% as.data.frame() %>% dplyr::rename(sample = Var1, celltype = Var2, ratio = Freq) 
# sample_clus_perc %>% write.table(file=paste0(outdir,'/CellsRatio/',prefix,'_each_sample_celltype_perc.csv'),sep=',',row.names = F, quote=F)
# group_clus_perc <- group_clus_perc %>% as.data.frame() %>% dplyr::rename(group = Var1, celltype = Var2, ratio = Freq) 
# group_clus_perc %>% write.table(file=paste0(outdir,'/CellsRatio/',prefix,'_each_group_celltype_perc.csv'),sep=',',row.names = F, quote=F)
# clus_sample_perc <- clus_sample_perc %>% as.data.frame() %>% dplyr::rename(celltype = Var1, sample = Var2, ratio = Freq) 
# clus_sample_perc %>% write.table(file=paste0(outdir,'/CellsRatio/',prefix,'_each_celltype_sample_perc.csv'),sep=',',row.names = F, quote=F)
# clus_group_relative_perc <- clus_group_relative_perc %>% as.data.frame() %>% dplyr::rename(celltype = Var1, group = Var2, ratio = Freq)
# clus_group_relative_perc %>% write.table(file=paste0(outdir,'/CellsRatio/',prefix,'_each_celltype_group_perc.csv'),sep=',',row.names = F, quote=F)
# #clus_group_absolute_perc <- clus_group_absolute_perc %>% as.data.frame() %>% dplyr::rename(celltype = Var1, group = Var2, ratio = Freq) %>% mutate(cell_type = prefix)
# #clus_group_absolute_perc %>% write.table(file=paste0(outdir,'/CellsRatio/',prefix,'_each_celltype_group_absolute_perc.csv'),sep=',',row.names = F, quote=F)
# clus_tissue_perc <- clus_tissue_perc %>% as.data.frame() %>% dplyr::rename(celltype = Var1, tissue = Var2, ratio = Freq)
# clus_tissue_perc %>% write.table(file=paste0(outdir,'/CellsRatio/',prefix,'_each_celltype_tissue_perc.csv'),sep=',',row.names = F, quote=F)
# tissue_clus_perc <- tissue_clus_perc %>% as.data.frame() %>% dplyr::rename(tissue = Var1, celltype = Var2, ratio = Freq)
# tissue_clus_perc %>% write.table(file=paste0(outdir,'/CellsRatio/',prefix,'_each_tissue_celltype_perc.csv'),sep=',',row.names = F, quote=F)


# # 1st: stack position, by celltype
# barplot1 <- clus_group_relative_perc %>%
#   ggplot(aes(x = celltype, y = ratio, fill = group)) +
#   geom_bar(position = 'stack', stat = 'identity') +
#   labs(x = "Cell Type", y = "Cell Ratios") +
#   theme(panel.background=element_rect(fill='transparent', color='black'),
#         legend.key=element_rect(fill='transparent', color='transparent'),
#         axis.text = element_text(color="black")) +
#   scale_y_continuous(expand=c(0.001,0.001)) +
#   # scale_fill_manual(values=colour1) +
#   guides(fill = guide_legend(keywidth = 1, keyheight = 1,ncol=1,title = 'Group')) +
#   custom_theme +
#   ggtitle("Percentage of cells in each group for different cell types")

# ggsave(filename = paste0(outdir,'/CellsRatio/',prefix,'_each_celltype_group_percent_stack.pdf'), barplot1, dpi = 320,width = 12, height = 9)
# ggsave(filename = paste0(outdir,'/CellsRatio/',prefix,'_each_celltype_group_percent_stack.png'), barplot1, dpi = 320, width = 12, height = 9, type="cairo-png")

# # 2nd: dodge position, by celltype
# barplot2 <- clus_group_relative_perc %>%
#   ggplot(aes(x = celltype, y = ratio, fill = group)) +
#   geom_bar(position = 'dodge', stat = 'identity') +
#   labs(x = "Cell Type", y = "Cell Ratios") +
#   theme(panel.background=element_rect(fill='transparent', color='black'),
#         legend.key=element_rect(fill='transparent', color='transparent'),
#         axis.text = element_text(color="black")) +
#   scale_y_continuous(expand=c(0.001,0.001)) +
#   # scale_fill_manual(values=colour1) +
#   guides(fill = guide_legend(keywidth = 1, keyheight = 1,ncol=1,title = 'Group')) +
#   custom_theme +
#   ggtitle("Percentage of cells in each group for different cell types")

# ggsave(filename = paste0(outdir,'/CellsRatio/',prefix,'_each_celltype_group_percent_dodge.pdf'), barplot2, dpi = 320,width = 12, height = 9)
# ggsave(filename = paste0(outdir,'/CellsRatio/',prefix,'_each_celltype_group_percent_dodge.png'), barplot2, dpi = 320, width = 12, height = 9, type="cairo-png")

# # 3rd: stack position, by sample
# barplot3 <- sample_clus_perc %>%
#   ggplot(aes(x = sample, y = ratio, fill = celltype)) +
#   geom_bar(position = 'stack', stat = 'identity') +
#   labs(x = "Cell Type", y = "Cell Ratios") +
#   theme(panel.background=element_rect(fill='transparent', color='black'),
#         legend.key=element_rect(fill='transparent', color='transparent'),
#         axis.text = element_text(color="black")) +
#   scale_y_continuous(expand=c(0.001,0.001)) +
#   # scale_fill_manual(values=colour1) +
#   guides(fill = guide_legend(keywidth = 1, keyheight = 1,ncol=1,title = 'Cell Type')) +
#   custom_theme +
#   ggtitle("Percentage of cells in each cell type for different samples")

# ggsave(filename = paste0(outdir,'/CellsRatio/',prefix,'_each_sample_celltype_percent.pdf'), barplot3, dpi = 320,width = 12, height = 9)
# ggsave(filename = paste0(outdir,'/CellsRatio/',prefix,'_each_sample_celltype_percent.png'), barplot3, dpi = 320, width = 12, height = 9, type="cairo-png")

# barplot4 <- clus_tissue_perc %>%
#   ggplot(aes(x = celltype, y = ratio, fill = tissue)) +
#   geom_bar(position = 'stack', stat = 'identity') +
#   labs(x = "Cell Type", y = "Cell Ratios") +
#   theme(panel.background=element_rect(fill='transparent', color='black'),
#         legend.key=element_rect(fill='transparent', color='transparent'),
#         axis.text = element_text(color="black")) +
#   scale_y_continuous(expand=c(0.001,0.001)) +
#   # scale_fill_manual(values=colour1) +
#   guides(fill = guide_legend(keywidth = 1, keyheight = 1,ncol=1,title = 'Tissue')) +
#   custom_theme +
#   ggtitle("Percentage of cells in each cell type for different tissues")

# ggsave(filename = paste0(outdir,'/CellsRatio/',prefix,'_each_celltype_tissue_percent.pdf'), barplot4, dpi = 320,width = 12, height = 9)
# ggsave(filename = paste0(outdir,'/CellsRatio/',prefix,'_each_celltype_tissue_percent.png'), barplot4, dpi = 320, width = 12, height = 9, type="cairo-png")

# barplot5 <- tissue_clus_perc %>%
#   ggplot(aes(x = tissue, y = ratio, fill = celltype)) +
#   geom_bar(position = 'stack', stat = 'identity') +
#   labs(x = "Tissue", y = "Cell Type") +
#   theme(panel.background=element_rect(fill='transparent', color='black'),
#         legend.key=element_rect(fill='transparent', color='transparent'),
#         axis.text = element_text(color="black")) +
#   scale_y_continuous(expand=c(0.001,0.001)) +
#   # scale_fill_manual(values=colour1) +
#   guides(fill = guide_legend(keywidth = 1, keyheight = 1,ncol=1,title = 'Cell Type')) +
#   custom_theme +
#   ggtitle("Percentage of cells in each tissue for different cell types")

###--------------------------------------------------------------------------------------------
###----------- 3rd. Find markers in each cell type
object.combined@misc$final_markers <- FindAllMarkers(object = object.combined, only.pos = F, min.pct = 0.25)

combined.markers <- object.combined@misc$final_markers
colnames(combined.markers)[7] <- 'Gene_name'
combined.markers <- combined.markers[,c(7,6,2,1,5,3,4)]

# calculate the average expression of all the DEGs in each celltype and save the results into txt/csv.
all_markers <- combined.markers$Gene_name
gene_cell_exp <- object.combined %>%
  Seurat::AverageExpression(assays = "RNA", features = all_markers, group.by = celltype)
gene_cell_exp <- gene_cell_exp$RNA %>% as.data.frame() %>% rownames_to_column(var = "Gene_name")
combined.markers.avg.df <- combined.markers %>% left_join(gene_cell_exp, by = "Gene_name")

combined.markers.avg.df %>% write.table(paste0(outdir, '/Final_Marker/avg_exp_allmarkers.csv'), quote = F, row.names = F, col.names = T, sep = ",")
combined.markers.avg.df %>% write.table(paste0(outdir, '/Final_Marker/avg_exp_allmarkers.txt'), quote = F, row.names = F, col.names = T, sep = "\t")

# select top 10 markers based on the avg_log2FC values
# save all the related statistics of the top 10 markers for each celltype into txt/csv
combined.markers.avg.df %>% group_by(cluster) %>% top_n(10, avg_log2FC) %>%
  write.table(file=paste0(outdir,'/Final_Marker/','top10_marker_each_celltype.txt'),quote=F,row.names=F,col.names=T,sep='\t')
combined.markers.avg.df %>% group_by(cluster) %>% top_n(10, avg_log2FC) %>%
  write.table(file=paste0(outdir,'/Final_Marker/','top10_marker_each_celltype.csv'),quote=F,row.names=F,col.names=T,sep=',')

# select top 4  markers
top_markers <- combined.markers %>% group_by(cluster) %>% top_n(4, avg_log2FC)
top_gene <- top_markers$Gene_name %>% unique()
clus<-unique(top_markers$cluster)

# violin and tsne graph
for (each in clus){
  top_clus <- top_markers %>% filter(cluster == each)
  top_gene_name <- top_clus$Gene_name
  top2_features<-(top_clus$Gene_name)[1:2]

  name <- gsub(" ","_",each)

  plot1_marker <- FeaturePlot(object.combined, features = top2_features, blend = TRUE)
  ggsave(paste0(outdir,'/Final_Marker/',prefix,'_',name,'_co-expression.pdf'), plot1_marker, width = 12, height = 9)
  ggsave(paste0(outdir,'/Final_Marker/',prefix,'_',name,'_co-expression.png'), plot1_marker, width = 12, height = 9, type = 'cairo-png')

  plot2_marker <- VlnPlot(object.combined,top_gene_name,ncol =2,pt.size = 0) + xlab("Cluster") + ylab("log(UMI)")
  ggsave(paste0(outdir,'/Final_Marker/',prefix,'_',name,'_violin.pdf'), plot2_marker, width = 12, height = 9)
  ggsave(paste0(outdir,'/Final_Marker/',prefix,'_',name,'_violin.png'), plot2_marker, width = 12, height = 9, type = 'cairo-png')

  plot3_marker <- FeaturePlot(object.combined,top_gene_name,cols = c("grey", "blue"), reduction = "tsne")
  ggsave(paste0(outdir,'/Final_Marker/',prefix,'_',name,'_tsne.pdf'), plot3_marker, width = 12, height = 9)
  ggsave(paste0(outdir,'/Final_Marker/',prefix,'_',name,'_tsne.png'), plot3_marker, width = 12, height = 9, type="cairo-png")

  plot4_marker <- FeaturePlot(object.combined,top_gene_name,cols = c("grey", "blue"), reduction = "umap")
  ggsave(paste0(outdir,'/Final_Marker/',prefix,'_',name,'_umap.pdf'), plot4_marker,width = 12, height = 9)
  ggsave(paste0(outdir,'/Final_Marker/',prefix,'_',name,'_umap.png'), plot4_marker,width = 12, height = 9, type="cairo-png")
}

# calculate the correlation between the markers
AverageExp_top4markers <- AverageExpression(object.combined, assays = "RNA",features=unique(top_markers$Gene_name),group.by = celltype)
coorda <- corr.test(as.matrix(AverageExp_top4markers$RNA),as.matrix(AverageExp_top4markers$RNA),method="spearman")
plot_heatmap <- pheatmap(coorda$r)
ggsave(paste0(outdir,"/Final_Marker/", prefix, '_correlation_top4_markers.pdf'), plot_heatmap, width = 12, height = 9)
ggsave(paste0(outdir,"/Final_Marker/", prefix, '_correlation_top4_markers.png'), plot_heatmap, width = 12, height = 9, type = 'cairo-png')

object.combined <- ScaleData(object = object.combined, features = top_gene)
features_heatmap <- DoHeatmap(object = object.combined, features = top_gene, disp.min = -3, disp.max = 3,label = FALSE)
ggsave(paste0(outdir, '/Final_Marker/', prefix, '_top4markers_Heatmap.pdf'), features_heatmap, width = 12, height = 9)
ggsave(paste0(outdir, '/Final_Marker/', prefix, '_top4markers_Heatmap.png'), features_heatmap, width = 12, height = 9, type = 'cairo-png')

plot_dot <- DotPlot(object.combined, features = top_gene, assay = "RNA") +
  theme(axis.text.x = element_text(size = 9, angle = 90, vjust = 1, hjust=1),
        plot.margin=unit(rep(1.5,4),'lines'))
ggsave(paste0(outdir, '/Final_Marker/', prefix, '_top4markers_Dotplot.pdf'), plot_dot, width = 12, height = 9)
ggsave(paste0(outdir, '/Final_Marker/', prefix, '_top4markers_Dotplot.png'), plot_dot, width = 12, height = 9, type = 'cairo-png')

