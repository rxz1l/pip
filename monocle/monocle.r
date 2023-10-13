# load packages ------
library(monocle)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(optparse)
library(purrr)

# parameters ------
option_list = list(
    make_option(c("-m", "--monocle_rds"), type="character", default = NULL, action = "store", help="monocle rds file name"),
    make_option(c("-s", "--sample_rds"), type="character", default = NULL, action = "store", help="sample rds file name"),
    make_option(c("-c", "--cluster"), type="character", default = NULL, action = "store", help="cluster file name"),
    make_option(c("-o", "--output"), type="character", default = NULL, action = "store", help="output dir name"),
    make_option(c("-t", "--column"), type="character", default = NULL, action = "store", help="the column to trajectory plot")
)
opt = parse_args(OptionParser(option_list=option_list))
print(opt)
monocle_rds = opt$monocle_rds
sample_rds = opt$sample_rds
cluster = opt$clsuster
output = opt$output
column = opt$column

# load data ------
data = readRDS(monocle_rds)
sc_data = readRDS(sample_rds)
clus<-read.csv(cluster,header=TRUE)
clus$Cluster <- as.factor(clus$Cluster)
m<-intersect(clus$Barcode,colnames(sc_data))
seurat_match<- subset(x = sc_data ,cells= m)

# tarjectory plot 
columns = str_split(column, ",")
walk(columns, \(x) {
    p_cluster = plot_cell_trajectory(data, color_by = x, cell_size = 0.1,show_backbone =TRUE) + facet_wrap(x)
    p_cluster_all = plot_cell_trajectory(data, color_by = x, cell_size = 0.1,show_backbone =TRUE)
    ggsave(p_cluster, file = str_glue("{output}/cell_trajectory_{x}.pdf"), width = 10, height = 7)
    ggsave(p_cluster, file = str_glue("{output}/cell_trajectory_{x}.png"), width = 10, height = 7)
    ggsave(p_cluster_all, file = str_glue("{output}/cell_trajectory_{x}_all.pdf"), width = 10, height = 7)
    ggsave(p_cluster_all, file = str_glue("{output}/cell_trajectory_{x}_all.png"), width = 10, height = 7)
})
# p_cluster = plot_cell_trajectory(data, color_by = "Celltype", cell_size = 0.1,show_backbone =TRUE) + facet_wrap(~Celltype) 
# p_cluster_all = plot_cell_trajectory(data, color_by = "Celltype", cell_size = 0.1,show_backbone =TRUE)
# ggsave(p_cluster, file = "/public/work/Project/Single_cell/atri_xbc_0921_rxz/04.monocle/result/cell_trajectory_cluster.pdf")
# ggsave(p_cluster, file = "/public/work/Project/Single_cell/atri_xbc_0921_rxz/04.monocle/result/cell_trajectory_cluster.png")
# ggsave(p_cluster_all, file = "/public/work/Project/Single_cell/atri_xbc_0921_rxz/04.monocle/result/cell_trajectory_cluster_all.pdf", width = 10, height = 7)
# ggsave(p_cluster_all, file = "/public/work/Project/Single_cell/atri_xbc_0921_rxz/04.monocle/result/cell_trajectory_cluster_all.png", width = 10, height = 7)

# p_batch = plot_cell_trajectory(data, color_by = "batch", cell_size = 0.1, show_backbone = TRUE) + facet_wrap("batch") 
    
# p_batch_all = plot_cell_trajectory(data, color_by = "batch", cell_size = 0.1, show_backbone = TRUE) 
# ggsave(p_batch, file = "/public/work/Project/Single_cell/atri_xbc_0921_rxz/04.monocle/result/cell_trajectory_batch.pdf")
# ggsave(p_batch, file = "/public/work/Project/Single_cell/atri_xbc_0921_rxz/04.monocle/result/cell_trajectory_batch.png")
# ggsave(p_batch_all, file = "/public/work/Project/Single_cell/atri_xbc_0921_rxz/04.monocle/result/cell_trajectory_batch_all.pdf")
# ggsave(p_batch_all, file = "/public/work/Project/Single_cell/atri_xbc_0921_rxz/04.monocle/result/cell_trajectory_batch_all.png")

# p_species = plot_cell_trajectory(data, color_by = "species", cell_size = 0.1, show_backbone = TRUE) + facet_wrap(~species) 
    
# p_species_all = plot_cell_trajectory(data, color_by = "species", cell_size = 0.1, show_backbone = TRUE) 
# ggsave(p_species, file = "/public/work/Project/Single_cell/atri_xbc_0921_rxz/04.monocle/result/cell_trajectory_species.pdf")
# ggsave(p_species, file = "/public/work/Project/Single_cell/atri_xbc_0921_rxz/04.monocle/result/cell_trajectory_species.png")
# ggsave(p_species_all, file = "/public/work/Project/Single_cell/atri_xbc_0921_rxz/04.monocle/result/cell_trajectory_species_all.pdf")
# ggsave(p_species_all, file = "/public/work/Project/Single_cell/atri_xbc_0921_rxz/04.monocle/result/cell_trajectory_species_all.png")

# pData(data)  = pData(data) %>% mutate(group = case_when(str_detect(batch, "2H") ~ "2H", str_detect(batch, "4H") ~ "4H", 
#                                                         str_detect(batch, "6H") ~ "6H", str_detect(batch, "9H") ~ "9H", 
#                                                         str_detect(batch, "DZ") ~ "DZ"))

# p_gruop = plot_cell_trajectory(data, color_by = "group", cell_size = 0.1, show_backbone = TRUE) + facet_wrap(~group) 
    
# p_group_all = plot_cell_trajectory(data, color_by = "group", cell_size = 0.1, show_backbone = TRUE) 
# ggsave(p_gruop, file = "06.monocle_no_soupx/cell_trajectory_group.pdf")
# ggsave(p_gruop, file = "06.monocle_no_soupx/cell_trajectory_group.png")
# ggsave(p_group_all, file = "06.monocle_no_soupx/cell_trajectory_group_all.pdf")
# ggsave(p_group_all, file = "06.monocle_no_soupx/cell_trajectory_group_all.png")



p_Pseudotime = plot_cell_trajectory(data, color_by = "Pseudotime", cell_size = 0.1, show_backbone = TRUE) 
ggsave(p_Pseudotime, file = str_glue("{output}/cell_trajectory_Pseudotime.pdf"))
ggsave(p_Pseudotime, file = str_glue("{output}/cell_trajectory_Pseudotime.png"))

p_State = plot_cell_trajectory(data, color_by="State", cell_size = 0.1) + facet_wrap(~State)
p_State_all = plot_cell_trajectory(data, color_by="State", cell_size = 0.1)
ggsave(p_State, file = str_glue("{output}/cell_trajectory_state.pdf"))
ggsave(p_State, file = str_glue("{output}/cell_trajectory_state.png"))
ggsave(p_State_all, file = str_glue("{output}/cell_trajectory_state_all.pdf"))
ggsave(p_State_all, file = str_glue("{output}/cell_trajectory_state_all.png"))


###  choose genes that define a cell's progress
disp_table <- dispersionTable(data)


seurat_match <- FindVariableFeatures(seurat_match)
data <- setOrderingFilter(data,seurat_match@assays$RNA@var.features)

p_order_gene = plot_ordering_genes(data)
ggsave(p_order_gene, file = str_glue("{output}/ordering_genes.pdf"))
ggsave(p_order_gene, file = str_glue("{output}/ordering_genes.png"))

write.table(pData(data), file = str_glue("{output}/cell_Pseudotime.csv", row.names = F, quote = F,sep=','))


tsne_xy <- t(reducedDimS(data))
colnames(tsne_xy) <- c('Component1','Component2')
tsne_xy <- cbind(rownames(tsne_xy),tsne_xy)
colnames(tsne_xy)[1] <- 'Barcode'             
write.table(tsne_xy, file = paste(output,'cell_reduction.csv',sep='/'), row.names = F, quote = F,sep=',') 


### Analyzing Branches in Single-Cell Trajectories
BEAM_res <- BEAM(data, branch_point = 1, cores = 10, progenitor_method = "duplicate")  ## branch_point: alternative
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
write.table(BEAM_res, file = paste(output,'gene_related_to_branch.txt',sep='/'), row.names = F, quote = F)
heatmap_gene<-row.names(BEAM_res)[order(BEAM_res$qval)][1:50]

save(heatmap_gene, file = str_glue("{output}/heatmap_gene.Rdata"))
## num_clusters: Number of clusters for the heatmap of branch genes
load("06.monocle_no_soupx/heatmap_gene.Rdata")

# pdf(file = paste("/public/work/Project/Single_cell/atri_xbc_0921_rxz/04.monocle/result","branch_dependent_gene_heatmap.pdf",sep='/'), width = 15, height = 5)
p_branched_heatmap = plot_genes_branched_heatmap(data[heatmap_gene,],
                                          branch_point = 1,
                                          num_clusters = 6,
                                          cores = 4,
                                          use_gene_short_name = T,
                                          show_rownames = T,
                                          return_heatmap = T)

# dev.off()
ggsave(p_branched_heatmap$ph_res, file = str_glue("{output}/branch_dependent_gene_heatmap.pdf"))
ggsave(p_branched_heatmap$ph_res, file = str_glue("{output}/branch_dependent_gene_heatmap.png"))

branched_genes <- row.names(BEAM_res)[order(BEAM_res$qval)][1:5]
save(branched_genes, file = str_glue("{output}/branched_genes.Rdata"))
# load("06.monocle_no_soupx/branched_genes.Rdata")
# pdf(file = paste(outdir,"genes_branched_pseudotime.pdf",sep='/'), width = 9, height = 5)
p5 = plot_genes_branched_pseudotime(data[branched_genes[1:2],],
                       branch_point = 1,
                       color_by = "Cluster",
                       ncol = 1) +
                       theme(plot.title = element_text(hjust = 0.5),legend.position = "right")
# dev.off()
ggsave(p5, file = str_glue("{output}/genes_branched_pseudotime.png"))
ggsave(p5, file = str_glue("{output}/genes_branched_pseudotime.pdf"))

# pdf(file = paste(outdir,"genes_in_pseudotime.pdf",sep='/'), width = 9, height = 5)
p6 = plot_genes_in_pseudotime(data[branched_genes[1:2],],
                       color_by = "Cluster",
                       ncol = 1) +
                       theme(plot.title = element_text(hjust = 0.5),legend.position = "right")

# dev.off()
ggsave(p6, file = str_glue("{output}/genes_in_pseudotime.pdf"))
ggsave(p6, file = str_glue("{output}/genes_in_pseudotime.png"))
