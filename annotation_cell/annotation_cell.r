###
 # @Author: renxz 409368950@qq.com
 # @Date: 2023-09-11 16:01:58
 # @LastEditors: renxz 409368950@qq.com
 # @LastEditTime: 2023-12-27 17:36:05
 # @FilePath: /pipline/annotation_cell/annotation_cell.r
 # @Description: annotation cell
 # @
 # @Copyright (c) 2023 by ${git_name_email}, All Rights Reserved. 
###

## load packages ----
library(Seurat)
library(AUCell)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(optparse)

## 传参 ----
option_list <- list(  #构建参数列表
    make_option(c("-d", "--data_path"), type = "character", default = NULL, help = "数据路径"),  
    make_option(c("-l", "--markerlist"), type = "character", default = NULL, help = "markerlist路径"),
    make_option(c("-o", "--out_dir"), type = "character", default = NULL, help = "输出路径"))
args = parse_args(OptionParser(option_list = option_list))

data_path = args$data_path
markerlist = args$markerlist
out_dir = args$out_dir

## load data ----
data = readRDS(data_path)
markerlist = read.csv(markerlist)
markerlist = markerlist %>% as.list() %>% sapply(unique)

## AUCell ------
### 提取矩阵
exprMatrix = data@assays$RNA@data

### 表达矩阵转rank矩阵
cells_rankings = AUCell_buildRankings(exprMatrix)

### 给每个基因集分配活性阈值
set.seed(123)
AUC = function(markerlist, out_dir, data, ...){
        dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

        cells_AUC = AUCell_calcAUC(markerlist, cells_rankings, nCores = 10, aucMaxRank=nrow(cells_rankings)*0.05)

        cells_assignment = AUCell_exploreThresholds(cells_AUC, plotHist=F, assign=TRUE)

        aucMat = getAUC(cells_AUC) %>% data.frame() %>% rownames_to_column(var = "celltype") %>% pivot_longer(-celltype, names_to = "samplename", values_to = "AUCCell")

        cluster = data@meta.data %>% dplyr::select(seurat_clusters) %>% rownames_to_column(var = "samplename") %>% mutate(samplename = str_replace_all(samplename, "-", "\\."))

        cluster_auc = inner_join(cluster, aucMat, by = "samplename")

        cluster_celltype = cluster_auc %>% dplyr::group_by(celltype, seurat_clusters) %>% dplyr::summarise(AUCCell_median = mean(AUCCell)) 

        plot_data = cluster_celltype %>% pivot_wider(names_from = seurat_clusters, values_from = AUCCell_median) %>% column_to_rownames(var = "celltype")

        plot_data1 = t(scale(t(plot_data),scale = T,center = T))
        pheatmap(plot_data1, filename = str_glue("{out_dir}/celltype.pdf"))
        
        plot_data2 = scale(plot_data,scale = T,center = T)
        pheatmap(plot_data2, filename = str_glue("{out_dir}/cluster.pdf"))

        return(cells_AUC)
}

## pheatmap 
cells_AUC = AUC(markerlist = markerlist, data = data, out_dir = out_dir)

## umap plot
auc_umapplot = function(markerlist, cells_AUC, data, out_dir, ...) {
        dir.create(str_glue("{out_dir}/umap/AUCell"),recursive = T)
        umapdata = getAUC(cells_AUC) %>% data.frame() %>% t()
        rownames(umapdata) = rownames(umapdata) %>% str_replace_all("\\.", "-")
        data@meta.data = cbind(data@meta.data, umapdata)

        walk(colnames(umapdata), function(x) {
                p = FeaturePlot(data, reduction = "umap", features = x)
                ggsave(p, filename = str_glue("{out_dir}/umap/AUCell/{x}.pdf"))
                ggsave(p, filename = str_glue("{out_dir}/umap/AUCell/{x}.png"))
        })
}

auc_umapplot(markerlist = markerlist, cells_AUC, data, out_dir = out_dir)

# dot plot ------
dotplot = function(markerlist, data, out_dir, ...) {
    out_dir = file.path(out_dir, "dotpot")
    dir.create(out_dir, recursive = T)
    walk(names(markerlist), \(x) {
        tryCatch({
            p = DotPlot(data, features = markerlist[[x]]) + coord_flip()
            ggsave(p, filename = str_glue("{out_dir}/{x}.pdf"))
            ggsave(p, filename = str_glue("{out_dir}/{x}.png"))
        }, error = function(e) {
            print("None of the requested variables were found")
        })
    })
}

dotplot(markerlist, data, out_dir)
