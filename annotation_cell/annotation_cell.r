###
 # @Author: renxz 409368950@qq.com
 # @Date: 2023-09-11 16:01:58
 # @LastEditors: renxz 409368950@qq.com
 # @LastEditTime: 2023-09-11 17:50:24
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

## load data ----
data = readRDS()
markerlist = 

## umap polt ----
### addmodelscore ----
data_AddModuleScore = AddModuleScore(data, features = markerlist, name = names(markerlist))
celltype = data_AddModuleScore %>% head %>% colnames %>% .[8:length(.)]

walk(celltype, function(x){
        p = FeaturePlot(data_AddModuleScore, reduction = "umap", features = x)
        ggsave(p, filename = str_glue("/public/work/Project/Single_cell/qiaomai_renxiaozhen/11.CellAnnotation/umap/AddModuleScore/{x}.pdf"))
        ggsave(p, filename = str_glue("/public/work/Project/Single_cell/qiaomai_renxiaozhen/11.CellAnnotation/umap/AddModuleScore/{x}.png"))
})

### AUC plot ----
#### 提取矩阵
exprMatrix = data@assays$RNA@data

#### 表达矩阵转rank矩阵
cells_rankings = AUCell_buildRankings(exprMatrix)

#### 给每个基因集分配活性阈值
set.seed(123)
cells_AUC = AUCell_calcAUC(markerlist, cells_rankings, nCores = 10, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment = AUCell_exploreThresholds(cells_AUC, plotHist=F, assign=TRUE)

umapdata = getAUC(cells_AUC) %>% data.frame() %>% t()
rownames(umapdata) = rownames(umapdata) %>% str_replace("\\.", "-")
data@meta.data = cbind(data@meta.data, umapdata)

walk(names(markerlist), function(x) {
        p = FeaturePlot(data, reduction = "umap", features = x)
        ggsave(p, filename = str_glue("/public/work/Project/Single_cell/qiaomai_renxiaozhen/11.CellAnnotation/umap/AUCell/{x}.pdf"))
        ggsave(p, filename = str_glue("/public/work/Project/Single_cell/qiaomai_renxiaozhen/11.CellAnnotation/umap/AUCell/{x}.png"))
})

