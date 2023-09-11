###
 # @Author: renxz 409368950@qq.com
 # @Date: 2023-09-08 14:28:11
 # @LastEditors: renxz 409368950@qq.com
 # @LastEditTime: 2023-09-11 11:23:01
 # @FilePath: /pipline/hdWGCNA/hdwgcna.r
 # @Description: single cell wgcna 
 # @
 # @Copyright (c) 2023 by ${git_name_email}, All Rights Reserved. 
###

# system setup ----------------------------------------------
##single-cell analysis package
library(Seurat)
library(harmony)
library(UCell)

##plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)
library(igraph)
library(corrplot)
library(optparse)

##co-expression network analysis packages
library(WGCNA)
library(hdWGCNA)

##using the cowplot theme for ggplot
theme_set(theme_cowplot())

##set random seed for reproducibility
set.seed(12345)

# parameter setting --------------------------------------------------------------------
option_list <- list(
  make_option(c("-a", "--gene_select_approach"), type = "character", default = "fraction",
              action = "store", help = "How to select genes? Select variable, fraction, all, or custom."
  ),
  make_option(c("-b", "--fraction"), type = "integer", default = 0.05,
              action = "store", help = "A numeric that determines the minimum cells that a gene must be expressed in order to be included. "
  ),
  make_option(c("-c", "--input_gene_list"), type = "character", default = NULL,
              action = "store", help = "the path to the input gene list, only used if fraction = custom!"
  ),
  make_option(c("-d", "--input_data"), type = "character", default = NULL,
              action = "store", help = "the path to the input data"
  ), 
  make_option(c("-e", "--group_by"), type = "character", default = NULL,
              action = "store", help = "A character vector of Seurat metadata column names representing groups for which metacells will be computed.
                                        If you specify multiple groups, please use comma (,) to separate them."
  ), 
  make_option(c("-f", "--ident_group"), type = "character", default = NULL,
              action = "store", help = "set the Idents of the metacell seurat object"
  ), 
  make_option(c("-g", "--k"), type = "integer", default = 25,
              action = "store", help = "set the Idents of the metacell seurat object"
  ), 
  make_option(c("-h", "--max_shared"), type = "integer", default = 15,
              action = "store", help = "set the Idents of the metacell seurat object"
  ), 
  make_option(c("-i", "--min_cells"), type = "integer", default = 100,
              action = "store", help = "the minimum number of cells in a particular grouping to construct metacells"
  ), 
  make_option(c("-j", "--outdir"), type = "character", default = NULL,
              action = "store", help = "the path to the output directory"
  ), 
  make_option(c("-k", "--group_by_vars"), type = "character", default = NULL,
              action = "store", help = "groups to harmonize, the sample column of data"
  ),  
)

## 解析参数
opt = parse_args(OptionParser(option_list = option_list))
str(opt)
gene_select_approach = opt$gene_select_approach
fraction = opt$fraction
input_gene_list = opt$input_gene_list
input_data = opt$input_data
group_by = opt$group_by
ident_group = opt$ident_group
k = opt$k
max_shared = opt$max_shared
min_cells = opt$min_cells
outdir = opt$outdir
group_by_vars = opt$group_by_vars



# wgcna ----------------------------------------------------------------------------------
## load data
data = readRDS(input_data)

## Set up Seurat object for WGCNA ------------------------------
### gene_select_function :fraction, custom, variable
if(gene_select_approach == "fraction"){
    ####fraction
    seurat_obj <- SetupForWGCNA(
        data,
        gene_select = "fraction", # the gene selection approach
        fraction = fraction, # fraction of cells that a gene needs to be expressed in order to be included
        wgcna_name = "tutorial" # the name of the hdWGCNA experiment
    )
}else if(gene_select_approach == "custom"){
    ####custom
    seurat_obj <- SetupForWGCNA(
        data,
        gene_select = "custom", # the gene selection approach
        gene_list = input_gene_list, # custom gene_list
        wgcna_name = "tutorial" # the name of the hdWGCNA experiment
    )
}else if(gene_select_approach == "variable"){
    ####variable
    data = FindVariableFeatures(data)
    seurat_obj <- SetupForWGCNA(
        data,
        gene_select = "variable", # the gene selection approach
        wgcna_name = "tutorial" # the name of the hdWGCNA experiment
    )
}

## Construct metacells-----------------------------------------------
#> 使用K-最邻近（KNN）算法识别源自相同生物来源样本的小群相似细胞的聚集体（元细胞（metacell）），计算这些聚集体中的细胞的平均或求和表达。
### construct metacells  in each group
#> 作者指出细胞类型中数目较少时，元细胞聚合方法效果不好,可有选择的排除。
#> MetacellsByGroups中参数min_cells来排除小于指定单元数的组，最好设定为大于K，不然若有小于K的分组就会报错`Error in combn(nrow(cell_sample), 2) : n < m`,若报错, 可选择调整参数K，min_cells
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = group_by, # specify the columns in seurat_obj@meta.data to group by
  k = k, # nearest-neighbors parameter 20~75 Smaller k numbers can be used for small data sets
  max_shared = max_shared, # maximum number of shared cells between two metacells
  ident.group = ident_group, # set the Idents of the metacell seurat object
  min_cells = min_cells # the minimum number of cells in a particular grouping to construct metacells
)

### normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)

### extraction Metacell data from seurat_obj
# metacell_obj <- GetMetacellObject(seurat_obj)

### Process the Metacell Seurat Object 
# seurat_obj <- NormalizeMetacells(seurat_obj)
# seurat_obj <- FindVariableFeatures(seurat_obj)
# seurat_obj <- ScaleMetacells(seurat_obj, features=VariableFeatures(seurat_obj))
# seurat_obj <- RunPCAMetacells(seurat_obj, features=VariableFeatures(seurat_obj))
# seurat_obj <- RunHarmonyMetacells(seurat_obj, group.by.vars='sample')
# seurat_obj <- RunUMAPMetacells(seurat_obj, reduction='harmony', dims=1:15)

# p1 <- DimPlotMetacells(seurat_obj, group.by='manual_pred_celltype') + umap_theme() + ggtitle("manual_pred_celltype")
# p2 <- DimPlotMetacells(seurat_obj, group.by='sample') + umap_theme() + ggtitle("sample")

# p1 | p2
celltype = data@meta.data %>% pull(group_by)
## WGCNA analysis by cell grouping ----
walk(celltype, function(x){
  ## Co-expression network analysis -----------------------------------------
  ### Set up the expression matrix
  #> 使用参数`group.by`指定要分析的细胞分组, 使用参数`group_name`指定要分析的分组中的子集 
  seurat_obj <- SetDatExpr(
    seurat_obj,
    group_name = x, # the name of the group of interest in the group.by column
    group.by = group_by, # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
    assay = 'RNA', # using RNA assay
    slot = 'data' # using normalized data
  )

  ### Select a soft threshold
  #### Test different soft powers:
  seurat_obj <- TestSoftPowers(
    seurat_obj,
    group_name = x, # the name of the group of interest in the group.by column
    group.by = group_by # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  )

  #### plot the results:
  plot_list <- PlotSoftPowers(seurat_obj)

  #### assemble with patchwork
  dir.create(str_glue("{outdir}/{x}/figure"), recursive = T)
  p = wrap_plots(plot_list, ncol=2)
  ggsave(p, filename = str_glue("{outdir}/{x}/figure/SoftPowers.pdf"), width=12, height=6)
  ggsave(p, filename = str_glue("{outdir}/{x}/figure/SoftPowers.png"), width=12, height=6)

  #### table of TestSoftPowers
  power_table <- GetPowerTable(seurat_obj)
  select_soft_power = power_table %>% filter(SFT.R.sq > 0.8) %>% slice_min(Power) %>% pull(Power)

  ### Construct a co-expression network:
  seurat_obj <- ConstructNetwork(
    seurat_obj, soft_power=select_soft_power,
    setDatExpr=FALSE,
    tom_name = x, # name of the topoligical overlap matrix written to disk
    tom_outdir = str_glue("{outdir}/{x}/TOM")
  )

  ### Visualization of co-expression networks
  pdf(str_glue('{outdir}/{x}/figure//Visualization of co-expression networks.pdf'))
  PlotDendrogram(seurat_obj, main= str_glue('{x} hdWGCNA Dendrogram'))
  dev.off()

  ### Modular gene
  dir.create(str_glue("{outdir}/{x}/table"), recursive = T)
  modular_gene = seurat_obj@misc$tutorial$wgcna_modules
  write.csv(modular_gene, str_glue('{outdir}/{x}/table/wgcna_modules.csv'))

  # ### inspect the topoligcal overlap matrix (TOM)
  # TOM <- GetTOM(seurat_obj)
  # write.csv(TOM, str_glue('{outdir}/table/{companion cell}/wgcna_modules.csv'))

  ## Module Eigengenes and Connectivity ----
  ### Compute harmonized module eigengenes
  #> Module Eigengenes (MEs) 是一种常用的指标，用于总结整个共表达模块的基因表达谱。
  #> 模块特征基因是通过对包含每个模块的基因表达矩阵的子集执行主成分分析 (PCA) 来计算的, 这些 PCA 矩阵中的第一个PC 就是 ME。
  #### need to run ScaleData first or else harmony throws an error:
  seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))

  #### compute all MEs in the full single-cell dataset
  seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars = group_by_vars ###	groups to harmonize by
  ) #时间较长

  #### harmonized module eigengenes:
  hMEs <- GetMEs(seurat_obj) #每个细胞对于每个模块的特征值
  write.csv(hMEs, str_glue('{outdir}/{x}/table/hMEs.csv'))

  #### module eigengenes:
  MEs <- GetMEs(seurat_obj, harmonized=FALSE)
  write.csv(MEs, str_glue('{outdir}/{x}/table/MEs.csv'))

  ### Compute module connectivity
  #### compute eigengene-based connectivity (kME):
  seurat_obj <- ModuleConnectivity(
    seurat_obj,
    group.by = group_by, group_name = x
  )

  #### plot genes ranked by kME for each module
  num = seurat_obj@misc$tutorial$wgcna_modules$module %>% table %>% length
  vec = sapply(1:num, function(x) x[num %% x == 0]) %>% unlist() %>% na.omit()
  ncol = vec[((vec %>% length()) %/% 2 + 1)]

  p <- PlotKMEs(seurat_obj, ncol = ncol)
  ggsave(p, filename = str_glue('{outdir}/{x}/figure/kMEs.png'))
  ggsave(p, filename = str_glue('{outdir}/{x}/figure/kMEs.pdf'))

  ### Getting the module assignment table
  #### get the module assignment table:
  modules <- GetModules(seurat_obj)
  write.csv(modules, str_glue('{outdir}/{x}/table/modules.csv'))

  #### get hub genes
  hub_df <- GetHubGenes(seurat_obj, n_hubs = n_hubs)
  write.csv(hub_df, str_glue('{outdir}/{x}/table/hub_genes.csv'))

  ### compute gene scoring for the top 25 hub genes by kME for each module
  #### with Seurat method
  seurat_obj <- ModuleExprScore(
    seurat_obj,
    n_genes = 25,
    method='Seurat'
  )

  #### with UCell method
  seurat_obj <- ModuleExprScore(
    seurat_obj,
    n_genes = 25,
    method='UCell'
  )

  ### saveRDS
  saveRDS(seurat_obj, str_glue("{outdir}/{x}/hdWGCNA_object.rds"))

  ## Basic Visualization ----
  ### Module Feature Plots
  #### make a featureplot of hMEs for each module
  plot_list <- ModuleFeaturePlot(
    seurat_obj,
    features='hMEs', # plot the hMEs
    order=TRUE # order so the points with highest hMEs are on top
  )

  #### stitch together with patchwork
  p = wrap_plots(plot_list, ncol = ncol)
  ggsave(p, filename = str_glue('{outdir}/{x}/figure/featureplot.png'))
  ggsave(p, filename = str_glue('{outdir}/{x}/figure/featureplot.pdf'))

  #### make a featureplot of hub scores for each module
  plot_list <- ModuleFeaturePlot(
    seurat_obj,
    features='scores', # plot the hub gene scores
    order='shuffle', # order so cells are shuffled
    ucell = TRUE # depending on Seurat vs UCell for gene scoring
  )

  #### stitch together with patchwork
  p = wrap_plots(plot_list, ncol = ncol)
  ggsave(p, filename = str_glue('{outdir}/{x}/figure/feature_score_plot.png'))
  ggsave(p, filename = str_glue('{outdir}/{x}/figure/feature_score_plot.pdf'))

  ### Module Correlations
  png(str_glue('{outdir}/{x}/figure/Correlations_plot.png'))
  ModuleCorrelogram(seurat_obj)
  dev.off()

  pdf(str_glue('{outdir}/{x}/figure/Correlations_plot.pdf'))
  ModuleCorrelogram(seurat_obj)
  dev.off()

  ### Seurat plotting functions
  #### get hMEs from seurat object
  MEs <- GetMEs(seurat_obj, harmonized = TRUE)
  mods <- colnames(MEs)
  mods <- mods[mods != 'grey']

  #### add hMEs to Seurat meta-data:
  seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

  #### plot with Seurat's DotPlot function
  p <- DotPlot(seurat_obj, features=mods, group.by = group_by)

  #### flip the x/y axes, rotate the axis labels, and change color scheme:
  p <- p +
    coord_flip() +
    RotatedAxis() +
    scale_color_gradient2(high='red', mid='grey95', low='blue')

  #### plot output
  ggsave(p, filename = str_glue('{outdir}/{x}/figure/seurat_dotplot.png'))
  ggsave(p, filename = str_glue('{outdir}/{x}/figure/seurat_dotplot.pdf'))


  ### Network Visualization 
  ModuleNetworkPlot(seurat_obj, outdir = str_glue("{outdir}/{x}/figure/ModuleNetworks"))

  ### Combined hub gene network plots
  #### hubgene network
  pdf(str_glue("{outdir}/{x}/figure/ModuleNetworks/HubGeneNetworkPlot.pdf"))
  HubGeneNetworkPlot(
    seurat_obj,
    n_hubs = 3, n_other=5,
    edge_prop = 0.75,
    mods = 'all'
  )
  dev.off()

  png(str_glue("{outdir}/{x}/figure/ModuleNetworks/HubGeneNetworkPlot.png"))
  HubGeneNetworkPlot(
    seurat_obj,
    n_hubs = 3, n_other=5,
    edge_prop = 0.75,
    mods = 'all'
  )
  dev.off()

  ### Applying UMAP to co-expression networks
  seurat_obj <- RunModuleUMAP(
    seurat_obj,
    n_hubs = 10, # number of hub genes to include for the UMAP embedding
    n_neighbors=15, # neighbors parameter for UMAP
    min_dist=0.1 # min distance between points in UMAP space
  )

  png(str_glue("{outdir}/{x}/figure/ModuleNetworks/ModuleUMAPPlot.png"))
  ModuleUMAPPlot(
    seurat_obj,
    edge.alpha=0.25,
    sample_edges=TRUE,
    edge_prop=0.1, # proportion of edges to sample (20% here)
    label_hubs=2 ,# how many hub genes to plot per module?
    keep_grey_edges=FALSE
  )
  dev.off()

  pdf(str_glue("{outdir}/{x}/figure/ModuleNetworks/ModuleUMAPPlot.pdf"))
  ModuleUMAPPlot(
    seurat_obj,
    edge.alpha=0.25,
    sample_edges=TRUE,
    edge_prop=0.1, # proportion of edges to sample (20% here)
    label_hubs=2 ,# how many hub genes to plot per module?
    keep_grey_edges=FALSE
  )
  dev.off()
})
