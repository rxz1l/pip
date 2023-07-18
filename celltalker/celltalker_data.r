#加载包
suppressMessages({
library(celltalker)
library(conflicted)
library(Seurat)
library(dplyr)
library(stringr)
library(magrittr)
library(biomaRt)
library(purrr)
library(optparse)
library(randomcoloR)
library(PlantPhoneDB)
library(circlize)
library(igraph)
library(ggplot2)
library(pheatmap)
library(tidyverse)
library(ggplotify)
})
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")

#弦图函数
myCCI_network <- function(interaction_count,mycolor, vertex.label.cex=1.5, edge.label.cex=1, title="",edgeLabel=TRUE){
    color <- data.frame(Ligands_cell=unique(interaction_count$Ligands_cell),color=mycolor)
    interaction_count <- interaction_count %>%
        inner_join(color)
    net <- graph_from_data_frame(interaction_count)
    karate_groups <- cluster_optimal(net)
    coords <- layout_in_circle(net, order= order(membership(karate_groups))) 
    E(net)$width  <- E(net)$Number/50 
    V(net)$color <- mycolor[get.vertex.attribute(net)$name]
    #E(net)$color <- mycolor
    if(edgeLabel){
        E(net)$label <- E(net)$Number
    }
    pic <- plot(net, edge.arrow.size=0.5, 
     edge.curved=0.1,
     vertex.label.color="black",
     layout = coords,
  edge.label.cex= edge.label.cex,
     vertex.label.cex=vertex.label.cex,main=title)
}


#传参
option_list <- list(  #构建参数列表
    make_option(c("-d", "--data"), type = "character", default = NULL, help = "数据路径"),  
    make_option(c("-l", "--ligand_receptor_pairs"), type = "character", default = NULL, help = "配受体矩阵"),
    make_option(c("-s", "--split"), type = "character", default = NULL, help = "数据按照某一指标切分（如tissue）,若不需要切分可不输入该l参数"),
    make_option(c("-g", "--metadata_grouping"), type = "character", default = NULL, help = "Serurat数据中定义细胞类型的列名"),
    make_option(c("-n", "--number_cells_required"), type = "integer", default = 100, help = "进行配体/受体相互作用分析,每个分组所需的细胞数。默认为100"),
    make_option(c("-i", "--min_expression"), type = "integer", default = 1000, help = "在相互作用分析中考虑配体或受体的最小表达计数,默认为1000"),
    make_option(c("-a", "--max_expression"), type = "integer", default = 20000, help = "在相互作用分析中考虑配体或受体的最大表达计数,默认为2000"),
    make_option(c("-t", "--scramble_times"), type = "integer", default = 10, help = "配体/受体相互作用的随机排列次数,默认为10"),
    make_option(c("-P", "--PlantPhoneDB"), type = "logical", default = FALSE, help = "是否使用plantphoneDB进行分析，默认为不使用(TRUE/FALSE)"),
    make_option(c("-m", "--method"), type = "character", default = "Average", help = "plantPhoneDB的用法，默认为Average(LRscore/WeightProduct/Average/Product)"),
    make_option(c("-c", "--min.pct"), type = "integer", default = 0.1, help = "只测试在两个群体中任何一个群体中最小比例的细胞中检测到的基因。默认值是0.1"),
    make_option(c("-e", "--iterations"), type = "integer", default = 100, help = "如果方法设置为Average或Product，则使用排列测试通过随机打乱聚类标签来计算配体-受体相互作用得分。默认值为100"),
    make_option(c("-o", "--output_dir"), type = "character", default = system("pwd",intern = TRUE), help = "输出路径,默认为当前目录")
) 
args = parse_args(OptionParser(option_list = option_list))

#读取数据
data = readRDS(args$data)
Idents(data) = args$metadata_grouping
ligand_receptor_pairs = eval(parse(text = load(args$ligand_receptor_pairs)))

#celltalk
##split
if(!is.null(args$split)){
    dir.create(str_glue("{args$output_dir}/celltalk_output/{args$split}"), recursive = T)
    ser_split = SplitObject(data, split = args$split)
    split = data@meta.data %>% select(args$split) %>% table() %>% names()
    data_split = map_dfr(split,function(tissue){
        hca_bm_interactions<- celltalk(input_object=ser_split[[tissue]],
            metadata_grouping=args$metadata_grouping,
            ligand_receptor_pairs=ligand_receptor_pairs,
            number_cells_required=args$number_cells_required,
            min_expression=args$min_expression,
            max_expression=args$max_expression,
            scramble_times=args$scramble_times)
        ## Identify top statistically significant interactions
        top_stats <- hca_bm_interactions %>%
            mutate(fdr=p.adjust(p_val,method="fdr")) %>%
            filter(fdr<0.05) %>%
            group_by(cell_type1) %>%
            top_n(3,interact_ratio) %>%
            ungroup() 
        ## Generate a circos plot
        colors_use <- randomColor(length(unique(ser_split[[tissue]]@meta.data %>% pull(args$metadata_grouping))))
        pdf(str_glue("{args$output_dir}/celltalk_output/{args$split}/{tissue}.pdf"))
        tryCatch({
            circos_plot(ligand_receptor_frame=top_stats,
                cell_group_colors=colors_use,
                ligand_color="blue",
                receptor_color="red",
                cex_outer=0.5,
                cex_inner=0.4)
            },error = function(e){
                message("细胞类型太多了，无法绘制弦图")
            })
        dev.off()
        result <- hca_bm_interactions %>%
            mutate(fdr=p.adjust(p_val,method="fdr")) %>%
            filter(fdr<0.05)
        result %<>% mutate(!!args$split := tissue)
    })
    save(data_split, file = str_glue("{args$output_dir}/celltalk_output/{args$split}/data_split.Rdata"))
    write.csv(data_split, file = str_glue("{args$output_dir}/celltalk_output/{args$split}/data_split.csv"))
}
##all
dir.create(str_glue("{args$output_dir}/celltalk_output/all"), recursive = T)
hca_bm_interactions<- celltalk(input_object=data,
    metadata_grouping=args$metadata_grouping,
    ligand_receptor_pairs=ligand_receptor_pairs,
    number_cells_required=args$number_cells_required,
    min_expression=args$min_expression,
    max_expression=args$max_expression,
    scramble_times=args$scramble_times)
###Identify top statistically significant interactions
top_stats_all <- hca_bm_interactions %>%
    mutate(fdr=p.adjust(p_val,method="fdr")) %>%
    filter(fdr<0.05) %>%
    group_by(cell_type1) %>%
    top_n(3,interact_ratio) %>%
    ungroup()
###Generate a circos plot
colors_use <- randomColor(length(unique(data@meta.data %>% pull(args$metadata_grouping))))
pdf(str_glue("{args$output_dir}/celltalk_output/all/all.pdf"))
tryCatch({
    circos_plot(ligand_receptor_frame=top_stats_all,
        cell_group_colors=colors_use,
        ligand_color="blue",
        receptor_color="red",
        cex_outer=0.5,
        cex_inner=0.4)
    },error = function(e){
        message("细胞类型太多了，无法绘制弦图")
    }
)
dev.off()
data_all = hca_bm_interactions %>%
    mutate(fdr=p.adjust(p_val,method="fdr")) %>%
    filter(fdr<0.05)
save(data_all, file = str_glue("{args$output_dir}/celltalk_output/all/data_all.Rdata"))
write.csv(data_all, file = str_glue("{args$output_dir}/celltalk_output/all/data_all.csv"))
#plantphoneDB
##配受体数据
if(args$PlantPhoneDB){
    LR_pair = ligand_receptor_pairs %>% select(ligand, receptor) %>% rename(Ligands = ligand, Receptors = receptor)
    ##SCT
    data = SCTransform(data, verbose = FALSE)

    if(!is.null(args$split)){
        dir.create(str_glue("{args$output_dir}/plantphone_output/{args$split}"), recursive = T)
        #拆分数据
        ser_split = SplitObject(data, split = args$split)
        split = data@meta.data %>% select(args$split) %>% table() %>% names()
        
        #在split中循环
        plantphone_data_split = map_dfr(split,function(tissue){
            ## 计算得分LRScore
            LRp = LRscore(ser_split[[tissue]]@assays$SCT@data, LRdb=LR_pair, cluster = Idents(data), min.pct = args$min.pct,iterations=args$iterations, method=args$method)
            
            ## 筛选显著的配售体对
            LRp_sig <- LRp %>%
                filter(Pvalue<0.05)

            ## 统计Autocrine和Paracrine个数
            interaction_count <- LRp_sig%>%
                group_by(Ligands_cell,Receptors_cell) %>%
                summarise(Number=n(),.groups = 'drop')
            Autocrine <- interaction_count[interaction_count$Ligands_cell == interaction_count$Receptors_cell,]
            Paracrine <- interaction_count[interaction_count$Ligands_cell != interaction_count$Receptors_cell,]
            
            #绘制弦图
            colors_use <- randomColor(length(unique(ser_split[[tissue]]@meta.data %>% pull(args$metadata_grouping))))
            pdf(str_glue("{args$output_dir}/plantphone_output/{args$split}/{tissue}_Paracrine.pdf"))
                tryCatch({
                    CCI_circle(Paracrine, colors_use)
                },error = function(e){
                    message("细胞类型太多了，无法绘制弦图")
                })
            dev.off()
            pdf(str_glue("{args$output_dir}/plantphone_output/{args$split}/{tissue}_Autocrine.pdf"))
                tryCatch({
                    CCI_circle(Autocrine, colors_use)
                },error = function(e){
                    message("细胞类型太多了，无法绘制弦图")
                })
            dev.off()
            pdf(str_glue("{args$output_dir}/plantphone_output/{args$split}/{tissue}_all.pdf"))
                tryCatch({
                    myCCI_network(interaction_count, colors_use, vertex.label.cex=4, edge.label.cex=1)
                },error = function(e){
                    message("细胞类型太多了，无法绘制弦图")
                })
            dev.off()

            #绘制热图
            pdf(str_glue("{args$output_dir}/plantphone_output/{args$split}/{tissue}_heatmap.pdf"))
                tryCatch({
                    heatmap_count(interaction_count,text_size=5,number_size=2,decimal=4,title="")
                },error = function(e){
                    message("细胞类型太多了，无法绘制弦图")
                })
            dev.off()
            
            #输出数据
            LRp_sig %<>% mutate(!!args$split := tissue)
        })
        save(plantphone_data_split, file = str_glue("{args$output_dir}/plantphone_output/{args$split}/plantphone_data_split.Rdata"))
        write.csv(plantphone_data_split, file = str_glue("{args$output_dir}/plantphone_output/{args$split}/plantphone_data_split.csv"))
    }

    #all
    dir.create(str_glue("{args$output_dir}/plantphone_output/all"), recursive = T)
    LRP_all = LRscore(data@assays$SCT@data, LRdb=LR_pair, cluster = Idents(data), min.pct = args$min.pct,iterations=args$iterations, method=args$method)

    ## 筛选显著的配售体对
    LRp_all_sig <- LRP_all %>%
                filter(Pvalue<0.05)

    ## 统计Autocrine和Paracrine个数
    interaction_count <- LRp_all_sig %>%
        group_by(Ligands_cell,Receptors_cell) %>%
        summarise(Number=n(),.groups = 'drop')
    Autocrine <- interaction_count[interaction_count$Ligands_cell == interaction_count$Receptors_cell,]
    Paracrine <- interaction_count[interaction_count$Ligands_cell != interaction_count$Receptors_cell,]
            
    #绘制弦图
    colors_use <- randomColor(length(unique(data@meta.data %>% pull(args$metadata_grouping))))
    pdf(str_glue("{args$output_dir}/plantphone_output/all/ALL_Paracrine.pdf"))
    tryCatch({
        CCI_circle(Paracrine, colors_use)
        },error = function(e){
        message("细胞类型太多了，无法绘制弦图")
        }
    )
    dev.off()

    pdf(str_glue("{args$output_dir}/plantphone_output/all/ALL_Autocrine.pdf"))
    tryCatch({
        CCI_circle(Autocrine, colors_use)
        },error = function(e){
        message("细胞类型太多了，无法绘制弦图")
        }
    )
    dev.off()

    pdf(str_glue("{args$output_dir}/plantphone_output/all/ALL.pdf"))
    tryCatch({
        myCCI_network(interaction_count, colors_use, vertex.label.cex=1, edge.label.cex=0.5)
        },error = function(e){
        message("细胞类型太多了，无法绘制弦图")
        }
    )
    dev.off()

    #绘制热图
    pdf(str_glue("{args$output_dir}/plantphone_output/all/all_heatmap.pdf"))
    #tryCatch({
    heatmap_count(interaction_count,text_size=4,number_size=1,decimal=3,title="")
        #},error = function(e){
        #message("细胞类型太多了，无法绘制弦图")
        #}
    #)
    dev.off()      
    #输出数据
    plantphone_data_all = LRp_all_sig
    save(plantphone_data_all, file = str_glue("{args$output_dir}/plantphone_output/all/plantphone_data_all.Rdata"))
    write.csv(plantphone_data_all, file = str_glue("{args$output_dir}/plantphone_output/all/plantphone_data_all.csv"))
}