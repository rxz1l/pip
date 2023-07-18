library(tidyverse)
library(ggalluvial)
library(PlantPhoneDB)
library(circlize)
library(igraph)
library(ggplot2)
library(pheatmap)
library(tidyverse)
library(ggplotify)
#######两列########################
plotdata = data_tissue %>% data.frame() %>% dplyr::select(c("interaction_pairs","tissue"))

plotdata2 <- to_lodes_form(plotdata[,1:ncol(plotdata)],
                           axes = 1:ncol(plotdata),
                           id = "value")
col<- rep(c('#c3e7d1', '#e58484', '#bf033b', '#f00a36', '#ed3b21', '#ff6908', '#ffc719',
        '#598c14', '#335238', '#4a8594', '#085ef3', '#dbe0e3'), 50)
# pdf("3.pdf")
ggplot(plotdata2, aes(x = x, fill=stratum, label=stratum,
               stratum = stratum, alluvium  = value))+#数据
  geom_flow(width = 0.3,#连线宽度
            curve_type = "sine",#曲线形状，有linear、cubic、quintic、sine、arctangent、sigmoid几种类型可供调整
            alpha = 0.5,#透明度
            color = 'white',#间隔颜色
            size = 0.1)+#间隔宽度
  geom_stratum(width = 0.28)+#图中方块的宽度    
  geom_text(stat = 'stratum', size = 2, color = 'black')+
  scale_fill_manual(values = col)+#自定义颜色
  theme_void()+#主题（无轴及网格线）
  theme(legend.position = 'none')
  geom_alluvium()
# dev.off()


#######将细胞相互作用分为两列放置到两边################################
##将细胞对平分，将一半的细胞对放到第三列
talk1 = plotdata$interaction_pairs %>% unique() %>% .[1:33]
talk2 = plotdata$interaction_pairs %>% unique() %>% .[34:67]
plotdata3 = plotdata %>% filter(interaction_pairs %in% talk1) %>% mutate(interaction_pairs1 = NA)
plotdata4 = plotdata %>% filter(interaction_pairs %in% talk2) %>% 
  rename(interaction_pairs1 = interaction_pairs) %>% 
  mutate(interaction_pairs = NA) %>% 
  dplyr::select(interaction_pairs, tissue, interaction_pairs1)
data = rbind(plotdata3, plotdata4)
##ggalluvial绘制流程
plotdata5 <- to_lodes_form(data[,1:ncol(data)],
                           axes = 1:ncol(data),
                           id = "value")
plotdata6 = na.omit(plotdata5)
col<- rep(c('#c3e7d1', '#e58484', '#bf033b', '#f00a36', '#ed3b21', '#ff6908', '#ffc719',
        '#598c14', '#335238', '#4a8594', '#085ef3', '#dbe0e3'), 50)
# pdf("3.pdf")
ggplot(plotdata6, aes(x = x, fill=stratum, label=stratum,
               stratum = stratum, alluvium  = value))+#数据
  geom_flow(width = 0.3,#连线宽度
            curve_type = "sine",#曲线形状，有linear、cubic、quintic、sine、arctangent、sigmoid几种类型可供调整
            alpha = 0.5,#透明度
            color = 'white',#间隔颜色
            size = 0.1)+#间隔宽度
  geom_stratum(width = 0.28)+#图中方块的宽度    
  geom_text(stat = 'stratum', size = 2, color = 'black')+
  scale_fill_manual(values = col)+#自定义颜色
  theme_void()+#主题（无轴及网格线）
  theme(legend.position = 'none')
  geom_alluvium()
# dev.off()

#######plantphoneDB绘图调整################
####w网络图
CCI_network <- function(interaction_count,mycolor, vertex.label.cex=1.5, edge.label.cex=1,edge.width = 50, title="",edgeLabel=TRUE){
    color <- data.frame(Ligands_cell=unique(interaction_count$Ligands_cell),color=mycolor)
    interaction_count <- interaction_count %>%
        inner_join(color)
    net <- graph_from_data_frame(interaction_count)
    karate_groups <- cluster_optimal(net)
    coords <- layout_in_circle(net, order= order(membership(karate_groups))) 
    E(net)$width  <- E(net)$Number/edge.width
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
###弦图
circos_plot = function(interaction_count, mycolor, cex = 1){
    circos.clear()
    circos.par(start.degree = 90, gap.degree = 4, track.margin = c(-0.1, 
        0.1), points.overflow.warning = FALSE)
    par(mar = rep(0, 4))
    arr.col <- interaction_count %>% select(Ligands_cell, Receptors_cell) %>% 
        unique() %>% mutate(color = "black") %>% as.data.frame()
    chordDiagram(x = interaction_count, grid.col = mycolor, transparency = 0.25, 
        directional = 1, direction.type = c("arrows", "diffHeight"), 
        diffHeight = -0.04, annotationTrack = "grid", annotationTrackHeight = c(0.05, 
            0.1), symmetric = TRUE, link.sort = TRUE, link.arr.col = arr.col, 
        link.arr.length = 0.3, link.largest.ontop = TRUE)
    circos.trackPlotRegion(track.index = 1, bg.border = NA, panel.fun = function(x, 
        y) {
        xlim = get.cell.meta.data("xlim")
        sector.index = get.cell.meta.data("sector.index")
        circos.text(x = mean(xlim), y = 3.2, labels = sector.index, 
            facing = "bending", cex = cex)
    })
}