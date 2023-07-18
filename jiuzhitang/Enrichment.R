#/public/work/Project/Single_cell/SA2022040201_junke2brain_dsh/02.seurat/Brain/Sample_diff/default/bin/run.Astrocytes.r

#/public/work/Personal/chenzhiqiang/pipline/pip_public/pipeline/Annotation/scripts/rpt_mask_fasta.pl
#repeat改小写的

#library(patchwork)
#library(dplyr)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(ggplot2)
library(clusterProfiler)
#library(DO.db)
library(argparse)
library(tidyverse)
library(pathview)
library(dplyr)
library(ggvenn)


####function
enrich_function = function(genelist, species, cell_type, ud, outdir){
  if (species=='GRCh38'){
  library(org.Hs.eg.db)
  gs = bitr(genelist$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  ego.all = enrichGO(gene=gs$ENTREZID, OrgDb = org.Hs.eg.db,ont= "ALL",pAdjustMethod = "BH",pvalueCutoff= 1,qvalueCutoff= 1,readable= TRUE,minGSSize=1)
  kk <- enrichKEGG(gene= gs$ENTREZID, organism = 'hsa',pvalueCutoff = 1,qvalueCutoff=1,minGSSize=1)
  kkk <- kk
}else if (species=='mm10'){
  library(org.Mm.eg.db)
  gs = bitr(genelist$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
  ego.all = enrichGO(gene=gs$ENTREZID, OrgDb = org.Mm.eg.db,ont= "ALL",pAdjustMethod = "BH",pvalueCutoff= 1,qvalueCutoff= 1,readable= TRUE,minGSSize=1)    
  kk <- enrichKEGG(gene= gs$ENTREZID, organism = 'mmu',pvalueCutoff = 1,qvalueCutoff=1,minGSSize=1)
  kkk <- kk
}

#!is.null(kk)&&dim(kk)[1]!= 0   ego.all
if (!is.null(ego.all)&&dim(ego.all)[1]!= 0){
  ego.all_pic <- as.data.frame(ego.all)
  ego.all_pic <- subset(ego.all_pic,select=-qvalue)
  ego.all_pic1 <- ego.all_pic
  names(ego.all_pic1) <- c('Category','ID','Description','GeneRatio','BgRatio','pvalue','padj','geneName','Count')

  write.table(ego.all_pic1, file=paste0(outdir,'/','/GO/',cell_type,str_glue('_{ud}_go.csv')),quote=F,sep=',',row.names = FALSE)
  write.table(ego.all_pic1, file=paste0(outdir,'/','/GO/',cell_type,str_glue('_{ud}_go.txt')),quote=F,sep=',',row.names = FALSE)

  if (nrow(subset(ego.all_pic,ONTOLOGY=='BP'))<10){bprows=nrow(subset(ego.all_pic,ONTOLOGY=='BP'))}
  if (nrow(subset(ego.all_pic,ONTOLOGY=='BP'))>=10){bprows=10}
  if (nrow(subset(ego.all_pic,ONTOLOGY=='CC'))<10){ccrows=nrow(subset(ego.all_pic,ONTOLOGY=='CC'))}
  if (nrow(subset(ego.all_pic,ONTOLOGY=='CC'))>=10){ccrows=10}
  if (nrow(subset(ego.all_pic,ONTOLOGY=='MF'))<10){mfrows=nrow(subset(ego.all_pic,ONTOLOGY=='MF'))}
  if (nrow(subset(ego.all_pic,ONTOLOGY=='MF'))>=10){mfrows=10}

  bp <- subset(ego.all_pic,ONTOLOGY=='BP')[1:bprows,c('ONTOLOGY','GeneRatio','Description','p.adjust','Count')]
  cc <- subset(ego.all_pic,ONTOLOGY=='CC')[1:ccrows,c('ONTOLOGY','GeneRatio','Description','p.adjust','Count')]
  mf <- subset(ego.all_pic,ONTOLOGY=='MF')[1:mfrows,c('ONTOLOGY','GeneRatio','Description','p.adjust','Count')]

  ego.all_picc <- rbind(bp,cc,mf)
  ego.all_picc <- dplyr::filter(ego.all_picc,  !is.na(ONTOLOGY))
  ratio <- matrix(as.numeric(unlist(strsplit(as.character(ego.all_picc$GeneRatio),"/"))),ncol=2,byrow=TRUE)
  ego.all_picc$GeneRatio <- ratio[,1]/ratio[,2]

  Description <- c()
    for (i in 1:nrow(ego.all_picc)) {
       if (nchar(as.character(ego.all_picc$Description[i])) >= 50) {
          vectors <- unlist(strsplit(as.character(ego.all_picc$Description[i]),' '))
          description <- paste(vectors[1],vectors[2],vectors[3],vectors[4],'...',sep=' ')
       if (description %in% Description) { description <- paste(description,'.',sep='')}
          Description[i] <- description}
       if (nchar(as.character(ego.all_picc$Description[i])) < 50) {
          Description[i] <- as.character(ego.all_picc$Description[i])}}
  ego.all_picc$Description <- factor(Description,levels=unique(Description))
  ego.all_picc <- ego.all_picc[order(ego.all_picc$p.adjust,decreasing=F),]
  p <- ggplot(ego.all_picc,aes(x=GeneRatio,y=Description,colour=p.adjust,size=Count))
  p <- p + geom_point() + scale_colour_gradientn(colours=rainbow(4),guide="colourbar") + expand_limits(color=seq(0,1,by=0.25)) + theme(axis.text=element_text(color="black",size=10))

  p <- p + theme_bw() + theme(panel.border=element_rect(colour="black"))

  ggsave(paste0(outdir,'/','/GO/',cell_type,str_glue('_{ud}_go.dot.pdf')),p)
  ggsave(paste0(outdir,'/','/GO/',cell_type,str_glue('_{ud}_go.dot.png')),p,type = 'cairo-png')

####################
  p1 <- ggplot(ego.all_picc,aes(x=Description,y=-log10(p.adjust),fill=ONTOLOGY))
  p1 <- p1 + geom_bar(stat='identity') + theme(plot.margin=unit(c(1,1,2,4),'lines'))
  p1 <- p1 + theme(panel.background=element_rect(fill="transparent"),axis.line=element_line())
  p1 <- p1 + theme(axis.text.x=element_text(hjust=1,angle=45,size=6))

  ggsave(paste0(outdir,'/','/GO/',cell_type,str_glue('_{ud}_go.bar.pdf')),p1)
  ggsave(paste0(outdir,'/','/GO/',cell_type,str_glue('_{ud}_go.bar.png')),p1,type = 'cairo-png')
}else{
  print(str_glue('GO {ud} Not enriched'))
}

################################


if (!is.null(kk)&&dim(kk)[1]!= 0){
  kk <- as.data.frame(kk)
  kk <- subset(kk,select=-qvalue)

  gene_list <- strsplit(kk$geneID,split='/')
  clustergene <- as.character(gs$ENTREZID)
  genename_vt <- as.character(gs$SYMBOL)
  names(genename_vt)<-clustergene
  geneName <- unlist(lapply(gene_list,FUN=function(x){paste(genename_vt[x],collapse='/')}))
  if (species=='GRCh38'){
  gss = bitr(gs$ENTREZID, fromType="ENTREZID", toType="ENSEMBL", OrgDb="org.Hs.eg.db")
  }else if (species=='mm10'){
  gss = bitr(gs$ENTREZID, fromType="ENTREZID", toType="ENSEMBL", OrgDb="org.Mm.eg.db")
  }
  clustergene_id <- as.character(gss$ENTREZID)
  geneid_vt <- as.character(gss$ENSEMBL)
  names(geneid_vt)<- clustergene_id
  geneid <- unlist(lapply(gene_list,FUN=function(x){paste(geneid_vt[x],collapse='/')}))

  if (species=='GRCh38'){
    kegg_frame <- read.delim('/public/work/Pipline/Single_RNA/standard_analysis/Transcriptome/08.Enrichment/data/Homo_sapiens/Homo_sapiens_Ensemble_90_hsa_kegg.txt',header=TRUE,sep='\t')
    kegg_frame <- na.omit(kegg_frame)
    kegg_uq <- unique(kegg_frame[,1:2])
    keggid_vt <- as.character(unique(kegg_frame[,1:2])[,2])
    keggid_abbr_vt <- c()
    for ( i in keggid_vt) { 
        keggid_abbr_vt <- c(keggid_abbr_vt,paste('hsa',":",i,sep=""))}

    names(keggid_vt) <- as.character(unique(kegg_frame[,1:2])[,1])
    names(keggid_abbr_vt) <- as.character(unique(kegg_frame[,1:2])[,2])
  }else if (species=='mm10'){
    kegg_frame <- read.delim('/public/work/Pipline/Single_RNA/standard_analysis/Transcriptome/08.Enrichment/data/MM10/Mus_musculus_Ensemble_90_mmu_kegg.txt',header=TRUE,sep='\t')
    kegg_frame <- na.omit(kegg_frame)
    kegg_uq <- unique(kegg_frame[,1:2])
    keggid_vt <- as.character(unique(kegg_frame[,1:2])[,2])
    keggid_abbr_vt <- c()
    for ( i in keggid_vt) { 
        keggid_abbr_vt <- c(keggid_abbr_vt,paste('mmu',":",i,sep=""))}
    names(keggid_vt) <- as.character(unique(kegg_frame[,1:2])[,1])
    names(keggid_abbr_vt) <- as.character(unique(kegg_frame[,1:2])[,2])
  }
  keggID <- unlist(lapply(gene_list,FUN=function(x){paste(keggid_abbr_vt[x],collapse='/')}))
  kk_1  <- data.frame(kk[,1:6],geneid,geneName,keggID,kk[,8,drop=F])
  names(kk_1) <- c('KEGGID','Description','GeneRatio','BgRatio','pvalue','padj','geneID','geneName','keggID','Count')
  write.table(kk_1, file=paste0(outdir,'/','/KEGG/',cell_type,str_glue('_{ud}_kegg.csv')) , sep = ',', row.names = FALSE, quote = FALSE)
  write.table(kk_1, file=paste0(outdir,'/','/KEGG/',cell_type,str_glue('_{ud}_kegg.txt')) , sep = ',', row.names = FALSE, quote = FALSE)
  kk.all_pic <- as.data.frame(kk)
  if (nrow(kk.all_pic)<20) {rows=nrow(kk.all_pic)}
  if (nrow(kk.all_pic)>=20) {rows=20}
  kk.all_pic <- kk.all_pic[1:rows,c('ID','GeneRatio','Description','p.adjust','Count')]
  ratio <- matrix(as.numeric(unlist(strsplit(as.character(kk.all_pic$GeneRatio),"/"))),ncol=2,byrow=TRUE)
  kk.all_pic$GeneRatio <- ratio[,1]/ratio[,2]

  Description <- c()
  for (i in 1:nrow(kk.all_pic)) {
     if (nchar(as.character(kk.all_pic$Description[i])) >= 50) {
        vectors <- unlist(strsplit(as.character(kk.all_pic$Description[i]),' '))
        description <- paste(vectors[1],vectors[2],vectors[3],vectors[4],'...',sep=' ')
        if (description %in% Description) {
            if (paste(description,'.',sep='') %in% Description) {
                description <- paste(description,'..',sep='')
            }}
        Description[i] <- description}
    if (nchar(as.character(kk.all_pic$Description[i])) < 50) {
        Description[i] <- as.character(kk.all_pic$Description[i])}}

    kk.all_pic$Description <- factor(Description,levels=unique(Description))

    kk.all_pic <- kk.all_pic[order(kk.all_pic$p.adjust,decreasing=F),]
    p <- ggplot(kk.all_pic,aes(x=GeneRatio,y=Description,colour=p.adjust,size=Count))
    p <- p + geom_point() +  scale_colour_gradientn(colours=rainbow(4),guide="colourbar") + expand_limits(color=seq(0,1,by=0.25)) +  theme(axis.text=element_text(color="black",size=10))
    p <- p + theme_bw() + theme(panel.border=element_rect(colour="black"))
    ggsave(paste0(outdir,'/','/KEGG/',cell_type,str_glue('_{ud}_kegg.dot.pdf')),p)
    ggsave(paste0(outdir,'/','/KEGG/',cell_type,str_glue('_{ud}_kegg_dot.png')),p,type = 'cairo-png')

    p1 <- ggplot(kk.all_pic,aes(x=Description,y=-log10(p.adjust)))
    p1 <- p1 + geom_bar(stat='identity',fill='steelblue') + theme(plot.margin=unit(c(1,1,2,4),'lines')) + theme(legend.position = "none")
    p1 <- p1 + theme(panel.background=element_rect(fill="transparent"),axis.line=element_line())
    p1 <- p1 + theme(axis.text.x=element_text(hjust=1,angle=45,size=6))
    ggsave(paste0(outdir,'/','/KEGG/',cell_type,str_glue('_{ud}_kegg.bar.pdf')),p1)
    ggsave(paste0(outdir,'/','/KEGG/',cell_type,str_glue('_{ud}_kegg_bar.png')),p1,type = 'cairo-png')


  ### pathview
  pathid = kk_1 %>% filter(nchar(KEGGID) > 7, padj < 0.05) %>% pull(KEGGID)
  dir.create(str_glue("{outdir}/KEGG/pathview_{ud}"), recursive = T)
  setwd(str_glue("{outdir}/KEGG/pathview_{ud}"))
  walk(pathid, function(pathid){
    gene = kk_1 %>% filter(KEGGID == pathid) %>% dplyr::select(geneID) %>% str_split("/",simplify = T) %>% 
      bitr(fromType = "ENSEMBL", toType = "ENTREZID", OrgDb="org.Hs.eg.db") %>% dplyr::select(ENTREZID) %>% 
      inner_join(gs, by = "ENTREZID") %>% inner_join(genelist, by = "SYMBOL") %>% dplyr::select(ENTREZID, logFC) %>% 
      dplyr::distinct(ENTREZID, .keep_all = T) %>% column_to_rownames(var = "ENTREZID")
    pathview(gene.data = gene, pathway.id = pathid, species = "hsa", out.suffix = "gene_fc", kegg.native = T, same.layer=F)
  })
}else{
  print('KEGG up Not enriched')
}
}

#####parser
parser = ArgumentParser()
parser$add_argument("--CellType")
parser$add_argument("--InputPath")
parser$add_argument("--outdir")
# parser$add_argument("--markers1")
# parser$add_argument("--markers2")
# parser$add_argument("--markers3")
parser$add_argument("--species")
parser$add_argument("--SampleNames")
args <- parser$parse_args()
str(args)

cell_type=args$CellType
outdir=args$outdir
InputPath=args$InputPath
species=args$species
SampleNames=args$SampleNames

dir.create(outdir)

# if (!dir.exists(paste0(outdir,'/GO'))){
#   dir.create(paste0(outdir,'/GO'))
# }

# if (!dir.exists(paste0(outdir,'/KEGG'))){
#   dir.create(paste0(outdir,'/KEGG'))
# }

#if (!dir.exists(paste0(outdir,'/KEGG/pathway_plot'))){
#  dir.create(paste0(outdir,'/KEGG/pathway_plot'))
#}

#cell_type='B_cells'
#outdir='/public/work/Project/Single_cell/shouqian_shu_liujinxu/shouhou3_20220604/heart'

#markers<-read.csv("/public/work/Project/Single_cell/shouqian_shu_liujinxu/shouhou3_20220604/heart/B_cells/DIFF/B_cells_diff.xls",header=T)
# markers1<-read.csv(markers1,header=T,sep='\t')
# markers2<-read.csv(markers2,header=T,sep='\t')
# markers3<-read.csv(markers3,header=T,sep='\t')


# markers1_up <- subset(markers1 , cell_ratio_diff > 0.4 )
# up1 <- markers1_up$Gene_name

# markers1_down <- subset(markers1 , cell_ratio_diff < -0.4 )
# down1 <- markers1_down$Gene_name

# markers2_up <- subset(markers2 ,  cell_ratio_diff > 0.4 )
# up2 <- markers2_up$Gene_name

# markers2_down <- subset(markers2 , cell_ratio_diff < -0.4 )
# down2 <- markers2_down$Gene_name

# markers3_up <- subset(markers3 , cell_ratio_diff > 0.4 )
# up3 <- markers3_up$Gene_name

# markers3_down <- subset(markers3 , cell_ratio_diff < -0.4 )
# down3 <- markers3_down$Gene_name
path_list =  InputPath %>% str_split(",") %>% .[[1]]
# print(path_list)

sample_list = SampleNames %>% str_split(",") %>% .[[1]]

walk(c("ratio", "average"), function(x){
  outdir = str_glue("{outdir}/{x}")
  dir.create(outdir, recursive = T)
  if(sample_list %>% length() > 1){

    ob.list_up = map2(sample_list, path_list, function(sample, path){
      a = read.csv(path, header=T, sep='\t')
      up = a %>% filter(!!sym(str_glue("cell_{x}_diff")) > 0.4) %>% dplyr::select(Gene_name, !!sym(str_glue("cell_{x}_diff"))) %>% 
        dplyr::rename(SYMBOL = Gene_name, logFC = !!sym(str_glue("cell_{x}_diff")))
      })
    names(ob.list_up) = sample_list
    # ob.list_up <- list(
    # ob.list_up[['T_MvsTHP1']] <- up1
    # ob.list_up[['R-S1-RAW-coA3LvsR-S2-RAW']] <- up2
    # ob.list_up[['Macrophage_Monocyte_cells_Treat_Model']] <- up3
    
    gene_up_list = lapply(sample_list, function(x) x = ob.list_up[[x]]$SYMBOL)
    names(gene_up_list) = sample_list
    venn_up <- ggvenn(gene_up_list, sample_list,set_name_color=rep("black", length(sample_list))) 
    #+ annotate("text", x = -0.8 , y = 1.2,label = name1,colour="black") +   annotate("text", x = 0.8 , y = -1.2,label = name2,colour="black")

    ggsave(paste0(outdir,'/up_venn.pdf'), venn_up ,height=8,width =8)
    ggsave(paste0(outdir,'/up_venn.png'), venn_up ,height=8,width =8,type = 'cairo-png')

    #########################################

    # ob.list_down <- list()
    # ob.list_down[['T_MvsTHP1']] <- down1
    # ob.list_down[['R-S1-RAW-coA3LvsR-S2-RAW']] <- down2
    # ob.list_down[['Macrophage_Monocyte_cells_Treat_Model']] <- down3
    path_list =  InputPath %>% str_split(",") %>% .[[1]]

    sample_list = SampleNames %>% str_split(",") %>% .[[1]]

    ob.list_down = map2(sample_list, path_list, function(sample, path){
      a = read.csv(path, header=T, sep='\t')
      up = a %>% filter(!!sym(str_glue("cell_{x}_diff")) < -0.4) %>% dplyr::select(Gene_name, !!sym(str_glue("cell_{x}_diff"))) %>% 
        dplyr::rename(SYMBOL = Gene_name, logFC = !!sym(str_glue("cell_{x}_diff")))
      })
    names(ob.list_down) = sample_list
    gene_down_list = lapply(sample_list, function(x) x = ob.list_down[[x]]$SYMBOL)
    names(gene_down_list) = sample_list

    venn_down <- ggvenn(gene_down_list,sample_list,set_name_color=rep("black", length(sample_list))) 
    #+ annotate("text", x = -0.8 , y = 1.2,label = name1,colour="black") +   annotate("text", x = 0.8 , y = -1.2,label = name2,colour="black")

    ggsave(paste0(outdir,'/down_venn.pdf'), venn_down ,height=8,width =8)
    ggsave(paste0(outdir,'/down_venn.png'), venn_down ,height=8,width =8,type = 'cairo-png')


    #up <- intersect(up1,up2)
    #down <- intersect(down1,down2)

    # up1_2 <- intersect(up1,up2)
    # up <- intersect(up1_2,up3)

    # down1_2 <- intersect(down1,down2)
    # down <- intersect(down1_2,down3)
    up = purrr::reduce(gene_up_list, intersect)

    down = purrr::reduce(gene_down_list, intersect)

    if (!is.null(up)){
      write.table(up,file = paste0(outdir,paste0('/','up_gene.xls')),sep='\t',quote = F,row.names =F,col.names='Gene_name')
    }

    if (!is.null(down)){
      write.table(down,file = paste0(outdir,paste0('/','down_gene.xls')),sep='\t',quote = F,row.names =F,col.names='Gene_name')
    }

    #markers_up <- subset(markers , cell_ratio_diff > 0.4 )
    #up <- markers_up$Gene_name

    #markers_down <- subset(markers , cell_ratio_diff < -0.4 )
    #down <- markers_down$Gene_name

    ####所有sample交集差异基因富集
    outdir_all = file.path(outdir, "sample_all")

    dir.create(file.path(outdir_all, "KEGG"),recursive = T)

    dir.create(file.path(outdir_all, "GO"),recursive = T)

    ###up
    up_fc = data.frame(SYMBOL = up, logFC = 1)
    enrich_function(genelist = up_fc, ud = "up", outdir = outdir_all, species = species, cell_type = cell_type)

    ###down
    down_fc = data.frame(SYMBOL = down, logFC = -1)
    enrich_function(genelist = down_fc, ud = "down", outdir = outdir_all, species = species, cell_type = cell_type)
  }

  ####each sample enrich
  walk(sample_list, function(sample){
    ###outdir
    outdir = file.path(outdir, sample)

    dir.create(file.path(outdir, "KEGG"), recursive = T)

    dir.create(file.path(outdir, "GO"), recursive = T)

    ###up
    enrich_function(genelist = ob.list_up[[sample]], ud = "up", outdir = outdir, species = species, cell_type = cell_type)

    ###down
    enrich_function(genelist = ob.list_down[[sample]], ud = "down", outdir = outdir, species = species, cell_type = cell_type)
  })
})
  






#data(korg)

#foldchange <- as.data.frame(markers$avg_log2FC)
#foldchange$SYMBOL <- markers$Gene_name
#colnames(foldchange) <- c('avg_log2FC','SYMBOL')

#if (dim(kkk)[1]== 0){
#  kk <- as.data.frame(kk)
#  sum <- kk
#}else if (dim(kk)[1]== 0){
#  kkk <- as.data.frame(kkk)
#  kkk <- subset(kkk,select=-qvalue)
#  sum <- kkk
#}else{
#  kk <- as.data.frame(kk)
#  kkk <- as.data.frame(kkk)
#  kkk <- subset(kkk,select=-qvalue)
#  sum <- rbind(kkk,kk)
#}

#ID <- sum$ID
#gene_sum <- c(up,down)

#number <- length(ID)
#setwd(paste0(outdir,'/','/KEGG/','pathway_plot'))

#if (species=='GRCh38'){
#  gss = bitr(gene_sum, fromType="SYMBOL", toType="ENSEMBL", OrgDb="org.Hs.eg.db")
#  aa  <-  inner_join(foldchange , gss, by = "SYMBOL")
#  aaa <- aa$avg_log2FC
#  names(aaa) <- aa$ENSEMBL
#  for (i in 1:number){
#    pathway_id=strsplit(ID[i],'hsa')[[1]][2]
#    pathview(gene.data = gene_sum , pathway.id = pathway_id, species = "hsa",out.suffix = pathway_id, kegg.native = T)
#}
#}else if (species=='mm10'){
#  gss = bitr(gene_sum, fromType="SYMBOL", toType="ENSEMBL", OrgDb="org.Mm.eg.db")  
#  aa  <-  inner_join(foldchange , gss, by = "SYMBOL")
#  aaa <- aa$avg_log2FC
#  names(aaa) <- aa$ENSEMBL
#  for (i in 1:number){
#    pathway_id=strsplit(ID[i],'mmu')[[1]][2]
#    pathview(gene.data = gene_sum , pathway.id = pathway_id, species = "mmu",out.suffix = pathway_id, kegg.native = T)
#}
#} 

