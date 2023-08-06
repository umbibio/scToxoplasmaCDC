library(Seurat)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(gam)
library(princurve)
library(parallel)
library(tidyverse)
library(MyEllipsefit)
library(sctransform)
library(openxlsx)
library(doParallel)
library(tidytext)
library(ggrepel)
library(dtwclust)
library(geomtextpath)
library(bigmemory)

## This script performes time course clustering using the rna expression of the markers within each transition group

clust.df <- function(tab , num.clust) {
  
  k <- num.clust
  sc.rna <- sc.rna.dtw.wide[,colnames(sc.rna.dtw.wide) %in% tab$gene_name ]
  sc.atac <- sc.atac.dtw.wide[,colnames(sc.atac.dtw.wide) %in% tab$gene_name]
  sc.atac<- sc.atac.dtw.wide[,colnames(sc.atac.dtw.wide) %in% colnames(sc.rna)]
  
  sc.rna.markers.hc_dtw <- dtwClustCurves(sc.rna, nclust = k)
  sc.atac.markers.hc_dtw <- dtwClustCurves(sc.atac, nclust = k)
  
  tab <- tab[tab$gene_name %in% colnames(sc.rna),]
  
  sc.rna.long <- inner_join(sc.rna.mu.scale, tab, by = c('GeneID' = 'gene_name'))
  sc.rna.long <- sc.rna.long %>%
    transmute(time = x, GeneID = GeneID, normExpr = expr, Name = ProductDescription)  %>% distinct()
  
  sc.atac.long <- inner_join(sc.atac.mu.scale, tab, by = c('GeneID' = 'gene_name'))
  sc.atac.long <- sc.atac.long %>%
    transmute(time = x, GeneID = GeneID, normExpr = expr, Name = ProductDescription) %>% distinct()
  
  
  sc.rna.clust.info <- data.frame(GeneID = colnames(sc.rna), cluster = cutree(sc.rna.markers.hc_dtw, k = k))
  sc.atac.clust.info <- data.frame(GeneID = colnames(sc.atac), cluster = cutree(sc.atac.markers.hc_dtw, k = k))
  
  sc.rna.long.clust <- inner_join(sc.rna.long, sc.rna.clust.info, by = 'GeneID')
  sc.atac.long.clust <- inner_join(sc.atac.long, sc.atac.clust.info, by = 'GeneID')
  
  sc.rna.sc.atac.joint <- inner_join(sc.rna.long.clust, sc.atac.long.clust,
                                     by = c("time", "GeneID", "Name"))
  colnames(sc.rna.sc.atac.joint) <- c("time", "GeneID", "scRNA", "Name", "cluster.RNA", "scATAC", "cluster.ATAC")
  
  sc.rna.sc.atac.joint.long <- sc.rna.sc.atac.joint %>%
    pivot_longer(-c('time', "GeneID", "Name", "cluster.RNA", "cluster.ATAC"),
                 names_to = 'data', values_to = 'normExpr')
  
  sc.rna.sc.atac.joint.long$cluster.RNA <- paste('C', sc.rna.sc.atac.joint.long$cluster.RNA)
  sc.rna.sc.atac.joint.long$cluster.ATAC <- paste('C', sc.rna.sc.atac.joint.long$cluster.ATAC)
  
  return(sc.rna.sc.atac.joint.long)
  
}

clust.atac.df <- function(tab, num.clust = num.clust){
  
  tab <- tab
  k <- num.clust
  
  sc.atac <- sc.atac.dtw.wide[,colnames(sc.atac.dtw.wide) %in% tab$gene_name]
  sc.atac.markers.hc_dtw <- dtwClustCurves(sc.atac, nclust = k)
  tab <- tab[tab$gene_name %in% colnames(sc.atac),]
  sc.atac.long <- inner_join(sc.atac.mu.scale, tab, by = c('GeneID' = 'gene_name'))
  sc.atac.long <- sc.atac.long %>%
    transmute(time = x, GeneID = GeneID, normExpr = expr, group = group, trans.cluster.rna) %>% 
    distinct()
  sc.atac.clust.info <- data.frame(GeneID = colnames(sc.atac), cluster = cutree(sc.atac.markers.hc_dtw, k = k))
  sc.atac.long.clust <- inner_join(sc.atac.long, sc.atac.clust.info, by = 'GeneID')
  colnames(sc.atac.long.clust) <- c("time", "GeneID", "normExpr", "group", "trans.cluster.rna", "cluster.ATAC")
  sc.atac.long.clust$cluster.ATAC <- paste('C', sc.atac.long.clust$cluster.ATAC)
  
  return(sc.atac.long.clust)
}

# facet by rna cluster
plot_rna_atac_trends <- function(sc.rna.sc.atac.joint.long.sub){
  p  <- ggplot(sc.rna.sc.atac.joint.long.sub, aes(x= time,y=normExpr)) +
    geom_path(aes(color = GeneID),alpha = 0.8, size = 0.8)+ 
    theme_bw(base_size = 14) +
    theme(legend.position = "right") +
    ylab('normExpr') + xlab('Time') +
    theme(axis.text.x = element_text(angle = 0, hjust = 1,face="bold", size = 15, colour = "black")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, face="bold",size = 15, colour = "black")) +
    theme(strip.background = element_rect(colour="black", fill="white",size=0.5, linetype="solid")) +
    theme(strip.text = element_text(size = 18, face="bold", angle = 0)) + 
    
    coord_cartesian(xlim = c(0,6.5)) + 
    facet_grid(cluster.RNA ~ data , scales = 'free', space = 'free') +
    #facet_grid(data ~ cluster.RNA, scales = 'free', space = 'free') +
    theme(
      plot.title = element_text(size=20, face = "bold.italic", color = 'red'),
      axis.title.x = element_text(size=18, face="bold", hjust = 1),
      axis.title.y = element_text(size=18, face="bold")
    ) + 
    theme(#legend.position = c(0.15, 0.85),
      legend.position = 'none',
      legend.title = element_text(colour="black", size=12, 
                                  face="bold"),
      legend.text = element_text(colour="black", size=12, 
                                 face="bold"))
  
  
  return(p)
  
}


plot_atac_trand <- function(sc.atac.long.clust){
  
  p  <- ggplot(sc.atac.long.clust, aes(x= time,y=normExpr)) +
    geom_path(aes(color = GeneID),alpha = 0.8, size = 0.8)+ 
    theme_bw(base_size = 14) +
    theme(legend.position = "right") +
    ylab('normExpr') + xlab('Time') +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 15, face="bold", colour = "black")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 15, face="bold", colour = "black")) +
    theme(strip.background = element_rect(colour="black", fill="white",size=0.5, linetype="solid")) +
    theme(strip.text = element_text(size = 18, face="bold", angle = 0)) + 
    
    coord_cartesian(xlim = c(0,6.5)) + 
    facet_grid(cluster.ATAC ~ ., scales = 'free', space = 'free') +
    theme(
      plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
      axis.title.x = element_text(size=18, face="bold", hjust = 1),
      axis.title.y = element_text(size=18, face="bold")
    ) + 
    theme(#legend.position = c(0.15, 0.85),
      legend.position = 'none',
      legend.title = element_text(colour="black", size=12, 
                                  face="bold"),
      legend.text = element_text(colour="black", size=12, 
                                 face="bold")) +
    theme(axis.ticks =  element_blank())
  
  return(p)
}

plot_rna_atac_trends.ord <- function(sc.rna.sc.atac.joint.long.sub){
  p  <- ggplot(sc.rna.sc.atac.joint.long.sub, aes(x= time,y=normExpr)) +
    geom_path(aes(color = GeneID),alpha = 0.8, size = 0.8)+ 
    theme_bw() +
    #theme_bw(base_size = 16) +
    theme(legend.position = "right") +
    ylab('normExpr') + xlab('Time') +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 22, face="bold", colour = "black")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 18, face="bold", colour = "black")) +
    theme(strip.background = element_rect(colour="black", fill="white",size=0.5, linetype="solid")) +
    theme(strip.text = element_text(size = 22, face="bold", angle = 0)) + 
    
    coord_cartesian(xlim = c(0,6.5)) + 
    facet_grid(cluster.RNA.ordered ~ data, scales = 'free', space = 'free') +
    theme(
      plot.title = element_text(size=20, face = "bold.italic", color = 'red'),
      axis.title.x = element_text(size=22, face="bold", hjust = 1),
      axis.title.y = element_text(size=22, face="bold")
    ) + 
    theme(#legend.position = c(0.15, 0.85),
      legend.position = 'none',
      legend.title = element_text(colour="black", size=12, 
                                  face="bold"),
      legend.text = element_text(colour="black", size=12, 
                                 face="bold"))
  
  
  return(p)
  
}

source('./util_funcs.R')

## IDs
prod.desc  <- read.xlsx('../Input_sub/toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input_sub/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')

## scDATA

rna_sub <- readRDS('../Input_sub/toxo_cdc/rds_ME49_59/S.O_intra_lables_pt.rds')
atac_sub <- readRDS('../Input_sub/toxo_cdc/rds_ME49_59/S.O_intra_atac_lables_pt.rds')


## Splines
sc.rna.spline.fits <- readRDS('../Input_sub/toxo_cdc/rds_ME49_59/sc_rna_spline_fits_all_genes.rds')
sc.atac.spline.fits <- readRDS('../Input_sub/toxo_cdc/rds_ME49_59/sc_atac_spline_fits_all_genes.rds')

## Turn the data into wide format (time by gene) and center & scale each gene
sc.rna.dtw.wide <- sc.rna.spline.fits %>% 
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  mutate_at(vars(matches("TGME")), ~scale(., center = T, scale = T)) %>%
  as.data.frame()

na.ind <- which(apply(sc.rna.dtw.wide, 2, function(x) any(is.na(x))))
if(length(na.ind)){
  sc.rna.dtw.wide <- sc.rna.dtw.wide[,-na.ind]
}

sc.atac.dtw.wide <- sc.atac.spline.fits %>% 
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  mutate_at(vars(matches("TGME")), ~scale(., center = T, scale = T)) %>%
  as.data.frame()

na.ind <- which(apply(sc.atac.dtw.wide, 2, function(x) any(is.na(x))))
if(length(na.ind)){
  sc.atac.dtw.wide <- sc.atac.dtw.wide[,-na.ind]
}

sc.rna.mu.scale <- sc.rna.dtw.wide %>% 
  pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')

sc.atac.mu.scale <- sc.atac.dtw.wide %>% 
  pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')


## 


rna.trans.marker.genes <- readRDS('../Input_sub/toxo_cdc/rds_ME49_59/rna_markers_rna_trns_sig_v2.rds')
rna.trans.marker.genes <- rna.trans.marker.genes %>% transmute(GeneID = gene, phase = cluster) %>% distinct()
rna.trans.marker.genes %>% group_by(phase) %>% summarise(n())

prod.desc  <- read.xlsx('../Input_sub/toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input_sub/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')
prod.desc <- left_join(prod.desc, TGGT1_ME49, by = c("GeneID" = "TGGT1"))
prod.desc <- prod.desc %>% 
  transmute(GeneID = gsub("_", "-", TGME49), ProductDescription = ProductDescription) %>% na.omit()

rna.trans.marker.genes <- left_join(rna.trans.marker.genes, prod.desc, by = "GeneID")
colnames(rna.trans.marker.genes)[1] <- "gene_name"

rna.trans.marker.genes.list <- split(rna.trans.marker.genes, rna.trans.marker.genes$phase) 

## cluster markers of each rna-transition (T1, T2, T3, T4)
trans.list <- c()
k <- 4

trans.list <- lapply(1:length(rna.trans.marker.genes.list), function(i) {
  my.df <- rna.trans.marker.genes.list[[i]]
  df <- clust.df(my.df, num.clust = k)
  df$group <- names(rna.trans.marker.genes.list)[i]
  df$data <- factor(df$data, levels = c("scRNA", "scATAC"))
  p1 <- plot_rna_atac_trends(df) +
    ggtitle(names(rna.trans.marker.genes.list[i]))
  p1
})

names(trans.list) <- names(rna.trans.marker.genes.list)



## order the rna clusters manually according to their peaks

rna.trans.marker.genes.list <- readRDS("../Input_sub/toxo_cdc/rds_ME49_59/rna_markers_rna_transitions_dtw_clust_list.rds")
rna.trans.data.list <- lapply(rna.trans.marker.genes.list, "[[", 1)
rna.trans.data <- do.call("rbind", rna.trans.data.list)

rna.trans.data <- rna.trans.data %>% 
  mutate(new.ord.clust = 
           case_when(group == "T1" & cluster.RNA == "C 1" ~ "D4",
                     group == "T1" & cluster.RNA == "C 2" ~ "D2", 
                     group == "T1" & cluster.RNA == "C 3" ~ "D3", 
                     group == "T1" & cluster.RNA == "C 4" ~ "D1",
                     group == "T2" & cluster.RNA == "C 1" ~ "D4",
                     group == "T2" & cluster.RNA == "C 2" ~ "D2", 
                     group == "T2" & cluster.RNA == "C 3" ~ "D1", 
                     group == "T2" & cluster.RNA == "C 4" ~ "D3", 
                     group == "T3" & cluster.RNA == "C 1" ~ "D4",
                     group == "T3" & cluster.RNA == "C 2" ~ "D3",
                     group == "T3" & cluster.RNA == "C 3" ~ "D2", 
                     group == "T3" & cluster.RNA == "C 4" ~ "D1", 
                     group == "T4" & cluster.RNA == "C 1" ~ "D3", 
                     group == "T4" & cluster.RNA == "C 2" ~ "D4", 
                     group == "T4" & cluster.RNA == "C 3" ~ "D1", 
                     group == "T4" & cluster.RNA == "C 4" ~ "D2", 
                     TRUE ~ "NA"))
rna.trans.data$cluster.RNA.ordered <- gsub("\\D", "C", rna.trans.data$new.ord.clust) 
rna.trans.data$data <- factor(rna.trans.data$data, levels = c("scRNA", "scATAC"))
rna.trans.data.list <- split(rna.trans.data, f= rna.trans.data$group)

saveRDS(rna.trans.data.list, "../Input_sub/toxo_cdc/rds_ME49_59/rna_markers_rna_transitions_dtw_clust_list_ordered.rds")


rna.trans.data.list <- readRDS("../Input_sub/toxo_cdc/rds_ME49_59/rna_markers_rna_transitions_dtw_clust_list_ordered.rds")
trans.plt <- c()
trans.plt <- lapply(1:length(rna.trans.data.list), function(i) {
  
  my.df <- rna.trans.data.list[[i]]
  
  p1 <- plot_rna_atac_trends.ord(my.df) 
  #   ggtitle(names(rna.trans.marker.genes.list[i]))
  p1
})

pp <- grid.arrange(grobs = trans.plt, ncol = 2)

ggsave( "../Output/toxo_cdc/ME49_59/figures_paper/rna_markers_rna_transitions_dtw_clust_list_ordered.pdf", 
        plot = pp,
        height = 16,width = 16, dpi = 300)


p1 <-trans.plt[[1]]
p1
ggsave( "../Output/toxo_cdc/ME49_59/figures_paper/T1_rna_clust_ordered.pdf",
        plot = p1,height = 8,width = 10, dpi = 300)
p2 <-trans.plt[[2]]
p2
ggsave( "../Output/toxo_cdc/ME49_59/figures_paper/T2_rna_clust_ordered.pdf",
        plot = p2,height = 8,width = 10, dpi = 300)

p3 <-trans.plt[[3]]
p3
ggsave( "../Output/toxo_cdc/ME49_59/figures_paper/T3_rna_clust_ordered.pdf",
        plot = p3,height = 8,width = 10, dpi = 300)

p4 <-trans.plt[[4]]
p4
ggsave( "../Output/toxo_cdc/ME49_59/figures_paper/T4_rna_clust_ordered.pdf",
        plot = p4,height = 8,width = 10, dpi = 300)


