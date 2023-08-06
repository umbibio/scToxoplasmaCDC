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


source('./util_funcs.R')
source('./util_funcs_YR.R')
## IDs
prod.desc  <- read.xlsx('../Input/toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')


## scDATA

rna_sub <- readRDS('../Input/toxo_cdc/rds_ME49_59/S.O_intra_lables_pt.rds')
atac_sub <- readRDS('../Input/toxo_cdc/rds_ME49_59/S.O_intra_atac_lables_pt.rds')


## Splines
sc.rna.spline.fits <- readRDS('../Input/toxo_cdc/rds_ME49_59/sc_rna_spline_fits_all_genes.rds')
sc.atac.spline.fits <- readRDS('../Input/toxo_cdc/rds_ME49_59/sc_atac_spline_fits_all_genes.rds')

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

## This function performs dynamic time warping clustering for both rna and atac profiles
## it gets a table of genes, 
## the table should include the product description (gene ID and gene Name) 
## you can specify the number of clusters you are looking for (cannot be less than 2)

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


## plot the expression and accessibility of genes within each cluster
# facet by rna cluster
plot_rna_atac_trends <- function(sc.rna.sc.atac.joint.long.sub){
  p  <- ggplot(sc.rna.sc.atac.joint.long.sub, aes(x= time,y=normExpr)) +
    geom_path(aes(color = GeneID),alpha = 0.8, size = 0.8)+ 
    theme_bw(base_size = 14) +
    theme(legend.position = "right") +
    ylab('normExpr') + xlab('Time') +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 18, face="bold", colour = "black")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 18, face="bold", colour = "black")) +
    theme(strip.background = element_rect(colour="black", fill="white",size=0.5, linetype="solid")) +
    theme(strip.text = element_text(size = 18, face="bold", angle = 0)) + 
    
    coord_cartesian(xlim = c(0,6.5)) + 
    facet_grid(cluster.RNA ~ data, scales = 'free', space = 'free') +
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

## plot the expression and accessibility of genes within each cluster
#facet by atac cluster
plot_atac_trends <- function(sc.rna.sc.atac.joint.long.sub){
  p  <- ggplot(sc.rna.sc.atac.joint.long.sub, aes(x= time,y=normExpr)) +
    geom_path(aes(color = GeneID),alpha = 0.8, size = 0.8)+ 
    theme_bw(base_size = 14) +
    theme(legend.position = "right") +
    ylab('normExpr') + xlab('Time') +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 18, face="bold", colour = "black")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 18, face="bold", colour = "black")) +
    theme(strip.background = element_rect(colour="black", fill="white",size=0.5, linetype="solid")) +
    theme(strip.text = element_text(size = 18, face="bold", angle = 0)) + 
    
    coord_cartesian(xlim = c(0,6.5)) + 
    facet_grid(cluster.ATAC ~ data, scales = 'free', space = 'free') +
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
##########################################################
# intersection of 4 CUT&RUN + DEG (KD_vs_WT, phase based)
##########################################################
#tab <- readRDS("../Input/toxo_cdc/rds_ME49_59/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v5.rds")
tab <- readRDS("../Input/toxo_cdc/rds_ME49_59/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v6.rds")
tab <- read.xlsx("../OutPut/toxo_cdc/ME49_59/tables/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v8.xlsx")
HC.peaks <- tab %>% 
  filter(intersection_CutRun_dataSets == "yes" & KD_vs_WT_phase_based %in% c('down_reg', 'up_reg', 'modulated') ) %>% 
  dplyr::select(chr, start_peak, end_peak, V4, V5, V6, TGME49, intersection_CutRun_dataSets, KD_vs_WT_phase_based,
                ProductDescription , Category ) %>% 
  distinct()
names(HC.peaks)[9] <- "dir"


### cluster up_reg genes
HC.peaks <- HC.peaks %>% filter(dir == "up_reg")
HC.peaks$gene_name <- gsub("_", "-", HC.peaks$TGME49) 
HC.peaks.clust <- clust.df(HC.peaks, num.clust = 3)
#HC.peaks.clust <- HC.peaks.clust[HC.peaks.clust$data == "scRNA",]
p <- plot_rna_atac_trends(HC.peaks.clust) 
p
ggsave('../Output/toxo_cdc/ME49_59/figures_paper/High_conf_peaks_up_reg_KD_vs_WT_phase_based_3_clust.pdf', 
       plot = p, width = 6, height = 6, dpi = 300)

## 
HC.peaks.clust.df <- HC.peaks.clust %>% dplyr::select(GeneID, Name, cluster.RNA) %>% distinct()
write.xlsx(HC.peaks.clust.df, "../Output/toxo_cdc/ME49_59/tables/High_conf_peaks_up_reg_KD_vs_WT_phase_based_3_clust.xlsx")

###  cluster down_reg genes
#tab <- readRDS("../Input/toxo_cdc/rds_ME49_59/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v5.rds")
tab <- readRDS("../Input/toxo_cdc/rds_ME49_59/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v6.rds")
tab <- read.xlsx("../OutPut/toxo_cdc/ME49_59/tables/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v8.xlsx")
#names(tab)[1] <- "chr"

HC.peaks <- tab %>% 
  filter(intersection_CutRun_dataSets == "yes" & KD_vs_WT_phase_based %in% c('down_reg', 'up_reg', 'modulated') ) %>% 
  dplyr::select(chr, start_peak, end_peak, V4, V5, V6,TGME49,intersection_CutRun_dataSets, 
                KD_vs_WT_phase_based,ProductDescription , Category ) %>% 
  distinct()
names(HC.peaks)[9] <- "dir"

HC.peaks <- HC.peaks %>% filter(dir == "down_reg")
colnames(HC.peaks) <- gsub("ProductDescription.x", "ProductDescription", colnames(HC.peaks))

HC.peaks$gene_name <- gsub("_", "-", HC.peaks$TGME49) 
HC.peaks.clust <- clust.df(HC.peaks, num.clust = 3)
#HC.peaks.clust <- HC.peaks.clust[HC.peaks.clust$data == "scRNA",]
p <- plot_rna_atac_trends(HC.peaks.clust) 
p

ggsave('../Output/toxo_cdc/ME49_59/figures_paper/High_conf_peaks_down_reg_KD_vs_WT_phase_based_3_clust.pdf', 
       plot = p, width = 6, height = 6, dpi = 300)
HC.peaks.clust.df <- HC.peaks.clust %>% dplyr::select(GeneID, Name, cluster.RNA) %>% distinct()
write.xlsx(HC.peaks.clust.df, "../Output/toxo_cdc/ME49_59/tables/High_conf_peaks_down_reg_KD_vs_WT_phase_based_3_clust.xlsx")


### cluster ribosomal genes 
#tab <- readRDS("../Input/toxo_cdc/rds_ME49_59/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v5.rds")
tab <- readRDS("../Input/toxo_cdc/rds_ME49_59/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v6.rds")
tab <- read.xlsx("../OutPut/toxo_cdc/ME49_59/tables/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v8.xlsx")
#names(tab)[1] <- "chr"

HC.peaks <- tab %>% 
  filter(intersection_CutRun_dataSets == "yes" & KD_vs_WT_phase_based %in% c('down_reg', 'up_reg', 'modulated') ) %>% 
  dplyr::select(chr, start_peak, end_peak, V4, V5, V6,TGME49,intersection_CutRun_dataSets, 
                KD_vs_WT_phase_based,ProductDescription , Category ) %>% 
  distinct()
#names(HC.peaks)[9] <- "dir"


HC.peaks <- HC.peaks %>% filter(Category == "ribosomal")
colnames(HC.peaks) <- gsub("ProductDescription.x", "ProductDescription", colnames(HC.peaks))

HC.peaks$gene_name <- gsub("_", "-", HC.peaks$TGME49) 
HC.peaks.clust <- clust.df(HC.peaks, num.clust = 3)
HC.peaks.clust.rna <- HC.peaks.clust[HC.peaks.clust$data == "scRNA",]

HC.peaks.clust.rna$cluster.RNA <- gsub("C 2", "C 1", gsub("C 3", "C 1", HC.peaks.clust.rna$cluster.RNA ))
p <- plot_rna_atac_trends(HC.peaks.clust.rna) 
p


ggsave('../Output/toxo_cdc/ME49_59/figures_paper/High_conf_peaks_WT_ribosomal_clust.pdf', 
       plot = p, width = 6, height = 4, dpi = 300)

HC.peaks.clust.atac <- HC.peaks.clust[HC.peaks.clust$data == "scATAC",]

p <- plot_atac_trends(HC.peaks.clust.atac) 
p
ggsave('../Output/toxo_cdc/ME49_59/figures_paper/High_conf_peaks_KD_vs_WT_ribosomal_3_clust.pdf', 
       plot = p, width = 6, height = 6, dpi = 300)


## summary table 
#tab <- readRDS("../Input/toxo_cdc/rds_ME49_59/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v5.rds")
#tab <- readRDS("../Input/toxo_cdc/rds_ME49_59/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v6.rds")
tab <- read.xlsx("../OutPut/toxo_cdc/ME49_59/tables/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v8.xlsx")
tab.down <- tab %>% 
  dplyr::select(intersection_CutRun_dataSets, TGME49, motif, KD_vs_WT_phase_based, ProductDescription, Category) %>%
  filter(intersection_CutRun_dataSets == "yes" & KD_vs_WT_phase_based =="down_reg") %>% 
  distinct()

# has at least one motif
tab.motif <- tab.down %>% filter(!is.na(motif)) %>% group_by(TGME49) %>% mutate(motif.list = list(unique(motif)))
tab.motif <- tab.motif %>% rowwise() %>%
  mutate(which_motif = ifelse(length(unlist(motif.list)) > 1 , "both" ,"one")) %>% distinct()
ind <- which(tab.motif$which_motif == "one")
tab.motif$which_motif[ind] <- tab.motif$motif[ind]


# has no motif
tab.no.motif <- tab.down %>% filter(is.na(motif)) %>% mutate(motif.list = "none", which_motif = "none")

tab.down.motif <- rbind(tab.motif, tab.no.motif)
table(tab.down.motif$which_motif, tab.down.motif$Category)
## number of both / 2

## up to here 

###########################################################
# intersection of 4 CUT&RUN + DEG (KD_vs_WT, Global)
##########################################################
tab <- readRDS("../Input/toxo_cdc/rds_ME49_59/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v4.rds")
names(tab)[1] <- "chr"

tab <- tab %>% 
  filter(intersection == "yes" & Global_KD_vs_WT %in% c('down_reg', 'up_reg', 'modulated') ) %>% 
  dplyr::select(chr, start_peak, end_peak, V4, V5, V11,gene_name,intersection, Global_KD_vs_WT, ProductDescription.x , Category ) %>% 
  distinct()
names(tab)[9] <- "dir"

direct <- "down_reg"
HC.peaks <- tab %>% filter(dir == direct)

colnames(HC.peaks) <- gsub("ProductDescription.x", "ProductDescription", colnames(HC.peaks))

HC.peaks$gene_name <- gsub("_", "-", HC.peaks$gene_name) 
HC.peaks.clust <- clust.df(HC.peaks, num.clust = 2)
HC.peaks.clust <- HC.peaks.clust[HC.peaks.clust$data == "scRNA",]
p <- plot_rna_atac_trends(HC.peaks.clust) + ggtitle(direct)
p

ggsave('../Output/toxo_cdc/ME49_59/figures_paper/High_conf_peaks_down_reg_Global_KD_vs_WT_2_clust.pdf', 
       plot = p, width = 6, height = 6, dpi = 300)

HC.peaks.clust.df.down <- HC.peaks.clust %>% dplyr::select(GeneID, Name, cluster.RNA) %>% distinct()
write.xlsx(HC.peaks.clust.df.down, "../Output/toxo_cdc/ME49_59/tables/High_conf_peaks_down_reg_Global_KD_vs_WT_2_clust.xlsx")


# up-reg
tab <- readRDS("../Input/toxo_cdc/rds_ME49_59/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v4.rds")
names(tab)[1] <- "chr"

tab <- tab %>% 
  filter(intersection == "yes" & Global_KD_vs_WT %in% c('down_reg', 'up_reg', 'modulated') ) %>% 
  dplyr::select(chr, start_peak, end_peak, V4, V5, V11,gene_name,intersection, Global_KD_vs_WT, ProductDescription.x , Category ) %>% 
  distinct()
names(tab)[9] <- "dir"

direct <- "up_reg"
HC.peaks <- tab %>% filter(dir == direct)

colnames(HC.peaks) <- gsub("ProductDescription.x", "ProductDescription", colnames(HC.peaks))

HC.peaks$gene_name <- gsub("_", "-", HC.peaks$gene_name) 
HC.peaks.clust <- clust.df(HC.peaks, num.clust = 2)
HC.peaks.clust <- HC.peaks.clust[HC.peaks.clust$data == "scRNA",]

p <- plot_rna_atac_trends(HC.peaks.clust) + ggtitle(direct)
p

ggsave('../Output/toxo_cdc/ME49_59/figures_paper/High_conf_peaks_up_reg_Global_KD_vs_WT_2_clust.pdf', 
       plot = p, width = 6, height = 6, dpi = 300)

HC.peaks.clust.df.up <- HC.peaks.clust %>% dplyr::select(GeneID, Name, cluster.RNA) %>% distinct()
write.xlsx(HC.peaks.clust.df, "../Output/toxo_cdc/ME49_59/tables/High_conf_peaks_up_reg_Global_KD_vs_WT_2_clust.xlsx")


up.down.global <- rbind(HC.peaks.clust.df.up, HC.peaks.clust.df.down)


## ribosomals (which are only present in down-reg list)
tab <- readRDS("../Input/toxo_cdc/rds_ME49_59/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v4.rds")
names(tab)[1] <- "chr"

tab <- tab %>% 
  filter(intersection == "yes" & Global_KD_vs_WT %in% c('down_reg', 'up_reg', 'modulated') ) %>% 
  dplyr::select(chr, start_peak, end_peak, V4, V5, V11,gene_name,intersection, Global_KD_vs_WT, ProductDescription.x , Category ) %>% 
  distinct()
names(tab)[9] <- "dir"

HC.peaks <- tab %>% filter(Category == "ribosomal")


colnames(HC.peaks) <- gsub("ProductDescription.x", "ProductDescription", colnames(HC.peaks))

HC.peaks$gene_name <- gsub("_", "-", HC.peaks$gene_name) 
HC.peaks.clust <- clust.df(HC.peaks, num.clust = 2)
HC.peaks.clust <- HC.peaks.clust[HC.peaks.clust$data == "scRNA",]
#p <- plot_rna_atac_trends(HC.peaks.clust) + ggtitle(direct)
p <- plot_rna_atac_trends(HC.peaks.clust) + ggtitle("ribosomals")
p



# ribosomals
write.xlsx(HC.peaks.clust, "../Output/toxo_cdc/ME49_59/tables/High_conf_peaks_down_reg_Global_KD_vs_WT_ribosomal_2_clust.xlsx")
ggsave('../Output/toxo_cdc/ME49_59/figures_paper/High_conf_peaks_down_reg_Global_KD_vs_WT_ribosomal_2_clust.pdf',
       plot = p, width = 6, height = 6, dpi = 300)



