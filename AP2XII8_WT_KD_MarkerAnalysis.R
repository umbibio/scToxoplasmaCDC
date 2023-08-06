
library(openxlsx)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(matrixStats)
library(tidyverse)
library(RColorBrewer)
library(sctransform)
library(glassoFast)
library(igraph)
library(ggraph)
library(graphlayouts)
library(Signac)
library(Seurat)
library(patchwork)
library(hdf5r)
library(GenomeInfoDbData)
library(GenomicRanges)
library(GenomicAlignments)
library(Biostrings)
library(rtracklayer)
library(GenomicFeatures)
library(Seurat)
library(plotly)

source('./util_funcs.R')


# new AP2XII-8 KD
S.O.integrated <- readRDS("../Input_sub/toxo_cdc/rds_ME49_59/S.O.integrated_rna_WT_AP2XII8KD_reference_rna_WT_transferred_lables_from_boot.rds")


## subset RNA
## plotting using integrated object just fo testing
## we do not use integrated wt and KD for plotting 
## dont forget to set the default assay to RNA when performing DEG
Idents(S.O.integrated) <- "new.spp"
S.O.integrated.rna <- subset(S.O.integrated, ident = "scRNA")
Idents(S.O.integrated.rna) <- "phase"
DimPlot(S.O.integrated.rna, split.by = "orig.ident", reduction = "pca",
        cols = c("G1a" = "#b6232a","G1b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee'))

Idents(S.O.integrated.rna ) <- "seurat_clausters_indiv"
DimPlot(S.O.integrated.rna, split.by = "orig.ident", reduction = "pca",
        label = T)

################################################################

##  description and name of genes which we add later to our genes of interest
prod.desc  <- read.xlsx('../Input_sub/toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input_sub/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')
prod.desc <- left_join(prod.desc, TGGT1_ME49, by = c("GeneID" = "TGGT1")) %>% na.omit()
MJ_annot <- read.xlsx("../Input_sub/Toxo_genomics/genes/MJ_annotation.xlsx")
MJ_annot <- MJ_annot %>% dplyr::select(!Product.Description)
prod.desc <- left_join(prod.desc, MJ_annot, by= "TGME49" )
prod.desc <- prod.desc %>% 
  mutate(Category = case_when(str_detect(pattern = "ribosomal" ,ProductDescription) ~ "ribosomal",
                                                       TRUE ~ "others"))
  

#################################################################
######### KD vs WT (ignore phases) - global #####################
#################################################################

Idents(S.O.integrated.rna) <- "orig.ident"
unique(Idents(S.O.integrated.rna))
DefaultAssay(S.O.integrated.rna) <- "RNA"
DEGs.KD.vs.WT <- FindAllMarkers(object = S.O.integrated.rna, only.pos = T, min.pct = 0, logfc.threshold = 0)

log2(1.5)

DEGs.KD.vs.WT.sig <- DEGs.KD.vs.WT %>% dplyr::filter(avg_log2FC > 0.58 & p_val_adj < 0.01)
DEGs.KD.vs.WT.sig$group[DEGs.KD.vs.WT.sig$cluster == "scRNA"] <- "down_reg" # down in KD
DEGs.KD.vs.WT.sig$group[DEGs.KD.vs.WT.sig$cluster == "scRNA_KD"] <- "up_reg" # up in KD
DEGs.KD.vs.WT.sig$comparison  <- "Global_KD_vs_WT" 
DEGs.KD.vs.WT.sig$comp.group <- paste(DEGs.KD.vs.WT.sig$group, DEGs.KD.vs.WT.sig$comparison, sep = "_")
DEGs.KD.vs.WT.sig$GeneID <- gsub("-", "_", DEGs.KD.vs.WT.sig$gene)
DEGs.KD.vs.WT.sig$dir[DEGs.KD.vs.WT.sig$cluster == "scRNA"] <- "activated"
DEGs.KD.vs.WT.sig$dir[DEGs.KD.vs.WT.sig$cluster == "scRNA_KD"] <- "repressed" 

DEGs.KD.vs.WT.sig.desc <- left_join(DEGs.KD.vs.WT.sig, prod.desc, by = c("GeneID" = "TGME49"))
saveRDS(DEGs.KD.vs.WT.sig.desc, "../Input_sub/toxo_cdc/rds_ME49_59/AP2XII8_KD_vs_WT_global_markers_sig_fc_1_5.rds")



DEGs.KD.vs.WT.sig.desc.wide <- DEGs.KD.vs.WT.sig.desc %>% dplyr::select(GeneID, comparison, dir)
DEGs.KD.vs.WT.sig.desc.wide <- DEGs.KD.vs.WT.sig.desc.wide %>% 
  pivot_wider(GeneID, names_from = comparison, values_from = dir)

write.xlsx(DEGs.KD.vs.WT.sig.desc, "../Output/toxo_cdc/ME49_59/tables/AP2XII8_KD_vs_WT_global_markers_sig_fc_1_5.xlsx")
saveRDS(DEGs.KD.vs.WT.sig.desc.wide, "../Input_sub/toxo_cdc/rds_ME49_59/AP2XII8_KD_vs_WT_global_markers_sig_fc_1_5_wide.rds")

## bar plot
DEGs.KD.vs.WT.sig.desc <- readRDS("../Input_sub/toxo_cdc/rds_ME49_59/AP2XII8_KD_vs_WT_global_markers_sig_fc_1_5.rds")

DEGs.KD.vs.WT.sig.sum <- DEGs.KD.vs.WT.sig.desc %>% 
  group_by(cluster, dir,Category) %>% summarise(num.deg = n()) 
DEGs.KD.vs.WT.sig.sum$Category <- factor(DEGs.KD.vs.WT.sig.sum$Category, levels = c("ribosomal", "others"))
DEGs.KD.vs.WT.sig.stat <- DEGs.KD.vs.WT.sig.desc %>% group_by(cluster, dir) %>% summarise(num.deg = n()) 
DEGs.KD.vs.WT.sig.stat$dir <- factor(DEGs.KD.vs.WT.sig.stat$dir, levels = c("activated", "repressed"))

p <- ggplot(DEGs.KD.vs.WT.sig.stat, aes(x=dir, y=num.deg, fill = dir)) +
  geom_bar(stat="identity")+
  geom_text(aes(label=num.deg), vjust=1.5, color="black", size=10, fontface = 'bold')+
  theme_bw()+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"))+
  theme(strip.background=element_rect(fill='white', color = 'black'),
        panel.spacing = unit(1.5, "lines"), 
        strip.text.x=element_text(angle=0, hjust=0.5,vjust=0.5, size = "14",face = 'bold'),
        plot.title = element_text(size=18, face = "bold.italic", color = 'black'),
        axis.title.x = element_text(size=16, face="bold"),
        axis.title.y = element_text(size=16, face="bold"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        strip.text = element_text(face = "bold", size = 12,  angle = 0), 
        strip.placement = "outside") +
  theme(legend.text = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold", size = 14))+
  theme(panel.spacing = unit(1.5, "lines")) 


p <- p + ggtitle("Global KD.vs.WT, fc = 1.5")

p

ggsave("../Output/toxo_cdc/ME49_59/figures_paper/AP2XII8_KD_vs_WT_global_markers_sig_fc_1_5.pdf",
       plot = p, height = 6, width = 6, dpi = 300)



###################################################################
################# KD vs WT (phase based) ##########################
###################################################################

# generate table of contrasts
makeMatchedContrasts <- function(S.O.integrated){
  
  objs <- as.character(unique(S.O.integrated@meta.data$phase.spp))
  
  contrasts <- data.frame(ref = objs, dummy = 1) %>% full_join( data.frame(query = objs, dummy = 1), by = 'dummy') %>% 
    mutate(ref.spp = gsub(":.*", "", ref), ref.phase = gsub(".*:", "", ref), 
           query.spp = gsub(":.*", "", query), query.phase = gsub(".*:", "", query))
  my.contrasts <- contrasts %>% dplyr::filter(ref.phase == query.phase & ref.spp != query.spp)
  
  return(my.contrasts)
  
}

# ignore warnings
contrasts <- makeMatchedContrasts(S.O.integrated.rna)
contrasts.groups <- contrasts %>% group_by(ref) %>% summarise(query = list(query))

contrasts.groups$phase <- gsub('.*:', '', contrasts.groups$ref)
contrasts.groups$ref.spp <- gsub(':.*', '', contrasts.groups$ref)


DefaultAssay(S.O.integrated.rna) <- "RNA"
Idents(S.O.integrated.rna) <- "phase.spp"
matched.DEGs <- mclapply(1:nrow(contrasts.groups), function(i){
  tmp <- FindMarkers(S.O.integrated.rna, ident.1 = contrasts.groups$ref[i],
                     only.pos = TRUE, ident.2 = c(unlist(contrasts.groups$query[i])), 
                     verbose = T, min.pct = 0, logfc.threshold = 0)
  tmp$ref <- contrasts.groups$ref[i]
  tmp$ref.spp <- contrasts.groups$ref.spp[i]
  tmp$phase <- contrasts.groups$phase[i]
  tmp$gene <- rownames(tmp)
  tmp$GeneID <- gsub('-', '_', tmp$gene)
  return(tmp)
})

KD.vs.WT.phase.marker <- bind_rows(matched.DEGs)
KD.vs.WT.phase.marker.sig <- KD.vs.WT.phase.marker %>% dplyr::filter(avg_log2FC > log2(1.5) & p_val_adj < 0.01) 

KD.vs.WT.phase.marker.sig$group[KD.vs.WT.phase.marker.sig$ref.spp == "scRNA"] <- "down_reg"
KD.vs.WT.phase.marker.sig$group[KD.vs.WT.phase.marker.sig$ref.spp == "scRNA_KD"] <- "up_reg" 
KD.vs.WT.phase.marker.sig$comparison  <- "KD_vs_WT_phase_based" 
KD.vs.WT.phase.marker.sig$comp.group <- paste(KD.vs.WT.phase.marker.sig$group, KD.vs.WT.phase.marker.sig$comparison, sep = "_")
KD.vs.WT.phase.marker.sig$GeneID <- gsub("-", "_", KD.vs.WT.phase.marker.sig$gene)
KD.vs.WT.phase.marker.sig$dir[KD.vs.WT.phase.marker.sig$ref.spp == "scRNA"] <- "activated"
KD.vs.WT.phase.marker.sig$dir[KD.vs.WT.phase.marker.sig$ref.spp == "scRNA_KD"] <- "repressed" 


KD.vs.WT.phase.marker.sig.desc <- left_join(KD.vs.WT.phase.marker.sig, prod.desc,
                                            by = c("GeneID" = "TGME49"))
saveRDS(KD.vs.WT.phase.marker.sig.desc, "../Input_sub/toxo_cdc/rds_ME49_59/AP2XII8_KD_vs_WT_markers_sig_phase_based_new_fc_1_5.rds")


## turn data to wide format 
# assig up-regulated, down-regulated, modulated
KD.vs.WT.phase.marker.sig.desc <- readRDS( "../Input_sub/toxo_cdc/rds_ME49_59/AP2XII8_KD_vs_WT_markers_sig_phase_based_new_fc_1_5.rds")
write.xlsx(KD.vs.WT.phase.marker.sig.desc, "../OutPut/toxo_cdc/ME49_59/tables/AP2XII8_KD_vs_WT_markers_sig_phase_based_new_fc_1_5.xlsx")

KD.vs.WT.phase.marker.sig.desc.wide <- KD.vs.WT.phase.marker.sig.desc %>% dplyr::select(GeneID, group, phase)
KD.vs.WT.phase.marker.sig.desc.wide <- KD.vs.WT.phase.marker.sig.desc.wide %>% pivot_wider(names_from =  phase, values_from = group)
KD.vs.WT.phase.marker.sig.desc.wide$comparison <- "KD_vs_WT_phase_based" 



KD.vs.WT.phase.wide <- KD.vs.WT.phase.marker.sig.desc.wide %>% 
  rowwise() %>% 
  mutate(across(everything(), ~na_if(.x, "NA")),
         regulation = case_when(all(c_across(C:S) == "up_reg", na.rm = T) ~ "up_reg",
                                all(c_across(C:S) == "down_reg", na.rm = T) ~ "down_reg",
                                all(c_across(C:S) %in% c("up_reg", "down_reg", NA)) ~ "modulated",
                                TRUE ~ "unknown")) %>% 
  ungroup() 

KD.vs.WT.phase.wide <- KD.vs.WT.phase.wide %>% 
  mutate(
  dir = case_when(regulation == "up_reg" ~ "repressed",
                  regulation == "down_reg" ~ "activated",
                  TRUE ~ "modulated"))


KD.vs.WT.phase.wide.desc <- left_join(KD.vs.WT.phase.wide, prod.desc, by = c("GeneID" = "TGME49"))

## proportion of ribosomal genes 
KD.vs.WT.phase.wide.desc <- KD.vs.WT.phase.wide.desc %>%
  mutate(Category = case_when(str_detect(pattern = "ribosomal" ,ProductDescription) ~ "ribosomal",
                              TRUE ~ "others"))
KD.vs.WT.sum <- KD.vs.WT.phase.wide.desc %>% dplyr::select(GeneID, ProductDescription, Category, dir)  %>%
  group_by(dir,Category) %>% summarise(total = n()) 

saveRDS(KD.vs.WT.phase.wide.desc, "../Input_sub/toxo_cdc/rds_ME49_59/AP2XII8_KD_vs_WT_markers_sig_phase_based_new_fc_1_5_WIDE.rds" )


# global vs phase based
venn.list <- list(global.genes = unique(DEGs.KD.vs.WT.sig.desc$GeneID), 
                  phase.based = unique(KD.vs.WT.phase.wide$GeneID) )
p <- ggVennDiagram(venn.list)
p
ggsave("../Output/toxo_cdc/ME49_59/figures_paper/global_and_phase_based_KD_vs_WT_venn_diag.pdf", 
       plot = p, height = 6, width = 6, dpi = 300)


## generate a bar plot 
KD.vs.WT.phase.marker.stat <- KD.vs.WT.phase.marker.sig.desc %>% group_by(ref.spp,dir, phase) %>% summarise(num.deg = n()) 
KD.vs.WT.phase.marker.stat$phase <- factor(KD.vs.WT.phase.marker.stat$phase, levels = c('G1a',"G1b" ,'S', 'M', 'C'))
KD.vs.WT.phase.marker.stat$dir <- factor(KD.vs.WT.phase.marker.stat$dir, levels = c("activated", "repressed"))
# KD.vs.WT.phase.marker.stat$group <- factor(KD.vs.WT.phase.marker.stat$group, levels = c("up_reg", "down_reg"))


p <- ggplot(KD.vs.WT.phase.marker.stat, aes(x=phase, y=num.deg, fill = phase)) +
  geom_bar(stat="identity")+
  geom_text(aes(label=num.deg), vjust=1.5, color="black", size=6, fontface = 'bold')+
  facet_grid(. ~ dir, scales = "free", space='free',
             labeller=label_wrap_gen(multi_line = TRUE))+
  theme_bw()+
  theme(
    axis.text.x = element_blank(),
    #axis.text.y = element_text(size = 10),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"))+
  theme(strip.background=element_rect(fill='white', color = 'black'),
        panel.spacing = unit(1.5, "lines"), 
        strip.text.x=element_text(angle=0, hjust=0.5,vjust=0.5, size = "14",face = 'bold'),
        plot.title = element_text(size=18, face = "bold.italic", color = 'black'),
        axis.title.x = element_text(size=16, face="bold"),
        axis.title.y = element_text(size=16, face="bold"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        strip.text = element_text(face = "bold", size = 12,  angle = 0), 
        strip.placement = "outside") +
  theme(legend.text = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold", size = 14))+
  theme(panel.spacing = unit(1.5, "lines")) 
#theme(legend.position = "none")

p <- p + ggtitle("phase based, KD.vs.WT, fc = 1.5")

p
ggsave("../Output/toxo_cdc/ME49_59/figures_paper/AP2XII8_sig_markers_KD_vs_WT_phase_based_fc_1_5.pdf", 
       plot = p, height = 6, width = 10, dpi = 300)

######

## marker analysis on buldging population

S.O.rna.KD <- readRDS('../Input_sub/toxo_cdc/rds_ME49_59/S.O.rna.AP2XII8.KD.new_transferred_lables_bootroyed.rds')
Idents(S.O.rna.KD) <- "seurat_clusters"
Idents(S.O.rna.KD) <- "phase"
DimPlot(S.O.rna.KD, label = T, dims = c(1,3)) + ggtitle("scRNA_KD")

Idents(S.O.rna.KD) <- "seurat_clusters"
Idents(S.O.rna.KD) <- "phase"
DimPlot(S.O.rna.KD, reduction = "pca", label = T, dims = c(1,2)) + ggtitle("scRNA_KD")


umap.KD <- FetchData(object = S.O.rna.KD, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "phase", "seurat_clusters"))
umap.KD$phase<- factor(umap.KD$phase)
plot_ly(umap.KD, 
        x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3,
        color =  umap.KD$phase, 
        size = 0.2)

umap.KD$seurat_clusters<- factor(umap.KD$seurat_clusters)
plot_ly(umap.KD, 
        x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3,
        color =  umap.KD$seurat_clusters, 
        size = 0.2)


pca.KD <- FetchData(object = S.O.rna.KD, vars = c("PC_1", "PC_2", "PC_3", "seurat_clusters", "phase"))
pca.KD$seurat_clusters <- factor(pca.KD$seurat_clusters)
plot_ly(pca.KD, 
        x = ~PC_1, y = ~PC_2, z = ~PC_3,
        color =  pca.KD$seurat_clusters, 
        size = 0.4)


## scRNA_KD:S1S2 vs scRNA_KD:S5

S.O.integrated <- readRDS("../Input_sub/toxo_cdc/rds_ME49_59/S.O.integrated_rna_WT_AP2XII8KD_reference_rna_WT_transferred_lables_from_boot.rds")

## subset RNA
Idents(S.O.integrated) <- "new.spp"
S.O.integrated.rna <- subset(S.O.integrated, ident = "scRNA")
DefaultAssay(S.O.integrated.rna) <- "RNA"
S.O.integrated.rna@meta.data$spp.phase.seurat.clust <- paste(S.O.integrated.rna@meta.data$orig.ident, 
                                                             S.O.integrated.rna@meta.data$phase.seurat.indiv,
                                                             sep = ":")
Idents(S.O.integrated.rna) <- "spp.phase.seurat.clust"

KD.S.phase.marker <- FindMarkers(S.O.integrated.rna, ident.1 = c("scRNA_KD:S_1", "scRNA_KD:S_2"), 
                                 ident.2 = "scRNA_KD:S_5",only.pos = FALSE, verbose = T, 
                                 min.pct = 0, logfc.threshold = 0)
KD.S.phase.marker$gene <- rownames(KD.S.phase.marker)
KD.S.phase.marker$GeneID <- gsub('-', '_', KD.S.phase.marker$gene)
KD.S.phase.marker$comparison <- "scRNA_KD:S1S2.vs.scRNA_KD:S5"
KD.S.phase.marker <- KD.S.phase.marker %>% mutate(group = ifelse(avg_log2FC > 0 ,"up_reg", "down_reg"))
KD.S.phase.marker <- KD.S.phase.marker %>% mutate(dir = ifelse(avg_log2FC > 0 ,"repressed", "activated"))

KD.S.phase.marker.sig <- KD.S.phase.marker %>% dplyr::filter(abs(avg_log2FC) > log2(1.5) & p_val_adj < 0.01) 
KD.S.phase.marker.sig <- left_join(KD.S.phase.marker.sig, prod.desc, by = c("GeneID" = "TGME49"))

saveRDS(KD.S.phase.marker.sig, "../Input_sub/toxo_cdc/rds_ME49_59/AP2XII8_KD_S1_S2_phase_vs_KD_S5_up_down_reg_markers.rds")


KD.S.phase.marker.sig.wide <- KD.S.phase.marker.sig %>% dplyr::select(GeneID, group, comparison)
KD.S.phase.marker.sig.wide <- KD.S.phase.marker.sig.wide %>% pivot_wider(names_from = comparison , values_from =group )

## bar plot ref.spp
KD.S.phase.marker.sig.stat <- KD.S.phase.marker.sig %>% group_by(comparison,group) %>% summarise(num.deg = n()) 
KD.S.phase.marker.sig.stat$group <- factor(KD.S.phase.marker.sig.stat$group, levels = c("up_reg", "down_reg"))

p <- ggplot(KD.S.phase.marker.sig.stat, aes(x=group, y=num.deg, fill = group)) +
  geom_bar(stat="identity")+
  geom_text(aes(label=num.deg), vjust=1.5, color="black", size=10, fontface = 'bold') + 
  theme_bw() + theme(
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.ticks = element_blank(),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"))+
  theme(legend.text = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold", size = 14)) + 
  theme(plot.title = element_text(size = 18, face = "bold"))
p <- p  + ggtitle(unique(KD.S.phase.marker.sig.stat$comparison))

p
ggsave("../Output/toxo_cdc/ME49_59/figures_paper/scRNA_KD_S1_S2_vs_scRNA_KD_S5_bar_plot.pdf", 
       plot = p, height = 6, width = 6, dpi = 300)



## cluster 1 & 2 vs everything in KD 

S.O.integrated.rna@meta.data$spp.seurat_clust <- paste(S.O.integrated.rna@meta.data$orig.ident, 
                                                       S.O.integrated.rna@meta.data$seurat_clausters_indiv, sep = ":")
Idents(S.O.integrated.rna) <- "spp.seurat_clust"

spp.clust <- sort(unique(S.O.integrated.rna@meta.data$spp.seurat_clust ))
spp.clust.KD <- spp.clust[grepl("KD", spp.clust)]

KD.S.phase.marker.v2 <- FindMarkers(S.O.integrated.rna, ident.1 = c("scRNA_KD:1", "scRNA_KD:2"),
                                 ident.2 = c(spp.clust.KD[c(1,4:10)]), 
                                 only.pos = FALSE, verbose = T, 
                                 min.pct = 0, logfc.threshold = 0)

KD.S.phase.marker.v2.sig <- KD.S.phase.marker.v2 %>% dplyr::filter(abs(avg_log2FC) > log2(1.5) & p_val_adj < 0.01)
KD.S.phase.marker.v2.sig$gene <- rownames(KD.S.phase.marker.v2.sig)
KD.S.phase.marker.v2.sig$GeneID <- gsub('-', '_', KD.S.phase.marker.v2.sig$gene)
KD.S.phase.marker.v2.sig <- KD.S.phase.marker.v2.sig %>% mutate(group = ifelse(avg_log2FC > 0 ,"up_reg", "down_reg"))
KD.S.phase.marker.v2.sig <- KD.S.phase.marker.v2.sig %>% mutate(dir = ifelse(avg_log2FC > 0 ,"repressed", "activated"))
KD.S.phase.marker.v2.sig$comparison <- "scRNA_KD:S1S2.vs.scRNA_KD:all"

KD.S.phase.marker.v2.sig <- left_join(KD.S.phase.marker.v2.sig, prod.desc, by = c("GeneID" = "TGME49"))

saveRDS(KD.S.phase.marker.v2.sig, "../Input_sub/toxo_cdc/rds_ME49_59/AP2XII8_KD_S1_S2_phase_vs_KD_all_up_down_reg_markers.rds")

KD.S.phase.marker.v2.wide <- KD.S.phase.marker.v2.sig %>% dplyr::select(GeneID, group, comparison)
KD.S.phase.marker.v2.wide <- KD.S.phase.marker.v2.wide %>% pivot_wider(names_from = comparison , values_from =group )


## bar plot 
KD.S.phase.marker.v2.stat <- KD.S.phase.marker.v2.sig %>% group_by(comparison,group) %>% summarise(num.deg = n()) 
KD.S.phase.marker.v2.stat$group <- factor(KD.S.phase.marker.sig.stat$group, levels = c("up_reg", "down_reg"))

p <- ggplot(KD.S.phase.marker.v2.stat, aes(x=group, y=num.deg, fill = group)) +
  geom_bar(stat="identity")+
  geom_text(aes(label=num.deg), vjust=1.5, color="black", size=10, fontface = 'bold') + 
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.ticks = element_blank(),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"))+
  theme(legend.text = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold", size = 14)) + 
  theme(plot.title = element_text(size = 18, face = "bold"))
p <- p  + ggtitle(unique(KD.S.phase.marker.v2.stat$comparison))
p
ggsave("../Output/toxo_cdc/ME49_59/figures_paper/scRNA_KD_S1_S2_vs_scRNA_KD_all_bar_plot.pdf", 
       plot = p, height = 6, width = 6, dpi = 300)


## scRNA_KD S1 and S2 vs scRNA_WT S

S.O.integrated.rna@meta.data$spp.phase.clust <- paste(S.O.integrated.rna@meta.data$orig.ident, 
                                                      S.O.integrated.rna@meta.data$phase.seurat, sep = ":")


new.meta <- S.O.integrated.rna@meta.data %>% 
  mutate(Type = case_when(orig.ident == 'scRNA'  ~ phase.spp,
                          orig.ident == 'scRNA_KD'  ~ spp.phase.clust))

S.O.integrated.rna <- AddMetaData(S.O.integrated.rna, new.meta)

DefaultAssay(S.O.integrated.rna) <- "RNA"
Idents(S.O.integrated.rna) <- "Type"
KD.vs.WT.S.phase <- FindMarkers(S.O.integrated.rna, ident.1 = c("scRNA_KD:S_1", "scRNA_KD:S_2"), 
                                ident.2 = "scRNA:S", only.pos = FALSE, verbose = T)

KD.vs.WT.S.phase$gene <- rownames(KD.vs.WT.S.phase)
KD.vs.WT.S.phase$GeneID <- gsub('-', '_', KD.vs.WT.S.phase$gene)
KD.vs.WT.S.phase$comparison <- "scRNA_KD:S1S2.vs.WT:S"
KD.vs.WT.S.phase <- KD.vs.WT.S.phase %>% mutate(group = ifelse(avg_log2FC > 0 ,"up_reg", "down_reg"))
KD.vs.WT.S.phase.sig <- KD.vs.WT.S.phase %>% dplyr::filter(abs(avg_log2FC) > log2(1.5) & p_val_adj < 0.01) 
KD.vs.WT.S.phase.sig <- left_join(KD.vs.WT.S.phase.sig, prod.desc, by = c("GeneID" = "TGME49"))
saveRDS(KD.vs.WT.S.phase.sig, "../Input_sub/toxo_cdc/rds_ME49_59/AP2XII8_KD_S1_S2_phase_vs_WT_S_up_down_reg_markers.rds")

KD.vs.WT.S.phase.sig.wide <- KD.vs.WT.S.phase.sig %>% dplyr::select(GeneID, group, comparison)
KD.vs.WT.S.phase.sig.wide <- KD.vs.WT.S.phase.sig.wide %>% pivot_wider(names_from = comparison , values_from =group )

## bar plot 
KD.vs.WT.S.phase.stat <- KD.vs.WT.S.phase %>% group_by(comparison, group) %>% summarise(num.deg = n()) 
KD.vs.WT.S.phase.stat$group <- factor(KD.vs.WT.S.phase.stat$group, levels = c("up_reg", "down_reg"))

p <- ggplot(KD.vs.WT.S.phase.stat, aes(x=group, y=num.deg, fill = group)) +
  geom_bar(stat="identity")+
  geom_text(aes(label=num.deg), vjust=1.5, color="black", size=10, fontface = 'bold') + 
  theme_bw() + theme(
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.ticks = element_blank(),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"))+
  theme(legend.text = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold", size = 14)) + 
  theme(plot.title = element_text(size = 18, face = "bold"))

p <- p  + ggtitle(unique(KD.vs.WT.S.phase.stat$comparison))
p
ggsave("../Output/toxo_cdc/ME49_59/figures_paper/scRNA_KD_S1_S2_vs_scRNA_WT_S_bar_plot.pdf", 
       plot = p, height = 6, width = 6, dpi = 300)


## KD:S1 and S2 vs everything in WT
## phase based marker analysis across KD and WT 

DefaultAssay(S.O.integrated.rna) <- "RNA"
Idents(S.O.integrated.rna) <- "Type"

Type <- unique(S.O.integrated.rna@meta.data$Type)
group2 <- Type[grepl("scRNA:",unique(Type) )]

KD.S1S2.vs.WT <- FindMarkers(S.O.integrated.rna, ident.1 = c("scRNA_KD:S_1", "scRNA_KD:S_2"), 
                                ident.2 = group2, only.pos = FALSE, verbose = T)

KD.S1S2.vs.WT$gene <- rownames(KD.S1S2.vs.WT)
KD.S1S2.vs.WT$GeneID <- gsub('-', '_', KD.S1S2.vs.WT$gene)
KD.S1S2.vs.WT$comparison <- "scRNA_KD:S1S2.vs.WT:all"
KD.S1S2.vs.WT <- KD.S1S2.vs.WT %>% mutate(group = ifelse(avg_log2FC > 0 ,"up_reg", "down_reg"))

KD.S1S2.vs.WT.sig <- KD.S1S2.vs.WT %>% dplyr::filter(abs(avg_log2FC) > log2(1.5) & p_val_adj < 0.01) 
KD.S1S2.vs.WT.sig <- left_join(KD.S1S2.vs.WT.sig, prod.desc, by = c("GeneID" = "TGME49"))

saveRDS(KD.S1S2.vs.WT.sig, "../Input_sub/toxo_cdc/rds_ME49_59/AP2XII8_KD_S1_S2_phase_vs_WT_all_phases_up_down_reg_markers.rds")


KD.S1S2.vs.WT.sig.wide <- KD.S1S2.vs.WT.sig %>% dplyr::select(GeneID, group, comparison)
KD.S1S2.vs.WT.sig.wide <- KD.S1S2.vs.WT.sig.wide %>% pivot_wider(names_from = comparison , values_from =group )


## bar plot 
KD.S1S2.vs.WT.sig.stat <- KD.S1S2.vs.WT.sig %>% group_by(comparison,group) %>% summarise(num.deg = n()) 
KD.S1S2.vs.WT.sig.stat$group <- factor(KD.S1S2.vs.WT.sig.stat$group, levels = c("up_reg", "down_reg"))

p <- ggplot(KD.S1S2.vs.WT.sig.stat, aes(x=group, y=num.deg, fill = group)) +
  geom_bar(stat="identity")+
  geom_text(aes(label=num.deg), vjust=1.5, color="black", size=10, fontface = 'bold') + 
  theme_bw() + theme(
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.ticks = element_blank(),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"))+
  theme(legend.text = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold", size = 14)) + 
  theme(plot.title = element_text(size = 18, face = "bold")) +
  ggtitle(unique(KD.S1S2.vs.WT.sig.stat$comparison))
  
p  
  

ggsave("../Output/toxo_cdc/ME49_59/figures_paper/scRNA_KD_S1_S2_vs_scRNA_WT_all_bar_plot.pdf", 
       plot = p, height = 6, width = 6, dpi = 300)



###############################################################
###############################################################
### ############## CUT&RUN and modulated genes ################

## description of genes to be added to the mega tab

prod.desc  <- read.xlsx('../Input_sub/toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input_sub/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')
prod.desc <- left_join(prod.desc, TGGT1_ME49, by = c("GeneID" = "TGGT1")) %>% na.omit()
MJ_annot <- read.xlsx("../Input_sub/Toxo_genomics/genes/MJ_annotation.xlsx")
MJ_annot <- MJ_annot %>% dplyr::select(!Product.Description)
prod.desc <- left_join(prod.desc, MJ_annot, by= "TGME49" )
prod.desc <- prod.desc %>% 
  mutate(Category = case_when(str_detect(pattern = "ribosomal" ,ProductDescription) ~ "ribosomal",
                              TRUE ~ "others"))


## Combine all DEG analysis performed 

DEGs.KD.vs.WT.sig.desc <- readRDS("../Input_sub/toxo_cdc/rds_ME49_59/AP2XII8_KD_vs_WT_global_markers_sig_fc_1_5.rds")
KD.vs.WT.phase.wide <- readRDS("../Input_sub/toxo_cdc/rds_ME49_59/AP2XII8_KD_vs_WT_markers_sig_phase_based_new_fc_1_5_WIDE.rds")
names(KD.vs.WT.phase.wide) <- gsub("regulation", "group", names(KD.vs.WT.phase.wide))
KD.S.phase.marker.sig <- readRDS("../Input_sub/toxo_cdc/rds_ME49_59/AP2XII8_KD_S1_S2_phase_vs_KD_S5_up_down_reg_markers.rds")
KD.S.phase.marker.v2.sig <- readRDS("../Input_sub/toxo_cdc/rds_ME49_59/AP2XII8_KD_S1_S2_phase_vs_KD_all_up_down_reg_markers.rds")
KD.vs.WT.S.phase.sig <- readRDS("../Input_sub/toxo_cdc/rds_ME49_59/AP2XII8_KD_S1_S2_phase_vs_WT_S_up_down_reg_markers.rds")
KD.S1S2.vs.WT.sig.wide <- readRDS("../Input_sub/toxo_cdc/rds_ME49_59/AP2XII8_KD_S1_S2_phase_vs_WT_all_phases_up_down_reg_markers.rds")


all.tab <- list(DEGs.KD.vs.WT.sig.desc, 
                KD.vs.WT.phase.wide,
                KD.S.phase.marker.sig, 
                KD.S.phase.marker.v2.sig,
                KD.vs.WT.S.phase.sig,
                KD.S1S2.vs.WT.sig.wide)

names(all.tab) <- c("global_KD_vs_WT",
                    "phase_based_KD_vs_WT", 
                    "KD_S1_S2_vs_KD_S5",
                    "KD_S1_S2_vs_KD_all", 
                    "KD_S1_S2_vs_WT_S", 
                    "KD_S1_S2_vs_WT_all")

all.tab.list <- lapply(1:length(all.tab), function(i){
  
  df <- all.tab[[i]]
  df <- df %>% dplyr::select(GeneID, comparison, group)
  
})
names(all.tab.list) <- c("global_KD_vs_WT",
                         "phase_based_KD_vs_WT", 
                         "KD_S1_S2_vs_KD_S5",
                         "KD_S1_S2_vs_KD_all", 
                         "KD_S1_S2_vs_WT_S", 
                         "KD_S1_S2_vs_WT_all")
all.tab.df <- do.call("rbind", all.tab.list)
all.tab.df.wide <- all.tab.df %>% pivot_wider(names_from = "comparison", values_from = "group")
all.tab.df.wide.desc <- left_join(all.tab.df.wide, prod.desc, by = c("GeneID" = "TGME49") )

modulated <- all.tab.df.wide.desc
saveRDS(modulated, "../Input_sub/toxo_cdc/rds_ME49_59/AP2XII8_KD_modulated_genes_all_comparisons.rds")
write.xlsx(modulated, "../OutPut/toxo_cdc/ME49_59/tables/AP2XII8_KD_modulated_genes_all_comparisons.xlsx")

# add cut and run
CutRun.peaks.motif.chip <- readRDS("../Input_sub/toxo_cdc/rds_ME49_59/union_peaks_cut_run_motif_NO_frag_filt_motif_chip.rds")
CutRun.modulated <- full_join(CutRun.peaks.motif.chip, modulated, by = c("gene_name" = "GeneID"))


