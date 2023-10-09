
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

S.O.integrated <- readRDS("../Input_sub/toxo_cdc/rds_ME49_59/S.O.integrated_rna_WT_AP2XII8KD_reference_rna_WT_transferred_lables_from_boot.rds")

prod.desc  <- read.xlsx('../Input_sub/toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input_sub/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')


# Marker analysis - KD vs WT (Global)

Idents(S.O.integrated.rna) <- "orig.ident"
unique(Idents(S.O.integrated.rna))
DefaultAssay(S.O.integrated.rna) <- "RNA"
DEGs.KD.vs.WT <- FindAllMarkers(object = S.O.integrated.rna, only.pos = T, min.pct = 0, logfc.threshold = 0)


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

  
# Marker analysis - KD vs WT (phase based) 

# Generate a table of contrasts

makeMatchedContrasts <- function(S.O.integrated){
  
  objs <- as.character(unique(S.O.integrated@meta.data$phase.spp))
  
  contrasts <- data.frame(ref = objs, dummy = 1) %>% full_join( data.frame(query = objs, dummy = 1), by = 'dummy') %>% 
    mutate(ref.spp = gsub(":.*", "", ref), ref.phase = gsub(".*:", "", ref), 
           query.spp = gsub(":.*", "", query), query.phase = gsub(".*:", "", query))
  my.contrasts <- contrasts %>% dplyr::filter(ref.phase == query.phase & ref.spp != query.spp)
  
  return(my.contrasts)
  
}

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


# Assig up-regulated, down-regulated, modulated

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




