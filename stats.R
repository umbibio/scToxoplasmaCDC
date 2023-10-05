
library(tidyr)
#stats <- read.xlsx('../Output/toxo_cdc/ME49_59/tables/all_genes_cyclic_timing.xlsx')
stats <- readRDS('../Input_sub/toxo_cdc/rds_ME49_59/all_genes_cyclic_timing.rds')
Intra.markers.sig <- readRDS('../Input_sub/toxo_cdc/rds_ME49_59/Intra_markers_sig.rds')
trans.markers.sig <- readRDS("../Input_sub/toxo_cdc/rds_ME49_59/rna_markers_rna_trns_sig_v2.rds")


stats.expressed <- stats %>% dplyr::filter(rna.expressed == 1)
nrow(stats.expressed) # expressed
stats.cyclic.rna <- stats %>% dplyr::filter(rna.expressed == 1, rna.cyclic == 1) 
nrow(stats.cyclic.rna) # of cyclic-expr
stats.cyclic.atac <- stats %>% dplyr::filter(rna.expressed == 1, atac.cyclic == 1) 
nrow(stats.cyclic.atac) # cyclic-atac
stats.cyclic.both <- stats %>% dplyr::filter(rna.expressed == 1, rna.cyclic == 1 & atac.cyclic == 1)
nrow(stats.cyclic.both) # cyclic - atac,expr



## DEGs - inferred phase
all(Intra.markers.sig$gene %in% stats$GeneID[stats$rna.cyclic == 1])
length(unique(Intra.markers.sig$gene)) 

Intra.markers.sig.sum <- Intra.markers.sig %>% group_by(cluster) %>% summarise(num.DEGs = n())
Intra.markers.sig.sum

# DEGs - inferred transition points (rna)
length(unique(trans.markers.sig$gene))
trans.markers.sig.sum <- trans.markers.sig %>% group_by(cluster) %>% summarise(num.DEGs = n())
trans.markers.sig.sum

## atac peaks
atac.peaks <- readRDS("../Input_sub/toxo_cdc/rds_ME49_59/atac_peaks.rds")
nrow(atac.peaks)

## cut and run 
cutRun.peaks <- readRDS("../Input_sub/toxo_cdc/rds_ME49_59/CutRunPeaks.rds")
nrow(cutRun.peaks)

cutRun.peaks.genes <- readRDS("../Input_sub/toxo_cdc/rds_ME49_59/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v8.rds")
cutRun.peaks.genes.all <- cutRun.peaks.genes %>% filter(assigned_to_CutRun_peaks == "yes") %>% distinct(TGME49)
nrow(cutRun.peaks.genes.all)

## cyclic AP2s
## According to Korosh periodic splines (smoothness not identified!!) AP2IV-2 is also cyclis
## However when I used 1.1 smoothed splines and run cyclic_genes_FFT.r only 31 AP2s are found to be cyclic.
## AP2IV-2 is not among those and AP2IX-10 (TGGT1_215895) which is cyclic is missed from the figure in paper.

stats.cyclic.AP2 <- stats %>% dplyr::filter( rna.cyclic == 1 & grepl("AP2 domain", ProductDescription))
length(unique(stats.cyclic.AP2$GeneID))
 


## Kourosh original cyclic file, we decided no to use this to report numbers
## this file gives us 33 cyclic AP2s
stats <- read.xlsx('../OutPut/toxo_cdc/ME49_59/tables/all_genes_cyclic_timing_KZ.xlsx')
"test"
