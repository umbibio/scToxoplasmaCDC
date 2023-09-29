library(tidyverse)
library(bedtoolsr)
library(openxlsx)
library(grid)
library(matrixStats)
library(tidyverse)
library(tidytext)
library(RColorBrewer)
library(parallel)
library(ComplexHeatmap)
library(circlize)
library(doParallel)
library(edgeR)
library(dtwclust)
library(geomtextpath)
library(bigmemory)
require(gridExtra)
library(grid)
library(ggVennDiagram)
library(ggVennDiagram)


source('./util_funcs.R')

#######################################################
########### Peak Gene Assignment (CUT&RUN) ############
#######################################################


# all cut and run samples old(AP2XII8) and new batch(multiple concentration of antibody)
# prod.desc  <- read.xlsx('../Input/toxo_genomics/genes/ProductDescription_GT1.xlsx')
# TGGT1_ME49 <- read.xlsx('../Input/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')
# prod.desc <- left_join(prod.desc, TGGT1_ME49, by = c("GeneID" = "TGGT1")) %>% na.omit()
# MJ_annot <- read.xlsx("../Input/Toxo_genomics/genes/MJ_annotation.xlsx")
# MJ_annot <- MJ_annot %>% dplyr::select(!Product.Description)
# prod.desc <- left_join(prod.desc, MJ_annot, by= "TGME49" )
# 


get_peak_genes_assign <- function(gtf, peaks, qval = qval){
  
  # sort cut and run peaks 
  CutRun <- peaks %>% filter(V9 > -log10(qval))
  CutRun <- CutRun %>% dplyr::select(V1, V2, V3, V7, V9)
  peaks.all.sort <- CutRun %>% arrange(V1, as.numeric(V2), as.numeric(V3))
  peaks.all.sort$V4 <- paste(paste(peaks.all.sort$V1, peaks.all.sort$V2, sep = ":"),peaks.all.sort$V3 ,sep = "-" )
  
  
  
  # prep gtf file
  gtf.filt <- gtf %>% dplyr::filter(!grepl('KE.*',gtf$V1))
  ## Remove the first Exon from transcripts.
  gtf.exon <- gtf.filt %>% dplyr::filter(V3 == 'exon')
  gtf.exon.sort <- gtf.exon %>% arrange(V1, V4, V5)
  parse.str <- strsplit(gtf.exon$V9, split = ' ')
  inds <- unlist(lapply(parse.str , function(x) which(grepl("gene_id", x)) + 1))
  gtf.exon$gene_name <- gsub(";", "", unlist(lapply(1:length(inds), function(i) parse.str[[i]][inds[[i]]])))
  gtf.exon$gene_name <- gsub("\"", "", gtf.exon$gene_name)
  gtf.exon <- gtf.exon %>% group_by(V9) %>% mutate(exon.ord = ifelse(V7 == '+', 1:n(), seq(n(), 1, by = -1)),
                                                   multiple.exon = ifelse(n() > 1, T, F))
  ## Remove the exon1, but keep the Intron 1 , build exon2ton
  gtf.exon.2Ton <- gtf.exon %>% mutate(V10 = ifelse(multiple.exon & V7 == '-', min(V4), 
                                                    ifelse(multiple.exon & V7 == '+', min(V5), V4)),
                                       V11 = ifelse(multiple.exon & V7 == '-', max(V4), 
                                                    ifelse(multiple.exon & V7 == '+', max(V5), V5))) %>%
    mutate(V4 = V10, V5 = V11) %>% 
    dplyr::select(-c(exon.ord,multiple.exon, V10, V11) ) %>% distinct()
  
  
  # peak-gene assignement
  
  ## Overlap with peaks and filter peaks that are entirely within the genes.
  ## Overlapping peaks with exon2Ton and check to see if it is entirely within the gene. 
  ## Then from the sorted peaks we throw out all peaks entirely within the gene & not overlapping with exon1, 
  ## These peaks should not be assigned to any peaks.
  
  options(bedtools.path = "/Users/kourosh.zarringhalam/miniconda3/bin/")
  
  peak.genes.ovrlp <- bedtoolsr::bt.intersect(a = peaks.all.sort, b = gtf.exon.2Ton, wo = T)
  peak.genes.filt <- peak.genes.ovrlp %>% dplyr::filter(V10  <= V2 & V11 >= V3)
  peak.filt <- peaks.all.sort[!(peaks.all.sort$V4 %in%  peak.genes.filt$V6), ]
  peak.filt.sort <- peak.filt %>% arrange(V1, as.numeric(V2), as.numeric(V3))
  peak.filt.sort <- peak.filt.sort %>% dplyr::select(c(V1, V2, V3, V4, everything()))
  
  
  ## filter gtf for transcripts only to get the coordinates of start and end of gene
  gtf.filt.trn <- gtf.filt %>% filter(V3 == "transcript")
  gtf.filt.trn$gene_name <- gsub("\\;.*", "", gsub("transcript_id ", "", gsub("-t.*", "", gtf.filt.trn$V9)))
  gtf.filt.trn$gene_name <- gsub("\"", "", gtf.filt.trn$gene_name)
  
  
  ## Filter for first exon coordinates (exon1 coordinates)
  tmp.neg <- gtf.exon %>% filter(V7 == "-") %>% group_by(V9) %>%  dplyr::slice(which.max(V5))
  tmp.pos <- gtf.exon %>% filter(V7 == "+") %>% group_by(V9) %>%  dplyr::slice(which.min(V5))
  gtf.exon1 <- bind_rows(tmp.pos, tmp.neg)
  gtf.exon1.sort <- gtf.exon1 %>% arrange(V1, V4, V5)
  
  
  ## Assign the peaks to nearest upstream gene (look at 5 closest in case of bi-directional)
  
  peaks.genes.dist <- bedtoolsr::bt.closest(a = peak.filt.sort, b = gtf.exon1.sort, D = "b", k = 5)
  parse.str2 <- strsplit(peaks.genes.dist$V15, split = ';')
  peaks.genes.dist$gene_name  <- unlist(lapply(parse.str2, '[[' , 3))
  peaks.genes.dist.trns <- left_join(peaks.genes.dist, gtf.filt.trn, by = "gene_name")
  
  ## V16 is the distance of the peak to the exon 1 
  ## we need to overcome the issue with the  ones with  dist = 0
  ## on pos strand V3.x (end of peak) should not exceed V5.y (end of transcript/exon_n)
  ## on neg strand V2.x (start of peak) is not less than V4.y (beggining of the transcript/exon_1)
  
  peaks.genes.dist.trns <- peaks.genes.dist.trns %>% filter(!(V18 == 0 & V13 == "+" & V3.x > V5.y))
  
  peaks.genes.dist.trns <- peaks.genes.dist.trns %>% filter(!(V18 == 0 & V13 == "-" & V2.x < V4.y))
  
  ## V16 <= 0 means the peak is at upstream 
  ## Find closest gene among top 5 that is upstreaam (min V16)
  peaks.genes.dist.trns <- peaks.genes.dist.trns %>% dplyr::filter(V18 <= 0)
  peaks.genes.dist.trns <- peaks.genes.dist.trns %>% group_by(V4.x) %>% 
    mutate(V19 = V18[which.min(abs(V18))])
  
  
  ## Filter the rest
  peaks.genes.dist.trns <- peaks.genes.dist.trns %>% dplyr::filter(V18 == V19)
  
  ## filter the ones that are too far (2000 bp)
  peaks.genes.dist.trns.filt <- peaks.genes.dist.trns %>% dplyr::filter(abs(V18) < 2000)
  
  
  # merge multiple peaks assigned to a single gene
  # the duplicated peaks are the bidirectioonal peaks 
  
  peak.genes <- peaks.genes.dist.trns.filt
  peak.genes <- peak.genes %>% dplyr::select(V1.x, V2.x, V3.x, V13, gene_name,  V5.x, V6.x) 
  peak.genes.bed.merged <- peak.genes %>% arrange(V2.x) %>% 
    group_by(gene_name) %>% mutate(start_peak = V2.x[which.min(V2.x)], end_peak = V3.x[which.max(V3.x)])  %>% 
    mutate(V4 = ".", V5 = ".")
  
  peak.genes.bed.merged.bed <- peak.genes.bed.merged %>% dplyr::select(V1.x, start_peak, end_peak, V4, V5, V13, gene_name, V5.x, V6.x) %>%
    distinct(gene_name, .keep_all = T)
  
  colnames(prod.desc) <- gsub("GeneID", "TGME49", colnames(prod.desc))
  
  peak.genes.bed.merged.bed <- left_join(peak.genes.bed.merged.bed, prod.desc, by = c("gene_name" = "TGME49"))
  peak.genes.bed.merged <- left_join(peak.genes.bed.merged, prod.desc, by = c("gene_name" = "TGME49"))
  
  # only peaks iinformation to be loaded into IGV
  peak.merged.bed <- peak.genes.bed.merged.bed %>% 
    ungroup() %>% dplyr::select(V1.x,start_peak, end_peak)
  
  
  return(list(peak.gene.all = peaks.genes.dist.trns.filt, 
              peak.gene.merged = peak.genes.bed.merged, 
              peak.gene.merged.bed = peak.genes.bed.merged.bed, 
              peak.gene.merged.IGV =  peak.merged.bed))
  
}


#############################################################################
## narrow peaks from macs2 for all cut and run samples 
#############################################################################


prod.desc  <- read.xlsx('../Input_sub/toxo_genomics/genes/ProductDescription_ME49.xlsx')
ribosomals <- prod.desc[grep('ribosomal', prod.desc$ProductDescription),]
ribosomal.proteins <- prod.desc[grep('ribosomal protein', prod.desc$ProductDescription),]
ribosomal.proteins <- ribosomal.proteins[-grep('KE', ribosomal.proteins$GenomicLocation),]
ribosomal.proteins <- ribosomal.proteins %>% dplyr::select(GeneID, ProductDescription) ## 137 Ribosomal Proteins

in.dir <- "../Input_sub/toxo_cdc/cutNrun/all_macs2_old_new_batch_NO_Filter/"
files <- list.files(in.dir, pattern = ".narrowPeak")
f.names <- gsub("\\.narrowPeak", "", files)

# read all called peaks (new AP2_TY (new) vs 4 controls + MiseqA_2mm (old) vs 4 controls)
qval <- 0.05
all.peaks <- list()
for (f in files) {
  tmp <- read.table(paste(in.dir, f, sep = ''), header=F, sep="\t", quote = NULL)
  tmp <- tmp %>% filter(!grepl("KE.*", V1)) 
  tmp <- tmp %>% filter(V9 > -log10(qval))
  all.peaks <- c(all.peaks, list(tmp))
}
names(all.peaks) <- f.names
all.peaks.tab <- do.call("rbind", all.peaks[1:4])

saveRDS(all.peaks.tab[1:4], "../Input_sub/toxo_cdc/rds_ME49_59/CutRunPeaks.rds")

# peak gene assignment 
gtf.file <- "../Genomes/ToxoDB-59_TgondiiME49.gtf"
gtf <- read.table(gtf.file, header = F, sep = '\t', quote = NULL)
qval <- 0.05

all.peak.genes <- lapply(1:length(files), function(i){
  ff <- f.names[i]
  peaks.df <- all.peaks[[i]]
  tmp <- get_peak_genes_assign(gtf = gtf, peaks = peaks.df, qval = qval)
  tmp <- tmp$peak.gene.merged %>% mutate(data = ff)
  return(tmp)
})
names(all.peak.genes) <- f.names
all.peak.genes.tab <- do.call("rbind", all.peak.genes[1:4])

## intensity 

all.peak.genes.tab.wide <- all.peak.genes.tab %>% dplyr::select(gene_name, V6.x, data, ProductDescription)
all.peak.genes.tab.wide <- all.peak.genes.tab.wide %>% group_by(gene_name, data) %>% mutate(n(), intensity =mean(V6.x))
all.peak.genes.tab.wide <- all.peak.genes.tab.wide %>% dplyr::select(gene_name, data, intensity, ProductDescription) %>% distinct() 

all.peak.genes.tab.wide <- all.peak.genes.tab.wide %>% pivot_wider(names_from = data, values_from = intensity)
all.peak.genes.tab.wide <- all.peak.genes.tab.wide %>% 
  rowwise() %>% mutate(intensity.mean = mean(c_across(where(is.numeric)), na.rm = T))
all.peak.genes.tab.wide <- all.peak.genes.tab.wide %>% 
  rowwise() %>% mutate(intersect = ifelse(sum(is.na(c_across(where(is.numeric)))) >= 1 , "no", "yes"))
all.peak.genes.tab.wide <- left_join(all.peak.genes.tab.wide, ribosomal.proteins, by = c("gene_name" = "GeneID") )
all.peak.genes.tab.wide <- all.peak.genes.tab.wide %>%
  mutate(group = ifelse(is.na(ProductDescription.y), "others", "ribosomal"))
all.peak.genes.tab.wide <- all.peak.genes.tab.wide %>%
  mutate(group2 =  ifelse(str_detect(group, "ribosomal") & intersect == "yes", 'ribo.intersect',
                          ifelse(str_detect(group, "ribosomal") & intersect == "no", 'ribo', 'others')))

write.xlsx(all.peak.genes.tab.wide, "../Input/toxo_cdc/rds_cutRun/cut_run_union_intensity_scores.xlsx")
saveRDS(all.peak.genes.tab.wide, "../Input_sub/toxo_cdc/rds_ME49_59/cut_run_union_intensity_scores.rds")


ggplot(tmp1, aes(x = intensity.mean)) + 
  #geom_histogram(aes(y = ..density..),
  #               colour = 1, fill = "white",binwidth = 1) +
  geom_density(aes(fill = intersect, color = intersect), alpha=.4) +
  theme_bw()

all.peak.genes.tab.wide 
p <- ggplot(all.peak.genes.tab.wide, aes(x = intensity.mean)) + 
  #geom_histogram(aes(y = ..density..),
  #               colour = 1, fill = "white",binwidth = 1) +
  geom_density(aes(fill = intersect, color = intersect), alpha=.4) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 16, face = 'bold', color = 'black'))+
  theme(axis.text.y = element_text(size = 16, face = 'bold', color = 'black'))+
  theme(axis.title = element_text(size = 20, face = "bold", color = "black")) +
  theme(strip.text = element_text(size = 14, face = "bold", colour = "black"))
p



all.peak.genes.tab.wide.sig <- all.peak.genes.tab.wide %>% filter(group == "ribosomal")
p <- ggplot(all.peak.genes.tab.wide.sig, aes(x = intensity.mean)) + 
  # geom_histogram(aes(y = ..density..),
  #                colour = 1, fill = "white",binwidth = 1) +
  geom_density(aes(fill = group2, color = group2), alpha=.4) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 16, face = 'bold', color = 'black'))+
  theme(axis.text.y = element_text(size = 16, face = 'bold', color = 'black'))+
  theme(axis.title = element_text(size = 20, face = "bold", color = "black")) +
  theme(strip.text = element_text(size = 14, face = "bold", colour = "black"))
p



all.peak.genes.ribo <- lapply(all.peak.genes, 
                         function(x) filter(x, str_detect(ProductDescription, "ribosomal protein")))


#############################

## concatenating peaks (union of 4 new data sets peaks)
in.dir <- "../Input_sub/toxo_cdc/cutNrun/all_macs2_old_new_batch_NO_Filter/"
files <- list.files(in.dir, pattern = ".narrowPeak")

qval <- 0.05
all.peaks <- list()
for (f in files) {
  tmp <- read.table(paste(in.dir, f, sep = ''), header=F, sep="\t", quote = NULL)
  tmp <- tmp %>% filter(!grepl("KE.*", V1)) 
  tmp <- tmp %>% filter(V9 > -log10(qval))
  all.peaks <- c(all.peaks, list(tmp))
}
names(all.peaks) <- gsub("\\.narrowPeak", "", files)


gtf.file <- "../Genomes/ToxoDB-59_TgondiiME49.gtf"
gtf <- read.table(gtf.file, header = F, sep = '\t', quote = NULL)
qval <- 0.05

# AP2_Ty vs 4 controls
all.peaks <- all.peaks[1:4]
all.peaks <- do.call("rbind", all.peaks)
peak.genes.union <- get_peak_genes_assign(gtf,all.peaks, qval)

saveRDS(peak.genes.union, "../Input_sub/toxo_cdc/rds_ME49_59/Union_all_new_peaks_0.05_qval_NO_frag_filt.rds")

## do motif search under cut and run peaks of 970 genes 
#  write bed/fasta files for motif search

peak.genes.union <- readRDS("../Input_sub/toxo_cdc/rds_ME49_59/Union_all_new_peaks_0.05_qval_NO_frag_filt.rds")

bed.dir <- "../Input_sub/toxo_cdc/cutNrun/all_macs2_old_new_batch_NO_Filter/BAMM_analysis_v2/allPeaksUnion/"
out.bed <- paste(bed.dir, 'peak_genes_union_all', '.bed', sep = "")
write.table(peak.genes.union$peak.gene.merged.bed[,1:6], out.bed,
            sep = "\t", quote = F, row.names = F, col.names = F)

fasta.dir <- "../Input/toxo_cdc/cutNrun/all_macs2_old_new_batch_NO_Filter/BAMM_analysis_v2/allPeaksUnion/"
out.fasta <- paste(fasta.dir,'peak_genes_union_all' , '.fasta', sep = "")
bedtoolsr::bt.getfasta(fi = '../Input/Toxo_genomics/genome/ToxoDB-59_TgondiiME49_Genome.fasta', 
                       bed = peak.genes.union$peak.gene.merged.bed[,1:6], fo = out.fasta, s = T)



###############################
## add motif info to the peak gene table 

peak.genes.union <- readRDS("../Input_sub/toxo_cdc/rds_ME49_59/Union_all_new_peaks_0.05_qval_NO_frag_filt.rds")
peak.genes <- peak.genes.union$peak.gene.merged.bed
peak.genes$assigned_to_CutRun_peaks <- "yes"
peak.genes <- peak.genes %>% dplyr::select(-c(ProductDescription, V5.x, V6.x))

## intensity scores
all.peak.genes.tab.wide <- read_rds("../Input_sub/toxo_cdc/rds_ME49_59/cut_run_union_intensity_scores.rds")
all.peak.genes.tab.wide <- all.peak.genes.tab.wide[,c(1,3,4,5,6,7,8)]
all.peak.genes <- left_join(all.peak.genes.tab.wide, peak.genes , by = "gene_name" )



## motif info output from BAMM motif finder 

in.dir <- '../Input_sub/toxo_cdc/cutNrun/all_macs2_old_new_batch_NO_Filter/BAMM_analysis_v2/allPeaksUnion/peak_genes_union_all_BaMMmotif/' ## 970 (union)
occ.files <- list.files(in.dir)[grep('occurrence', list.files(in.dir))]
in.file <- paste0(in.dir, occ.files)[-c(3,4)]
all.motifs.tab <- lapply(1:length(in.file), function(i){
  mofit <- gsub('peak_genes_union_all_', '', gsub('\\..*', '', occ.files[i]))
  
  occ.tab <- read.table(in.file[i], header = T, sep = '\t')
  locs <- strsplit(gsub('>', '', occ.tab$seq), split = ':')
  chr <- unlist(lapply(locs, `[[`, 1))
  pos_strd <- unlist(lapply(locs, `[[`, 2))
  strd <- gsub('\\)', '', gsub('.*\\(', '', pos_strd))
  tmp <- strsplit(gsub('\\(.*', '', pos_strd), split = '-')
  strt <- as.numeric(lapply(tmp, `[[`, 1))
  stp <- as.numeric(lapply(tmp, `[[`, 2))
  motif.strd <- occ.tab$strand
  tmp <- strsplit(occ.tab$start..end, split = '\\.')
  motif.strt <- as.numeric(lapply(tmp, `[[`, 1))
  motif.stp <- as.numeric(lapply(tmp, `[[`, 3))
  df <- data.frame(id = gsub('>', '', occ.tab$seq), seq.len = occ.tab$length,
                   chr = chr, start_peak = strt, end_peak = stp, strd = strd, 
                   motif.strt = motif.strt, motif.stp = motif.stp, motif.strd = motif.strd, 
                   pattern = occ.tab$pattern, motif = mofit)
  
  df <- df %>% mutate(motif.strt.rel = ifelse(motif.strd == '-', 2 * seq.len - motif.stp, motif.strt),
                      motif.stp.rel = ifelse(motif.strd == '-', 2 * seq.len - motif.strt, motif.stp))
  
  ## new 
  df <- df %>% mutate(motif.strt.abs = start_peak + motif.strt.rel,
                      motif.stp.abs = start_peak + motif.stp.rel)
  
  return(df)
})

all.motifs.tab <- bind_rows(all.motifs.tab)
all.motifs.tab <- all.motifs.tab %>% dplyr::select(id, motif, pattern)


all.peak.genes$id <- paste0(all.peak.genes$V1.x, ":", 
                            all.peak.genes$start_peak, "-", 
                            all.peak.genes$end_peak, "(", all.peak.genes$V13, ")")


all.peak.genes.motif <- left_join(all.peak.genes, all.motifs.tab, by = "id")


#######################################
## DEGs KD vs WT phase based
KD.vs.WT.phase.wide.desc <- readRDS("../Input_sub/toxo_cdc/rds_ME49_59/AP2XII8_KD_vs_WT_markers_sig_phase_based_new_fc_1_5_WIDE.rds" )
KD.vs.WT.phase.wide.desc <- KD.vs.WT.phase.wide.desc %>% dplyr::select(GeneID, regulation) %>% distinct()


cut.run.tab <- full_join(all.peak.genes.motif, KD.vs.WT.phase.wide.desc, 
                         by = c("gene_name" = "GeneID"))

prod.desc  <- read.xlsx('../Input_sub/toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input_sub/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')
prod.desc <- left_join(prod.desc, TGGT1_ME49, by = c("GeneID" = "TGGT1")) %>% na.omit()
MJ_annot <- read.xlsx("../Input_sub/Toxo_genomics/genes/MJ_annotation.xlsx")
MJ_annot <- MJ_annot %>% dplyr::select(!Product.Description)
prod.desc <- left_join(prod.desc, MJ_annot, by= "TGME49" )

prod.desc <- prod.desc %>% dplyr::select(TGME49, ProductDescription, new.name)


cut.run.tab <- left_join(cut.run.tab, prod.desc, by = c("gene_name" = "TGME49"))

colnames(cut.run.tab) <- c("TGME49", "AP2XII-8_Ty_S4_vs_AP2XII-8_IgG1_peaks" , "AP2XII-8_Ty_S4_vs_RH_IgG1_S1_peaks" , 
                           "AP2XII-8_Ty_S4_vs_RH_Negative_S2_peaks", "AP2XII-8_Ty_S4_vs_RH_Ty_S2_peaks", 
                           "intensity.mean" ,   "intersection_CutRun_dataSets", "chr", "start_peak" , "end_peak", 
                           "V4", "V5", "V6","Genomic_location", "assigned_to_CutRun_peaks" , "id" ,"motif", "pattern" ,
                           "KD_vs_WT_phase_based" ,"ProductDescription", "new.name"  )

cut.run.tab <- cut.run.tab %>% mutate(Category = ifelse(str_detect(ProductDescription, "ribosomal"), "ribosomal", "others"))

write.xlsx(cut.run.tab, "../Input_sub/toxo_cdc/rds_ME49_59/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v8.xlsx")
write.xlsx(cut.run.tab, "../OutPut/toxo_cdc/ME49_59/tables/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v8.xlsx")
saveRDS(cut.run.tab, "../Input_sub/toxo_cdc/rds_ME49_59/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v8.rds")

## modify supplementary table 
cut.run.tab <- readRDS("../Input_sub/toxo_cdc/rds_ME49_59/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v8.rds")
cut.run.tab <- cut.run.tab %>% mutate(MOTIF = case_when(motif == "motif_1" ~ "motif_2", 
                                                        motif == "motif_2" ~ "motif_1"))
cut.run.tab.supplmnt <- cut.run.tab %>% 
  dplyr::select(id, chr, start_peak, end_peak, V4, V5, V6,TGME49, assigned_to_CutRun_peaks,
         intersection_CutRun_dataSets, MOTIF, KD_vs_WT_phase_based, ProductDescription, new.name, Category)
write.xlsx(cut.run.tab.supplmnt, "../OutPut/toxo_cdc/ME49_59/tables/Supplement/cut_run_peak_gene_assignement_supplement.xlsx")
## motif table summary

tab <- read.xlsx("../Input_sub/toxo_cdc/rds_ME49_59/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v8.xlsx")
tab.down <- tab %>% 
  dplyr::select(intersection_CutRun_dataSets, TGME49,  motif, KD_vs_WT_phase_based, ProductDescription, Category) %>%
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

## number of genes with both motifs / 2



#####################

## atac cut&run peaks overlap
Tg_ATAC <- readRDS('../Input/toxo_cdc/rds_ME49_59/S.O_ATAC_peak.rds')

peak.regions <- data.frame(ATAC_peak_region = rownames(Tg_ATAC@assays$peaks@data))
peak.regions <- peak.regions %>% filter(!grepl("KE.*", ATAC_peak_region))
peak.regions.bed  <- data.frame(do.call(rbind, strsplit(peak.regions$ATAC_peak_region,"-")))
peak.regions.bed$X1 <- paste(peak.regions.bed$X1, peak.regions.bed$X2, sep = "_")
peaks.all.sort.atac <- peak.regions.bed  %>% dplyr::select(X1, X3, X4) %>%  arrange(X1, as.numeric(X3), as.numeric(X4))
names(peaks.all.sort.atac) <- c("V1", "V2", "V3")
peaks.all.sort.atac$V4 <- paste(paste(peaks.all.sort.atac$V1, peaks.all.sort.atac$V2, sep = ":"),peaks.all.sort.atac$V3 ,sep = "-" )
peak.genes.bed.merged.bed <- peaks.all.sort.atac

saveRDS(peak.genes.bed.merged.bed,"../Input/toxo_cdc/rds_ME49_59/atac_peaks.rds")


# overlap atac and cut RUN
options(bedtools.path = "/Users/kourosh.zarringhalam/miniconda3/bin/")
peak.cutRun.atac.ovrlp <- bedtoolsr::bt.intersect(a = cut.run.peaks, b = peak.genes.bed.merged.bed, wb = T)
saveRDS(peak.cutRun.atac.ovrlp, "../Input/toxo_cdc/rds_ME49_59/atac_cut_run_overlap.rds")

atac.peaks <- readRDS("../Input/toxo_cdc/rds_ME49_59/atac_peaks.rds")
cut.run.peaks <- readRDS("../Input/toxo_cdc/rds_ME49_59/CutRunPeaks.rds")
peak.cutRun.atac.ovrlp <- readRDS("../Input/toxo_cdc/rds_ME49_59/atac_cut_run_overlap.rds")
venn.plot <- draw.pairwise.venn(
  area1 = nrow(cut.run.peaks),
  area2 = nrow(atac.peaks),
  cross.area = nrow(peak.cutRun.atac.ovrlp),
  #category = c("ATAC", "C&R"),
  fill = c("#469C2C","#F19F39"),
  lty = rep("solid", 2),
  lwd = 6,
  col = c("darkgreen", "darkorange"),
  cex = 5.5,
  cat.cex = 3,
  ext.length = 0.9,
  ext.line.lwd = 2.5,
  
)
grid.draw(venn.plot)
dev.off()



## 

## overlap cut&run and atac genes
peak.genes.atac <- read.table("../Input/toxo_scATAC_MJ_ME49_59/peak_gene_assigned_final.bed")
peak.genes.cutRun <- readRDS("../Input/toxo_cdc/rds_ME49_59/Union_all_new_peaks_0.05_qval_NO_frag_filt.rds")
peak.genes.cutRun <- peak.genes.cutRun[1:7] 

ovlp <- inner_join(peak.genes.cutRun, peak.genes.atac, by = c("gene_name" = "V7"))
# venn.list <- list(atac.genes = unique(peak.genes.atac$V7),
#                   cutRun.genes = unique(peak.genes.cutRun$gene_name))
# ggVennDiagram(venn.list)

pdf(file = "../Output/toxo_cdc/ME49_59/figures_paper/cut_run_genes_overlap_atac_genes_venn.pdf",
    width = 12, height = 12)
venn.plot <- draw.pairwise.venn(
  area1 =length(unique(peak.genes.cutRun$gene_name)),
  area2 = length(unique(peak.genes.atac$V7)) ,
  cross.area = nrow(ovlp),
  #category = c("ATAC", "C&R"),
  fill = c("#469C2C","#F19F39"),
  lty = rep("solid", 2),
  lwd = 6,
  col = c("darkgreen", "darkorange"),
  cex = 5.5,
  cat.cex = 3,
  ext.length = 0.9,
  ext.line.lwd = 2.5,
  
)
grid.draw(venn.plot)
dev.off()

