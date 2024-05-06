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
files <- files[-3] # take out RH-Neg control which is from 500 million parasites
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
all.peaks.tab <- do.call("rbind", all.peaks[1:3])

## NEW
all.peaks.tab$CasevsControl <- gsub("_peak.*", "", all.peaks.tab$V4) 
all.peaks.tab$Case <- unlist(lapply(strsplit(all.peaks.tab$CasevsControl, "_vs_"), "[[", 1))
all.peaks.tab$Control <- unlist(lapply(strsplit(all.peaks.tab$CasevsControl, "_vs_"), "[[", 2))
all.peaks.tab <- all.peaks.tab %>% `rownames<-`(NULL)
names(all.peaks.tab)[1:10] <- c('chr', "strt", "end", "peak_name", "V5", "V6", "SignalValue", "-log10(pVal)", "-log10(qVal)", "peakScore")


saveRDS(all.peaks.tab, "../Input_sub/toxo_cdc/cutNrun_revision/rds/CutRunPeaks_3Controls_revision.rds")


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
all.peak.genes.tab <- do.call("rbind", all.peak.genes[1:3])
saveRDS(all.peak.genes.tab, "../Input_sub/toxo_cdc/cutNrun_revision/rds/CutRunPeaks_3Controls_genes_assigned_revision.rds")


#############################
## Ven diagram - overalp of peaks across 3 cut and run data 

all.peak.genes.tab <- readRDS("../Input_sub/toxo_cdc/cutNrun_revision/rds/CutRunPeaks_3Controls_genes_assigned_revision.rds")
all.peak.genes.tab.list <- split(all.peak.genes.tab, f = all.peak.genes.tab$data)
names(all.peak.genes.tab.list)
peaks.venn.list <- list(AP2XII8_vs_AP2XII8_IgG1 = unique(all.peak.genes.tab.list$`AP2XII-8_Ty_S4_vs_AP2XII-8_IgG1_peaks`$gene_name),
                        AP2XII8_vs_RH_IgG1 = unique(all.peak.genes.tab.list$`AP2XII-8_Ty_S4_vs_RH_IgG1_S1_peaks`$gene_name),
                        AP2XII8_vs_RH_Ty = unique(all.peak.genes.tab.list$`AP2XII-8_Ty_S4_vs_RH_Ty_S2_peaks`$gene_name))

p <- ggVennDiagram(peaks.venn.list, label = "both", label_alpha = 0.2,
                   label_color = "white", label_size = 14, set_size = 4)
p
ggsave("../Output/toxo_cdc/ME49_59/figures_paper/RevisionUpdate/all_new_cut_run_peak_gene_assigned_overlap_venn_qval_0.05_NO_filt_3Controls_revision.pdf", 
       plot = p, height = 12, width = 12, dpi = 300)


all.pesks.ribo <- all.peak.genes.tab %>% filter(str_detect(ProductDescription, "ribosomal protein"))
all.pesks.ribo.list <- split(all.pesks.ribo, f= all.pesks.ribo$data)
names(all.pesks.ribo.list)
ribo.peaks.venn.list <- list(AP2XII8_vs_AP2XII8_IgG1 = unique(all.pesks.ribo.list$`AP2XII-8_Ty_S4_vs_AP2XII-8_IgG1_peaks`$gene_name),
                             AP2XII8_vs_RH_IgG1 = unique(all.pesks.ribo.list$`AP2XII-8_Ty_S4_vs_RH_IgG1_S1_peaks`$gene_name),
                             AP2XII8_vs_RH_Ty = unique(all.pesks.ribo.list$`AP2XII-8_Ty_S4_vs_RH_Ty_S2_peaks`$gene_name))
p <- ggVennDiagram(ribo.peaks.venn.list, label = "both", label_alpha = 0.2,
                   label_color = "white", label_size = 14, set_size = 4)
p
ggsave("../Output/toxo_cdc/ME49_59/figures_paper/RevisionUpdate/all_new_cut_run_peak_gene_assigned_overlap_venn_qval_0.05_NO_filt_3Controls_revision_ribosomal.pdf", 
       plot = p, height = 12, width = 12, dpi = 300)

## intensity 

all.peak.genes.tab.wide <- all.peak.genes.tab %>% dplyr::select(gene_name, V5.x, data, ProductDescription)
all.peak.genes.tab.wide <- all.peak.genes.tab.wide %>% group_by(gene_name, data) %>% mutate(n(), intensity =mean(V5.x))
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

write.xlsx(all.peak.genes.tab.wide, "../Input_sub/toxo_cdc/cutNrun_revision/rds/cut_run_union_intensity_scores_3Controls_revision.xlsx")
saveRDS(all.peak.genes.tab.wide, "../Input_sub/toxo_cdc/cutNrun_revision/rds/cut_run_union_intensity_scores_3Controls_revision.rds")

all.peak.genes.tab.wide %>% head()
p <- ggplot(all.peak.genes.tab.wide, aes(x = intensity.mean)) + 
  #geom_histogram(aes(y = ..density..),
  #               colour = 1, fill = "white",binwidth = 1) +
  geom_density(aes(fill = intersect, color = intersect), alpha=.4) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 16, face = 'bold', color = 'black'))+
  theme(axis.text.y = element_text(size = 16, face = 'bold', color = 'black'))+
  theme(axis.title = element_text(size = 20, face = "bold", color = "black")) +
  theme(strip.text = element_text(size = 14, face = "bold", colour = "black"),
        plot.title = element_text(size = 16, face = "bold", color = "black")) +
  ggtitle("all genes")
p

ggsave("../OutPut/toxo_cdc/ME49_59/figures_paper/RevisionUpdate/Intensity_scores_density_plot_all_union_genes.pdf", 
       plot = p, width = 6, height = 6, dpi = 300)

all.peak.genes.tab.wide.sig <- all.peak.genes.tab.wide %>% filter(group == "ribosomal")
p <- ggplot(all.peak.genes.tab.wide.sig, aes(x = intensity.mean)) + 
  # geom_histogram(aes(y = ..density..),
  #                colour = 1, fill = "white",binwidth = 1) +
  geom_density(aes(fill = group2, color = group2), alpha=.4) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 16, face = 'bold', color = 'black'))+
  theme(axis.text.y = element_text(size = 16, face = 'bold', color = 'black'))+
  theme(axis.title = element_text(size = 20, face = "bold", color = "black")) +
  theme(strip.text = element_text(size = 14, face = "bold", colour = "black"),
        plot.title = element_text(size = 16, face = "bold", color = "black"))+
  ggtitle("ribosomal protein genes")
p
ggsave("../OutPut/toxo_cdc/ME49_59/figures_paper/RevisionUpdate/Intensity_scores_density_plot_ribosomal_genes.pdf", 
       width = 6, height = 6, dpi = 300)


all.peak.genes.ribo <- lapply(all.peak.genes, 
                              function(x) filter(x, str_detect(ProductDescription, "ribosomal protein")))



#############################

## concatenating peaks (union of 3 new data sets peaks)
in.dir <- "../Input_sub/toxo_cdc/cutNrun/all_macs2_old_new_batch_NO_Filter/"
files <- list.files(in.dir, pattern = ".narrowPeak")
files <- files[-3]

qval <- 0.05
all.peaks <- list()
for (f in files) {
  tmp <- read.table(paste(in.dir, f, sep = ''), header=F, sep="\t", quote = NULL)
  tmp <- tmp %>% filter(!grepl("KE.*", V1)) 
  tmp <- tmp %>% filter(V9 > -log10(qval))
  all.peaks <- c(all.peaks, list(tmp))
}
names(all.peaks) <- gsub("\\.narrowPeak", "", files)
all.peaks <- do.call("rbind", all.peaks)
## NEW
all.peaks$CasevsControl <- gsub("_peak.*", "", all.peaks$V4) 
all.peaks$Case <- unlist(lapply(strsplit(all.peaks$CasevsControl, "_vs_"), "[[", 1))
all.peaks$Control <- unlist(lapply(strsplit(all.peaks$CasevsControl, "_vs_"), "[[", 2))
all.peaks <- all.peaks %>% `rownames<-`(NULL)
#names(all.peaks)[1:10] <- c('chr', "strt", "end", "peak_name", "V5", "V6", "SignalValue", "-log10(pVal)", "-log10(qVal)", "peakScore")


gtf.file <- "../Genomes/ToxoDB-59_TgondiiME49.gtf"
gtf <- read.table(gtf.file, header = F, sep = '\t', quote = NULL)
qval <- 0.05

# AP2_Ty vs 3 controls

peak.genes.union <- get_peak_genes_assign(gtf,all.peaks, qval)


saveRDS(peak.genes.union, "../Input_sub/toxo_cdc/cutNrun_revision/rds/Union_all_new_peaks_0.05_qval_NO_frag_filt_3Controls.rds")

## IGV peaks
peak.genes.union <- readRDS("../Input_sub/toxo_cdc/cutNrun_revision/rds/Union_all_new_peaks_0.05_qval_NO_frag_filt_3Controls.rds")
write.table(peak.genes.union$peak.gene.merged.IGV, "../Input_sub/toxo_cdc/cutNrun_revision/all_macs2_old_new_batch_NO_Filter/peaks/peak_gene_cutRun_3Controls_final_only_peaks.bed",
            sep = "\t", quote = F, row.names = F, col.names = F)


## do motif search under cut and run peaks of 970 genes - done by KZ
#  write bed/fasta files for motif search

peak.genes.union <- readRDS("../Input_sub/toxo_cdc/cutNrun_revision/rds/Union_all_new_peaks_0.05_qval_NO_frag_filt_3Controls.rds")

bed.dir <- "../Input_sub/toxo_cdc/cutNrun_revision/all_macs2_old_new_batch_NO_Filter/BAMM_analysis_v2/allPeaksUnion/"
out.bed <- paste(bed.dir, 'peak_genes_union_all', '.bed', sep = "")
write.table(peak.genes.union$peak.gene.merged.bed[,1:6], out.bed,
            sep = "\t", quote = F, row.names = F, col.names = F)


fasta.dir <- "../Input_sub/toxo_cdc/cutNrun_revision/all_macs2_old_new_batch_NO_Filter/BAMM_analysis_v2/allPeaksUnion/"
out.fasta <- paste(fasta.dir,'peak_genes_union_all' , '.fasta', sep = "")
bedtoolsr::bt.getfasta(fi = '../Input_sub/Toxo_genomics/genome/ToxoDB-59_TgondiiME49_Genome.fasta', 
                       bed = peak.genes.union$peak.gene.merged.bed[,1:6], fo = out.fasta, s = T)


# ###############################
# ## add motif info to the peak gene table

peak.genes.union <- readRDS("../Input_sub/toxo_cdc/cutNrun_revision/rds/Union_all_new_peaks_0.05_qval_NO_frag_filt_3Controls.rds")
peak.genes <- peak.genes.union$peak.gene.merged.bed
peak.genes$assigned_to_CutRun_peaks <- "yes"
peak.genes <- peak.genes %>% dplyr::select(-c(ProductDescription, V5.x, V6.x))
 
# ## intensity scores
all.peak.genes.tab.wide <- readRDS("../Input_sub/toxo_cdc/cutNrun_revision/rds/cut_run_union_intensity_scores_3Controls_revision.rds")
all.peak.genes.tab.wide <- all.peak.genes.tab.wide[,c(1,3,4,5,6,7,8)]
all.peak.genes <- left_join(all.peak.genes.tab.wide, peak.genes , by = "gene_name" )



# ## motif info output from BAMM motif finder 
# 
#in.dir <- '../Input_sub/toxo_cdc/cutNrun/all_macs2_old_new_batch_NO_Filter/BAMM_analysis_v2/allPeaksUnion/peak_genes_union_all_BaMMmotif/' ## 970 (union)
#in.dir <- '../Input_sub/toxo_cdc/cutNrun_revision/all_macs2_old_new_batch_NO_Filter/BAMM_analysis_v2/allPeaksUnion/peak_genes_union_all_BaMMmotif/' ## 965 (union)
# occ.files <- list.files(in.dir)[grep('occurrence', list.files(in.dir))] # only motif 1 and motif 2
# in.file <- paste0(in.dir, occ.files)[-c(3,4)]


## motif info output from STREME motif finder 
all.motifs.tab <- readRDS("../Input_sub/Revision/all_motifs_streme_revision.rds")
peak.genes.motif <- all.motifs.tab[[which(grepl("peak_genes_union", names(all.motifs.tab)))]]

peak.genes.motif <- peak.genes.motif %>% 
  dplyr::select(motif_ALT_ID, motif_Sequence,motif_pvalue, seq_ID, data) 

peak.genes.motif <-peak.genes.motif %>% distinct(seq_ID, motif_ALT_ID, .keep_all = T) 
peak.genes.motif <- peak.genes.motif %>% filter(motif_ALT_ID == "STREME-2")
peak.genes.motif$motif_unique_ID <- gsub("STREME-2", "motif1", peak.genes.motif$motif_ALT_ID)
# unique(peak.genes.motif)

# peak.genes.motif$motif_unique_ID <- gsub("STREME", "Motif", gsub("-", "", peak.genes.motif$motif_ALT_ID))
# peak.genes.motif <- peak.genes.motif %>% mutate(MOTIF = case_when(motif_unique_ID == "Motif1" ~ "Motif2", 
#                                                                   motif_unique_ID == "Motif2" ~ "Motif1",
#                                                                   TRUE ~motif_unique_ID ))
peak.genes.motif <- peak.genes.motif %>% select(!motif_ALT_ID)

all.peak.genes$id <- paste0(all.peak.genes$V1.x, ":",
                            all.peak.genes$start_peak, "-",
                            all.peak.genes$end_peak, "(", all.peak.genes$V13, ")")
# 
# 
all.peak.genes.motif <- left_join(all.peak.genes, peak.genes.motif, by = c("id" = "seq_ID" ))
# 


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
cut.run.tab <- cut.run.tab %>% mutate(Category = ifelse(str_detect(ProductDescription, "ribosomal"), "ribosomal", "others"))

colnames(cut.run.tab) <- c("TGME49", "AP2XII-8_Ty_S4_vs_AP2XII-8_IgG1_peaks" , "AP2XII-8_Ty_S4_vs_RH_IgG1_S1_peaks" ,
                           "AP2XII-8_Ty_S4_vs_RH_Ty_S2_peaks",
                           "intensity.mean" ,   "intersection_CutRun_dataSets", "Ribosomal_description",
                           "chr", "start_peak" , "end_peak",
                           "V4", "V5", "V6","Genomic_location", "assigned_to_CutRun_peaks" , "id" , "pattern" ,"motif_pvalue","data","motif_id",
                           "KD_vs_WT_phase_based" ,"ProductDescription", "new.name" , "Category" )

# colnames(cut.run.tab) <- c("TGME49", "AP2XII-8_Ty_S4_vs_AP2XII-8_IgG1_peaks" ,
#                            "AP2XII-8_Ty_S4_vs_RH_IgG1_S1_peaks" ,"AP2XII-8_Ty_S4_vs_RH_Ty_S2_peaks",
#                            "intensity.mean" ,   "intersection_CutRun_dataSets", "Ribosomal_description",
#                            "Category", "chr", "start_peak" , "end_peak",
#                            "V4", "V5", "V6","Genomic_location", "assigned_to_CutRun_peaks" , 
#                            "KD_vs_WT_phase_based" ,"ProductDescription", "new.name"  )

write.xlsx(cut.run.tab, "../Input_sub/toxo_cdc/cutNrun_revision/rds/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v11_Streme_motif_revision.xlsx")
write.xlsx(cut.run.tab, "../OutPut/toxo_cdc/ME49_59/tables/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v11_Streme_motif_revision.xlsx")
saveRDS(cut.run.tab, "../Input_sub/toxo_cdc/cutNrun_revision/rds/cut_run_union_new_peaks_march_motif_modulated_genes_v11_Streme_motif_revision.rds")


#cut.run.tab <- cut.run.tab %>% mutate(Category = ifelse(str_detect(ProductDescription, "ribosomal"), "ribosomal", "others"))
# 
# write.xlsx(cut.run.tab, "../Input_sub/toxo_cdc/rds_ME49_59/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v9_revision.xlsx")
# write.xlsx(cut.run.tab, "../OutPut/toxo_cdc/ME49_59/tables/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v9_revision.xlsx")
# saveRDS(cut.run.tab, "../Input_sub/toxo_cdc/rds_ME49_59/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v9_revision.rds")
# 
# write.xlsx(cut.run.tab, "../Input_sub/toxo_cdc/rds_ME49_59/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v8.xlsx")
# write.xlsx(cut.run.tab, "../OutPut/toxo_cdc/ME49_59/tables/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v8.xlsx")
# saveRDS(cut.run.tab, "../Input_sub/toxo_cdc/rds_ME49_59/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v8.rds")




############## Update this part when KZ does the motif search ##################

## modify supplementary table 
#cut.run.tab <- readRDS("../Input_sub/toxo_cdc/rds_ME49_59/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v9_revision.rds")
cut.run.tab <- readRDS("../Input_sub/toxo_cdc/cutNrun_revision/rds/cut_run_union_new_peaks_march_motif_modulated_genes_v11_Streme_motif_revision.rds")

cut.run.tab.supplmnt <- cut.run.tab %>% 
  dplyr::select(id, chr, start_peak, end_peak, V4, V5, V6,TGME49, assigned_to_CutRun_peaks,
                intersection_CutRun_dataSets, MOTIF, KD_vs_WT_phase_based, ProductDescription, new.name, Category)
write.xlsx(cut.run.tab.supplmnt, "../OutPut/toxo_cdc/ME49_59/tables/Supplements_revision/cut_run_peak_gene_assignement_supplement_revision.xlsx")
#write.xlsx(cut.run.tab.supplmnt, "../OutPut/toxo_cdc/ME49_59/tables/Supplement/cut_run_peak_gene_assignement_supplement.xlsx")

## motif table summary Fig6 e

tab <- read.xlsx("../Input_sub/toxo_cdc/cutNrun_revision/rds/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v11_Streme_motif_revision.xlsx")

tab.down <- tab %>% 
  dplyr::select(intersection_CutRun_dataSets, TGME49,  motif_id, KD_vs_WT_phase_based, ProductDescription) %>%
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
Tg_ATAC <- readRDS('../Input_sub/toxo_cdc/rds_ME49_59/S.O_ATAC_peak.rds')

peak.regions <- data.frame(ATAC_peak_region = rownames(Tg_ATAC@assays$peaks@data))
peak.regions <- peak.regions %>% filter(!grepl("KE.*", ATAC_peak_region))
peak.regions.bed  <- data.frame(do.call(rbind, strsplit(peak.regions$ATAC_peak_region,"-")))
peak.regions.bed$X1 <- paste(peak.regions.bed$X1, peak.regions.bed$X2, sep = "_")
peaks.all.sort.atac <- peak.regions.bed  %>% dplyr::select(X1, X3, X4) %>%  arrange(X1, as.numeric(X3), as.numeric(X4))
names(peaks.all.sort.atac) <- c("V1", "V2", "V3")
peaks.all.sort.atac$V4 <- paste(paste(peaks.all.sort.atac$V1, peaks.all.sort.atac$V2, sep = ":"),peaks.all.sort.atac$V3 ,sep = "-" )
peak.genes.bed.merged.bed <- peaks.all.sort.atac

saveRDS(peak.genes.bed.merged.bed,"../Input_sub/toxo_cdc/rds_ME49_59/atac_peaks.rds")


## overlap cut&run and atac genes fig 5 b 
peak.genes.atac <- read.table("../Input_sub/toxo_scATAC_MJ_ME49_59/peak_gene_assigned_final.bed")
peak.genes.cutRun <- readRDS("../Input_sub/toxo_cdc/cutNrun_revision/rds/Union_all_new_peaks_0.05_qval_NO_frag_filt_3Controls.rds")
peak.genes.cutRun <- peak.genes.cutRun$peak.gene.merged.bed[1:7] 

ovlp <- inner_join(peak.genes.cutRun, peak.genes.atac, by = c("gene_name" = "V7"))
# venn.list <- list(atac.genes = unique(peak.genes.atac$V7),
#                   cutRun.genes = unique(peak.genes.cutRun$gene_name))
# ggVennDiagram(venn.list)
library(VennDiagram)
pdf(file = "../Output/toxo_cdc/ME49_59/figures_paper/RevisionUpdate/cut_run_genes_3Controls_overlap_atac_genes_venn_revision.pdf",
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

## regenerate fig 6 

## 6 b

tab <- read.xlsx("../Input_sub/toxo_cdc/cutNrun_revision/rds/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v11_Streme_motif_revision.xlsx")

dir.targets <- tab %>% filter(intersection_CutRun_dataSets == "yes" & KD_vs_WT_phase_based != "NA")
dim(dir.targets)

## DEGs KD vs WT phase based 
KD.vs.WT.phase <- tab %>% dplyr::select(TGME49, KD_vs_WT_phase_based) %>% 
  filter(KD_vs_WT_phase_based != "NA") %>% distinct()
dim(KD.vs.WT.phase)
DEGs <- KD.vs.WT.phase

# cutRun genes in intersection of 4 data
intrsct.peaks  <- tab %>% dplyr::select(TGME49,intersection_CutRun_dataSets) %>% 
  distinct() %>% filter(intersection_CutRun_dataSets == "yes")
dim(intrsct.peaks)

## overlap - venn
# DEGs <- KD.vs.WT.phase
# venn.list <- list(KD.vs.WT = unique(DEGs$TGME49),
#                   High.Conf.Peaks = intrsct.peaks$TGME49)
# 
# p <- ggVennDiagram(venn.list, set_size = 6, label_size = 8)
# p



pdf(file = "../Output/toxo_cdc/ME49_59/figures_paper/RevisionUpdate/KD_vs_WT_High_Conf_Peaks_overlaps.pdf",
    width = 12, height = 12)
venn.plot <- draw.pairwise.venn(
  area1 =length(unique(DEGs$TGME49)),
  area2 = length(unique(intrsct.peaks$TGME49)) ,
  cross.area = length(unique(dir.targets$TGME49)),
  #category = c("ATAC", "C&R"),
  fill = c("#A6CDFB","#BAC788"),
  lty = rep("solid", 2),
  lwd = 6,
  col = c("#2E4CA8", "darkgreen"),
  cex = 5.5,
  cat.cex = 3,
  ext.length = 0.9,
  ext.line.lwd = 2.5,
  
)
grid.draw(venn.plot)
dev.off()




### 
## motif search 
###############################################################
## write fasta file (cut and run regions) for high conf peaks 
## write the genes in excel file for GO term analysis on toxodb
###############################################################
#tab <- readRDS("../Input/toxo_cdc/rds_ME49_59/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v4.rds")
tab <- read.xlsx("../Input_sub/toxo_cdc/cutNrun_revision/rds/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v11_Streme_motif_revision.xlsx")

# phase based 
HC.peaks <- tab %>% 
  filter(intersection_CutRun_dataSets == "yes" & KD_vs_WT_phase_based %in% c('down_reg', 'up_reg', 'modulated') ) %>% 
  dplyr::select(chr, start_peak, end_peak, V4, V5, V6,TGME49,intersection_CutRun_dataSets, 
                KD_vs_WT_phase_based,ProductDescription  ) %>% 
  distinct()
names(HC.peaks)[9] <- "dir"


HC.peaks.list <- split(HC.peaks, f = HC.peaks$dir)

DEG.type <- "KD_vs_WT_phase_based"


#out.dir <- "../Input/toxo_cdc/cutNrun/all_macs2_old_new_batch_NO_Filter/BAMM_analysis/HighConfPeaks/bedFasta/"
#out.dir <- "../Input/toxo_cdc/cutNrun/all_macs2_old_new_batch_NO_Filter/BAMM_analysis/HighConfPeaks_globalDEG/bedFasta/"
out.dir <- "../Input_sub//toxo_cdc/cutNrun_revision/all_macs2_old_new_batch_NO_Filter/BAMM_analysis_v2/HighConfPeaks/bedFasta/"
cluster.bed.list <- lapply(1:length(HC.peaks.list), function(i){
  
  bed.name <- paste(names(HC.peaks.list)[i], paste(DEG.type, ".bed", sep = ""), sep = "_")
  fasta.name <- paste(names(HC.peaks.list)[i], paste(DEG.type, ".fasta", sep = ""), sep = "_")
  tmp <- HC.peaks.list[[i]] 
  cluster.bed <- tmp %>% select(chr, start_peak, end_peak ,V4, V5 ,  V6, TGME49,  dir)
  
  write.table(cluster.bed, paste(out.dir, bed.name, sep = ""), sep = "\t", quote = F, row.names = F, col.names = F)
  write.xlsx(cluster.bed, paste(out.dir, gsub('.bed', '.xlsx', bed.name), sep = "") )
  
  cluster.fasta <- bedtoolsr::bt.getfasta(fi = '../Input_sub/Toxo_genomics/genome/ToxoDB-59_TgondiiME49_Genome.fasta',
                                          bed = tmp,
                                          fo =paste(out.dir, fasta.name, sep = ""), s = T)
  return(cluster.bed)
})


##########


cut.run.tab <- read.xlsx("../Input_sub/toxo_cdc/cutNrun_revision/rds/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v11_Streme_motif_revision.xlsx")

