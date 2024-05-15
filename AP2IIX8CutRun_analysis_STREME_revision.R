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


prod.desc  <- read.xlsx('../Input_sub/rds_ME49_59/toxo_genomics/genes/ProductDescription_ME49.xlsx')
ribosomals <- prod.desc[grep('ribosomal', prod.desc$ProductDescription),]
ribosomal.proteins <- prod.desc[grep('ribosomal protein', prod.desc$ProductDescription),]
ribosomal.proteins <- ribosomal.proteins[-grep('KE', ribosomal.proteins$GenomicLocation),]
ribosomal.proteins <- ribosomal.proteins %>% dplyr::select(GeneID, ProductDescription) ## 137 Ribosomal Proteins

in.dir <- "../Input_sub/rds_ME49_59/cutNrun/all_macs2_old_new_batch_NO_Filter/"
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



# peak gene assignment 
gtf.file <- "../Input_sub/rds_ME49_59/Toxo_genomics/genome/ToxoDB-59_TgondiiME49.gtf"
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


## Ven diagram - overalp of peaks across 3 cut and run data 

all.peak.genes.tab <- readRDS("../Input_sub/rds_ME49_59/cutNrun_revision/rds/CutRunPeaks_3Controls_genes_assigned_revision.rds")
all.peak.genes.tab.list <- split(all.peak.genes.tab, f = all.peak.genes.tab$data)
names(all.peak.genes.tab.list)
peaks.venn.list <- list(AP2XII8_vs_AP2XII8_IgG1 = unique(all.peak.genes.tab.list$`AP2XII-8_Ty_S4_vs_AP2XII-8_IgG1_peaks`$gene_name),
                        AP2XII8_vs_RH_IgG1 = unique(all.peak.genes.tab.list$`AP2XII-8_Ty_S4_vs_RH_IgG1_S1_peaks`$gene_name),
                        AP2XII8_vs_RH_Ty = unique(all.peak.genes.tab.list$`AP2XII-8_Ty_S4_vs_RH_Ty_S2_peaks`$gene_name))

p <- ggVennDiagram(peaks.venn.list, label = "both", label_alpha = 0.2,
                   label_color = "white", label_size = 14, set_size = 4)
p



all.pesks.ribo <- all.peak.genes.tab %>% filter(str_detect(ProductDescription, "ribosomal protein"))
all.pesks.ribo.list <- split(all.pesks.ribo, f= all.pesks.ribo$data)
names(all.pesks.ribo.list)
ribo.peaks.venn.list <- list(AP2XII8_vs_AP2XII8_IgG1 = unique(all.pesks.ribo.list$`AP2XII-8_Ty_S4_vs_AP2XII-8_IgG1_peaks`$gene_name),
                             AP2XII8_vs_RH_IgG1 = unique(all.pesks.ribo.list$`AP2XII-8_Ty_S4_vs_RH_IgG1_S1_peaks`$gene_name),
                             AP2XII8_vs_RH_Ty = unique(all.pesks.ribo.list$`AP2XII-8_Ty_S4_vs_RH_Ty_S2_peaks`$gene_name))
p <- ggVennDiagram(ribo.peaks.venn.list, label = "both", label_alpha = 0.2,
                   label_color = "white", label_size = 14, set_size = 4)
p


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



## concatenating peaks (union of 3 new data sets peaks)
in.dir <- "../Input_sub/rds_ME49_59/cutNrun/all_macs2_old_new_batch_NO_Filter/"
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


gtf.file <- "../Input_sub/rds_ME49_59/Toxo_genomics/genome/ToxoDB-59_TgondiiME49.gtf"
gtf <- read.table(gtf.file, header = F, sep = '\t', quote = NULL)
qval <- 0.05

# AP2_Ty vs 3 controls

peak.genes.union <- get_peak_genes_assign(gtf,all.peaks, qval)


