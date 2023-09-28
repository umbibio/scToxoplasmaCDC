library(tidyverse)
library(bedtoolsr)
library(openxlsx)

## looked at the cut and run ribosomal genes at the intersection of 4 data sets 

## IDs

prod.desc  <- read.xlsx('../Input_sub/Toxo_genomics/genes/ProductDescription_ME49.xlsx')
ribosomals <- prod.desc[grep('ribosomal', prod.desc$ProductDescription),]

## NOTES: 
## 1. Some of the Ribosomal Proteins and particularly Ribosomal RNAs are located on contigs other than 
## the main Chromosomes. They will not show up in our analysis.
## 2. There is a discrepancy in product description between GT1 and ME49. Some Ribosomals have ID in only one
## 3. Many of the non-called peaks indeed have a peak in cut&run AP2XII-8_ty over AP2XII-8_IgG. 
## 4. Maybe we need to re-map the reads to the genome using the latest version of ME49 and it's annotation
## and not exclude other contigs, and revisit fine tuning the cut&run negative control to more accurately
## identify the peaks. 
## For now, we focus on the most-stringent and Ribosomal Proteins on the main chromosomes.
ribosomal.proteins <- prod.desc[grep('ribosomal protein', prod.desc$ProductDescription),]
ribosomal.proteins <- ribosomal.proteins[-grep('KE', ribosomal.proteins$GenomicLocation),]
ribosomal.proteins <- ribosomal.proteins %>% dplyr::select(GeneID, ProductDescription) ## 137 Ribosomal Proteins

## AP2XII-8 targets
## Direct targets of AP2XII8
#targs <- read.xlsx("../OutPut/toxo_cdc/ME49_59/tables/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v5.xlsx")
#targs <- read.xlsx("../OutPut/toxo_cdc/ME49_59/tables/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v6.xlsx")
targs <- read.xlsx("../OutPut/toxo_cdc/ME49_59/tables/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v8.xlsx")
## we dont use this as we need the union of cut and run ribosomal genes
##targs <- targs %>% filter(intersection_CutRun_dataSets == "yes")

targs <- full_join(targs, ribosomal.proteins, by = c('TGME49' = 'GeneID'))

cut_run_ribo_targs <- targs %>% dplyr::filter(!is.na(ProductDescription.y))

## Percent of Ribosomal Protein with cut&run peaks
cut_run_ribo_targs <- cut_run_ribo_targs %>% 
  dplyr::select("chr","start_peak", "end_peak", "V4", "V5", "V6", "TGME49", 
                'intersection_CutRun_dataSets','KD_vs_WT_phase_based', 'new.name', 'ProductDescription.y') %>% 
  distinct()


cut_run_ribo_targs.filt <- cut_run_ribo_targs[!is.na(cut_run_ribo_targs$chr),]
cut_run_ribo_targs.bed <- cut_run_ribo_targs.filt[,1:6]
cut_run_ribo_targs.bed.ext <- cut_run_ribo_targs.bed
cut_run_ribo_targs.bed.ext$start_peak <- cut_run_ribo_targs.bed$start_peak - 250 
cut_run_ribo_targs.bed.ext$end_peak <- cut_run_ribo_targs.bed$end_peak + 250

write.table(x = cut_run_ribo_targs.bed.ext, file = '../Input_KZ/toxo_cdc/cut_run_union_ribo_targs.bed', sep = '\t', quote = F, row.names = F, col.names = F)
bedtoolsr::bt.getfasta('../Input/toxo_genomics/genome/ToxoDB-59_TgondiiME49_Genome.fasta', '../Input_KZ/toxo_cdc/cut_run_union_ribo_targs.bed', s = T, fo = '../Input_KZ/toxo_cdc/cut_run_union_ribo_targs.fa')

## Running the above output on https://bammmotif.soedinglab.org/
## Analyzing the output of motif search

## From the developer:
#“The start and end possition here indicate pseudo-coordinates here. 
# The tool first append the (-) strand in a reverse complementary manner to the (+) strand, 
# seprated by a letter N, i.e. [+strand sequence]N[-strand sequence], and then scan the sequence for motifs. 
# Thus in your case, the motif is found on the - strand, which is 459 bp (2x885-1311) to 467 bp (2x885-1303) 
# away from the sequence start position (3756384).”

## Formula for negative strand: 
## motif.start = 2 * l - motif.stop
## motif.stop = 2 * l - motif.start
#in.dir <- '../Input_KZ/toxo_cdc/cut_run_ribo_targs_BaMMmotif_250_KZ/'
#in.dir <- '../Input_KZ/toxo_cdc/cut_run_stringent_ribo_targs_BaMMmotif/' ## 54 ribo cut and run (intersection)

in.dir <- '../Input_KZ/toxo_cdc/cut_run_union_ribo_targs_BaMMmotif/' ## 76 ribo cut run (union)
occ.files <- list.files(in.dir)[grep('occurrence', list.files(in.dir))]
in.file <- paste0(in.dir, occ.files)
all.motifs.tab <- lapply(1:length(in.file), function(i){
  mofit <- gsub('cut_run_union_ribo_targs_', '', gsub('\\..*', '', occ.files[i]))
  
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

cut_run_ribo_targs.filt$id <- paste0(cut_run_ribo_targs.filt$chr, ":",
                                as.numeric(cut_run_ribo_targs.filt$start_peak) - 250, "-",
                                as.numeric(cut_run_ribo_targs.filt$end_peak) + 250, "(",
                                cut_run_ribo_targs.filt$V6, ")")

cut_run_ribo_targs.filt <- cut_run_ribo_targs.filt %>% distinct()
all.motifs.tab.genes <- full_join(all.motifs.tab, cut_run_ribo_targs.filt, by = 'id')
all.motifs.tab.genes <- all.motifs.tab.genes %>% arrange(TGME49, motif.strt.rel, motif.stp.rel)
all.motifs.tab.genes <- all.motifs.tab.genes[!is.na(all.motifs.tab.genes$motif),]


length(unique(all.motifs.tab.genes$TGME49[!is.na(all.motifs.tab.genes$motif)])) 


## ribosomal genes 
AP2XII8_ribo_TATA <- all.motifs.tab.genes 
AP2XII8_ribo_TATA <- AP2XII8_ribo_TATA %>% 
  dplyr::select(TGME49, motif, pattern, intersection_CutRun_dataSets,KD_vs_WT_phase_based, new.name, 
                ProductDescription.y) %>%
  distinct()
AP2XII8_ribo_TATA %>% group_by(motif) %>% 
  summarise(num.genes = length(unique(TGME49)),
            percent.genes = length(unique(TGME49)) /length(unique(all.motifs.tab.genes$TGME49)))

write.xlsx(AP2XII8_ribo_TATA, "../OutPut/toxo_cdc/ME49_59/tables/AP2XII8_ribosomals_TATA_motif.xlsx")

## NEW
## TSS dist 
all.motifs.tab.genes.motif <- all.motifs.tab.genes %>% 
  select(id, TGME49, motif, chr.x, pattern, motif.strt.abs, motif.stp.abs, ProductDescription.y)

# ## distance of motifs relative to TSS
gtf.file <- "../Genomes/ToxoDB-59_TgondiiME49.gtf"
gtf <- read.table(gtf.file, header = F, sep = '\t', quote = NULL)
gtf.filt <- gtf %>% dplyr::filter(!grepl('KE.*', V1))
gtf.filt$gene_name <- gsub("\\;.*", "", gsub("transcript_id ", "", gsub("-t.*", "", gtf.filt$V9)))
gtf.filt$gene_name <- gsub("\"", "", gtf.filt$gene_name)
gtf.filt <- gtf.filt %>% filter(V3 == "transcript")

all.motifs.tab.genes.motif <- left_join(all.motifs.tab.genes.motif, gtf.filt, by = c("TGME49" = "gene_name") )

all.motifs.tab.genes.motif <-  all.motifs.tab.genes.motif %>% mutate(tss = ifelse((V7 == "+"), V4, V5))
all.motifs.tab.genes.motif <- all.motifs.tab.genes.motif %>% mutate(tss.dist = tss - motif.strt.abs)


p <- ggplot(all.motifs.tab.genes.motif, aes(x = tss.dist)) + 
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "white") +
  geom_density(alpha=.2, fill="#FF6666",lwd = 0.7,linetype = 1,colour = 2) +
  theme_bw() + 
  facet_grid(motif ~ .) +
  theme_bw()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,size = 16, face = 'bold', color = 'black'))+
  theme(axis.text.y = element_text(size = 16, face = 'bold', color = 'black'))+
  theme(axis.title = element_text(size = 20, face = "bold", color = "black")) +
  theme(strip.text = element_text(size = 14, face = "bold", colour = "black"))
p

ggsave("../OutPut/toxo_cdc/ME49_59/figures_paper/ribo_motifs/ribo_cut_run_union_motif_TSS_dist_BAMM.pdf", 
       plot = p,height = 6, width = 4, dpi = 300)

# x <- all.motifs.tab.genes %>% dplyr::filter(motif == 'motif_2') %>% group_by(TGME49) %>% 
#   mutate(both.strand = ifelse('+' %in% motif.strd & '-' %in% motif.strd, T, F)) %>% dplyr::filter(both.strand) %>%
#   arrange(TGME49, motif.strt.rel, motif.strt.rel) %>% dplyr::select(TGME49, motif.strt.rel, motif.strt.rel, motif.strd)

## Calculating a few stats

## Percent of Ribosomal Protein with cut&run peaks
## Total of 137 Ribosomal Proteins on the main chromosomes
## 54 of which (~40%) have cut&run peaks with high stringency
## Of the 76 with CUT&RUN peaks, 71 have motif

length(unique(cut_run_ribo_targs$TGME49[!is.na(cut_run_ribo_targs$chr)])) / length(unique(cut_run_ribo_targs$TGME49))
length(unique(all.motifs.tab.genes$TGME49[!is.na(all.motifs.tab.genes$motif)])) 

gg <- all.motifs.tab.genes %>% group_by(TGME49) %>% 
  summarise(num.motifs = n())

## Proportion of genes with various combination of motifs
mm <- all.motifs.tab.genes %>% dplyr::filter(!is.na(motif)) %>% group_by(TGME49) %>% 
  summarise(num.motifs = n(), motifs = list(motif))

nrow(mm)

## Proportion of genes with >= 2 motifs
sum(gg$num.motifs > 1) / nrow(gg)
p1 <- ggplot(data = gg, aes(x = num.motifs)) + 
  geom_histogram(aes(y=..density..), binwidth = 1, colour="black", fill="gray95")+
  geom_density(alpha=.2, fill="#FF6666") +
  theme_bw(base_size = 14) +
  ylab('density') + xlab('motif totals') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  theme(strip.text = element_text(size = 14, face="bold", angle = 0)) + 
  theme(
    plot.title = element_text(size=14, face = "bold.italic", color = 'black'),
    axis.title.x = element_text(size=14, face="bold", hjust = 1),
    axis.title.y = element_text(size=14, face="bold")
  ) + 
  theme(#legend.position = c(0.15, 0.85),
    legend.position = 'none',
    legend.title = element_text(colour="black", size=12, 
                                face="bold"),
    legend.text = element_text(colour="black", size=12, 
                               face="bold"))

p1

ggsave(filename="../Output/toxo_cdc/ME49_59/figures_paper/ribo_motifs/motif_totol_ribosomals_union.pdf", plot=p1,
       width = 4, height = 4, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)

## Proportion with the given motif
ss <- all.motifs.tab.genes %>% group_by(motif) %>% 
  summarise(num.genes = length(unique(TGME49)), 
            percent.genes = length(unique(TGME49)) /length(unique(all.motifs.tab.genes$TGME49)))


ss$motif <- gsub('motif_', 'm', ss$motif)
p2 <- ggplot(data = ss, aes(x = motif, y = percent.genes, colour = motif, fill = motif)) + 
  geom_bar(stat="identity")+
  theme_bw(base_size = 14) +
  geom_text(aes(label=round(percent.genes, 3)), vjust=1.6, color="black", size=6)+
  ylab('proportions') + xlab('motif') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  theme(strip.text = element_text(size = 14, face="bold", angle = 0)) + 
  theme(
    plot.title = element_text(size=14, face = "bold.italic", color = 'black'),
    axis.title.x = element_text(size=14, face="bold", hjust = 1),
    axis.title.y = element_text(size=14, face="bold")
  ) + 
  theme(#legend.position = c(0.15, 0.85),
    legend.position = 'none',
    legend.title = element_text(colour="black", size=12, 
                                face="bold"),
    legend.text = element_text(colour="black", size=12, 
                               face="bold"))





p2

ggsave(filename="../Output/toxo_cdc/ME49_59/figures_paper/ribo_motifs/motif_proportions_union.pdf", plot=p2,
       width = 4, height = 4, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)



## Calculating motif occurrence
x <- all.motifs.tab.genes %>% group_by(TGME49) %>% summarise(m1 = sum(motif == 'motif_1'),
                                                             m2 = sum(motif == 'motif_2'),
                                                             m3 = sum(motif == 'motif_3')
                                                             )

x


## off diagnoal are co-occurrence
mm <- as.matrix(x[,2:ncol(x)])
self_oc <- mm
mm[mm != 0] = 1
co.mat <- t(mm) %*% mm
## Diagnal is total genes including the motif. Update to self co-oc
self_oc[self_oc < 2] <- 0
self_oc[self_oc != 0] <- 1
diag(co.mat) <- colSums(self_oc)

co.mat.long <- data.frame(co.mat) %>% mutate(Motif1 = rownames(co.mat)) %>% 
  pivot_longer(-Motif1, names_to = 'Motif2', values_to = 'co_oc')
co.mat.long.filt <- co.mat.long %>% dplyr::filter(Motif1 <= Motif2)
co.mat.long.filt$scaled <- (co.mat.long.filt$co_oc - min(co.mat.long.filt$co_oc)) / 
  (max(co.mat.long.filt$co_oc) - min(co.mat.long.filt$co_oc))

ggheatmap <- ggplot(co.mat.long.filt, aes(Motif1, Motif2, fill = scaled))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       space = "Lab", 
                       midpoint = quantile(co.mat.long.filt$scaled,  probs =  0.6, na.rm = T), limit = c(0,1),
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 18, hjust = 0, face = 'bold', color = 'black'))+
  theme(axis.text.y = element_text(vjust = 1, 
                                   size = 18, hjust = 1,face = 'bold', color = 'black')) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    #panel.grid.major = element_blank(),
    panel.border = element_blank(),
    #panel.background = element_blank(),
    #axis.ticks = element_blank(),
    legend.position = 'none') +
  coord_fixed() + geom_tile(colour="white",size=4) + 
  geom_text(data = co.mat.long.filt, aes(Motif1, Motif2, label = co_oc), color = "black", 
            size = 8, fontface='bold')
# Print the heatmap
print(ggheatmap)

ggsave(filename="../Output/toxo_cdc/ME49_59/figures_paper/ribo_motifs/motif_co_occ_union.pdf", plot=ggheatmap,
       width = 4, height = 4, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)




## With strand info
y <- all.motifs.tab.genes  %>% group_by(TGME49) %>% summarise(m1pos = sum(motif == 'motif_1' & motif.strd == '+'),
                                          m1neg = sum(motif == 'motif_1' & motif.strd == '-'),
                                          m2pos = sum(motif == 'motif_2' & motif.strd == '+'),
                                          m2neg = sum(motif == 'motif_2' & motif.strd == '-'),
                                          m3pos = sum(motif == 'motif_3' & motif.strd == '+'),
                                          m3neg = sum(motif == 'motif_3' & motif.strd == '-')
                                          )

## Calculating motif occurrence on + and -
## Diagnal is total genes including the motif on the indicated strand
## off diagnoal are co-occurrence
## Seems like motif 1 co-occurres on both strands? homo-dimerization?

nn <- as.matrix(y[,2:ncol(y)])
self_oc_strd <- nn
nn[nn != 0] = 1
## off diagnoal are co-occurrence

co.mat.strd <- t(nn) %*% nn

## Diagnal is total genes including the motif. Update to self co-oc
self_oc_strd[self_oc_strd < 2] <- 0
self_oc_strd[self_oc_strd != 0] <- 1
diag(co.mat.strd) <- colSums(self_oc_strd)

co.mat.long.strd <- data.frame(co.mat.strd) %>% mutate(Motif1 = rownames(co.mat.strd)) %>% 
  pivot_longer(-Motif1, names_to = 'Motif2', values_to = 'co_oc')
co.mat.long.filt.strd <- co.mat.long.strd %>% dplyr::filter(Motif1 <= Motif2)
co.mat.long.filt.strd$scaled <- (co.mat.long.filt.strd$co_oc - min(co.mat.long.filt.strd$co_oc)) / 
  (max(co.mat.long.filt.strd$co_oc) - min(co.mat.long.filt.strd$co_oc))

ggheatmap <- ggplot(co.mat.long.filt.strd, aes(Motif1, Motif2, fill = scaled))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "gray95", 
                       space = "Lab", 
                       midpoint = quantile(co.mat.long.filt.strd$scaled,  probs =  0.6, na.rm = T), limit = c(0,1),
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 18, hjust = 0, face = 'bold', color = 'black'))+
  theme(axis.text.y = element_text(vjust = 1, 
                                   size = 18, hjust = 1,face = 'bold', color = 'black')) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    #panel.grid.major = element_blank(),
    panel.border = element_blank(),
    #panel.background = element_blank(),
    #axis.ticks = element_blank(),
    legend.position = 'none') +
  coord_fixed() + geom_tile(colour="white",size=4) + 
  geom_text(data = co.mat.long.filt.strd, aes(Motif1, Motif2, label = co_oc), color = "black", 
            size = 8, fontface='bold')
# Print the heatmap
print(ggheatmap)

ggsave(filename="../OutPut/toxo_cdc/ME49_59/figures_paper/ribo_motifs/motif_co_occ_strand_union.pdf", plot=ggheatmap,
       width = 4, height = 4, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)



## Successive motif distances
w <- all.motifs.tab.genes %>% ungroup() %>% 
  dplyr::select(TGME49, motif, pattern, motif.strt.rel, motif.stp.rel, motif.strd) %>% 
  arrange(TGME49, motif.strt.rel, motif.stp.rel) %>% 
  group_by(TGME49) %>% 
  summarise(m_i = list(motif[1:(length(motif)-1)]), m_ip1 = list(motif[2:length(motif)]),
            p_i = list(pattern[1:(length(pattern)-1)]), p_ip1 = list(pattern[2:length(pattern)]),
            m_i_strt = list(motif.strt.rel[1:(length(motif.strt.rel)-1)]), m_ip1_strt = list(motif.strt.rel[2:length(motif.strt.rel)]),
            m_i_stp = list(motif.stp.rel[1:(length(motif.stp.rel)-1)]), m_ip1_stp = list(motif.stp.rel[2:length(motif.stp.rel)]),
            m_i_strd = list(motif.strd[1:(length(motif.strd)-1)]), m_ip1_strd = list(motif.strd[2:length(motif.strd)])) %>% 
  unnest(cols = c(m_i, m_ip1, p_i, p_ip1, m_i_strt, m_ip1_strt, m_i_stp, m_ip1_stp, m_i_strd, m_ip1_strd)) %>%
  mutate(dist = m_ip1_strt - m_i_stp, 
         strand = ifelse(m_i_strd < m_ip1_strd, paste(m_i_strd, m_ip1_strd, sep = '/'), paste(m_ip1_strd, m_i_strd, sep = '/')),
         co_motif = ifelse(m_i < m_ip1, paste(m_i, m_ip1, sep = '/'), paste(m_ip1, m_i, sep = '/'))) %>%
  na.omit()
  
w$d.scaled <- (w$dist - min(w$dist)) / 
  (max(w$dist) - min(w$dist))

w$co_motif <- gsub('motif_', 'm', w$co_motif)
p <- ggplot(data = w, aes(x = co_motif, y = dist, color = co_motif, fill = co_motif)) + 
  geom_boxplot() + 
  geom_jitter(color="black", size=0.5, alpha=0.9) + 
  facet_grid(strand~.) + 
  ylab('dist(nucleotide)') + xlab("") + 
  theme_bw()+ # minimal theme
  theme(axis.text.x = element_text(angle = 90, vjust = 0, 
                                   size = 12, hjust = 1, face = 'bold', color = 'black'))+
  theme(axis.text.y = element_text(vjust = 1, 
                                   size = 12, hjust = 1,face = 'bold', color = 'black')) +
  theme(
    axis.title.y = element_text(size = 14, face = 'bold', color = 'black'),
    strip.text.y = element_text(size = 14, face = 'bold', color = 'black'),
    #axis.title.x = element_blank(),
    #axis.title.y = element_blank(),
    #panel.grid.major = element_blank(),
    #panel.border = element_blank(),
    #panel.background = element_blank(),
    #axis.ticks = element_blank(),
    legend.position = 'none') 
 
p
ggsave(filename="../OutPut/toxo_cdc/ME49_59/figures_paper/ribo_motifs/motif_dist_distrib_stringent.pdf", plot=p,
       width = 10, height = 5, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)


