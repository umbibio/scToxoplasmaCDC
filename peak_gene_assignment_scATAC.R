
source('./util_funcs.R')
source('./loadlb.R')


# ATAC peak Regions in bed format 

Tg_ATAC <- readRDS('../Input_sub/toxo_cdc/rds_ME49_59/S.O_ATAC_peak.rds')

peak.regions <- data.frame(ATAC_peak_region = rownames(Tg_ATAC@assays$peaks@data))
peak.regions <- peak.regions %>% filter(!grepl("KE.*", ATAC_peak_region))
peak.regions.bed  <- data.frame(do.call(rbind, strsplit(peak.regions$ATAC_peak_region,"-")))
peak.regions.bed$X1 <- paste(peak.regions.bed$X1, peak.regions.bed$X2, sep = "_")
peaks.all.sort <- peak.regions.bed  %>% dplyr::select(X1, X3, X4) %>%  arrange(X1, as.numeric(X3), as.numeric(X4))
names(peaks.all.sort) <- c("V1", "V2", "V3")
peaks.all.sort$V4 <- paste(paste(peaks.all.sort$V1, peaks.all.sort$V2, sep = ":"),peaks.all.sort$V3 ,sep = "-" )



# annotation file 

gtf.file <- "../Input_sub/Genomes/ToxoDB-59_TgondiiME49.gtf"
gtf <- read.table(gtf.file, header = F, sep = '\t', quote = NULL)
gtf.filt <- gtf %>% dplyr::filter(!grepl('KE.*', V1))


# remove the first Exon from transcripts.

gtf.exon <- gtf.filt %>% dplyr::filter(V3 == 'exon')
gtf.exon.sort <- gtf.exon %>% arrange(V1, V4, V5)
parse.str <- strsplit(gtf.exon$V9, split = ' ')
inds <- unlist(lapply(parse.str , function(x) which(grepl("gene_id", x)) + 1))
gtf.exon$gene_name <- gsub(";", "", unlist(lapply(1:length(inds), function(i) parse.str[[i]][inds[[i]]])))
gtf.exon$gene_name <- gsub("\"", "", gtf.exon$gene_name)
gtf.exon <- gtf.exon %>% group_by(V9) %>% mutate(exon.ord = ifelse(V7 == '+', 1:n(), seq(n(), 1, by = -1)),
                                                 multiple.exon = ifelse(n() > 1, T, F))

                                                  
# remove the exon1, but keep the Intron 1 , build exon2ton
                                                  
gtf.exon.2Ton <- gtf.exon %>% mutate(V10 = ifelse(multiple.exon & V7 == '-', min(V4), 
                                                  ifelse(multiple.exon & V7 == '+', min(V5), V4)),
                                     V11 = ifelse(multiple.exon & V7 == '-', max(V4), 
                                                  ifelse(multiple.exon & V7 == '+', max(V5), V5))) %>%
  mutate(V4 = V10, V5 = V11) %>% 
  dplyr::select(-c(exon.ord,multiple.exon, V10, V11) ) %>% distinct()


# ATAC Peak gene assignment

## Overlap with peaks and filter peaks that are entirely within the genes.
## Overlapping peaks with exon2Ton and check to see if it is entirely within the gene. 
## Then from the sorted peaks we throw out all peaks entirely within the gene & not overlapping with exon1, 
## These peaks should not be assigned to any peaks.


peak.genes.ovrlp <- bedtoolsr::bt.intersect(a = peaks.all.sort, b = gtf.exon.2Ton, wo = T)
peak.genes.filt <- peak.genes.ovrlp %>% dplyr::filter(V8  <= V2 & V9 >= V3)
peak.filt <- peaks.all.sort[!(peaks.all.sort$V4 %in%  peak.genes.filt$V4), ]
peak.filt.sort <- peak.filt %>% arrange(V1, as.numeric(V2), as.numeric(V3))


# filter gtf for transcripts only to get the coordinates of the start and end of gene
                                                  
gtf.filt.trn <- gtf.filt %>% filter(V3 == "transcript")
gtf.filt.trn$gene_name <- gsub("\\;.*", "", gsub("transcript_id ", "", gsub("-t.*", "", gtf.filt.trn$V9)))
gtf.filt.trn$gene_name <- gsub("\"", "", gtf.filt.trn$gene_name)


# filter for first exon coordinates (exon1 coordinates)
                                                  
tmp.neg <- gtf.exon %>% filter(V7 == "-") %>% group_by(V9) %>%  dplyr::slice(which.max(V5))
tmp.pos <- gtf.exon %>% filter(V7 == "+") %>% group_by(V9) %>%  dplyr::slice(which.min(V5))
gtf.exon1 <- bind_rows(tmp.pos, tmp.neg)
gtf.exon1.sort <- gtf.exon1 %>% arrange(V1, V4, V5)


# assign the peaks to the nearest upstream gene (look at the 5 closest in case of bi-directional)

peaks.genes.dist <- bedtoolsr::bt.closest(a = peak.filt.sort, b = gtf.exon1.sort, D = "b", k = 5)
parse.str2 <- strsplit(peaks.genes.dist$V13, split = ';')
peaks.genes.dist$gene_name  <- unlist(lapply(parse.str2, '[[' , 3))
peaks.genes.dist.trns <- left_join(peaks.genes.dist, gtf.filt.trn, by = "gene_name")

                                                  
## V16 is the distance of the peak to the exon 1 
## we need to overcome the issue with the  ones with  dist = 0
## on pos strand V3.x (end of peak) should not exceed V5.y (end of transcript/exon_n)
## on neg strand V2.x (start of peak) is not less than V4.y (beggining of the transcript/exon_1)

                                                  
peaks.genes.dist.trns <- peaks.genes.dist.trns %>% filter(!(V16 == 0 & V11 == "+" & V3.x > V5.y))

peaks.genes.dist.trns <- peaks.genes.dist.trns %>% filter(!(V16 == 0 & V11 == "-" & V2.x < V4.y))

## V16 <= 0 means the peak is at upstream 
## Find the closest gene among top 5 that is upstream (min V16)
                                                  
peaks.genes.dist.trns <- peaks.genes.dist.trns %>% dplyr::filter(V16 <= 0)
peaks.genes.dist.trns <- peaks.genes.dist.trns %>% group_by(V4.x) %>% 
  mutate(V17 = V16[which.min(abs(V16))])


# filter the rest
peaks.genes.dist.trns <- peaks.genes.dist.trns %>% dplyr::filter(V16 == V17)

# filter the ones that are too far (3000 bp)
peaks.genes.dist.trns.filt <- peaks.genes.dist.trns %>% dplyr::filter(abs(V16) < 3000)


# prep bed format 
# merge multiple peaks assigned to a single gene

peak.genes <- read.xlsx("../Input_sub/toxo_scATAC_MJ_ME49_59/scATAC_peaks_genes_assigned.xlsx")
peak.genes <- peak.genes %>% dplyr::select(V1.x, V2.x, V3.x, V11, gene_name) 
peak.genes.bed.merged <- peak.genes %>% arrange(V2.x) %>% 
  group_by(gene_name) %>% mutate(start_peak = V2.x[which.min(V2.x)], end_peak = V3.x[which.max(V3.x)])  %>% 
  mutate(V4 = ".", V5 = ".")

peak.genes.bed.merged.bed <- peak.genes.bed.merged %>% dplyr::select(V1.x, start_peak, end_peak, V4, V5, V11, gene_name) %>%
  distinct(gene_name, .keep_all = T)


