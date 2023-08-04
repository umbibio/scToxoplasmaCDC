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


source('./util_funcs.R')

## Read scRAN-Seq data and convert TGGT1 Ids to TGME49
S.O <- readRDS('../Input_KZ/toxo_cdc/rds_ME49_59/S.O.intra_lables.rds')

## IDs
prod.desc  <- read.xlsx('../Input/toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')

## Map to ME49
counts = S.O@assays$RNA@counts
rownames(counts) <- gsub('_', '-', TGGT1_ME49$TGME49[match(gsub('-', '_', rownames(counts)), TGGT1_ME49$TGGT1)])

S.O.ME49 <- CreateSeuratObject(counts = counts)
S.O.ME49$orig.ident <- 'scRNA'
S.O.ME49 <- AddMetaData(S.O.ME49, S.O@meta.data)
Idents(S.O.ME49) <- 'phase'

S.O.ME49@meta.data$phase <- factor(S.O.ME49@meta.data$phase, levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))
S.O.ME49 <- prep_S.O(S.O.ME49)
Idents(S.O.ME49) <- 'phase'
DimPlot(S.O.ME49, reduction = 'pca')


## prepare genome

ME49.fasta <- readDNAStringSet("../Input/toxo_genomics/genome/ToxoDB-59_TgondiiME49_Genome.fasta")
chrs <- names(ME49.fasta)[grep("TGME49_chr", names(ME49.fasta))]

chr.len <- data.frame(chr = gsub(" ", "", unlist(lapply(strsplit(chrs, split = '\\|'), `[[`, 1))),
                      len = as.numeric(gsub('length=', '', unlist(lapply(strsplit(chrs, split = '\\|'), `[[`, 4)))))

txdb <- makeTxDbFromGFF(file="../Input/toxo_genomics/genome/ToxoDB-59_TgondiiME49_filter.gtf",
                        dataSource="Toxodb",
                        organism="Toxoplasma")

trans_biotypes <- select(txdb, keys=keys(txdb, "TXID"), 
                         columns = "TXTYPE", keytype =  "TXID")

genome(txdb) <- 'ME49'

tx_trans <- exonsBy(txdb, by = "tx", use.names = TRUE)
tx_names <- names(tx_trans)
num.exons <- lapply(tx_trans, function(x) length(x))
tx_names <- rep(tx_names, unlist(num.exons))
tx_trans <- unlist(tx_trans)
tx_trans$tx_id <- tx_names
tx_trans$gene_id <- gsub('-t.*', '', tx_trans$tx_id)
tx_trans$gene_name <- tx_trans$gene_id
tx_trans$type <- 'exon'
tx_trans$gene_biotype <- 'protein_coding'
tx_trans$exon_name <- tx_trans$exon_rank

tmp <- chr.len$len[match(names(seqlengths(txdb)), chr.len$chr)]
names(tmp) <- names(seqlengths(txdb))
seqlengths(tx_trans) <- tmp
seqlevels(tx_trans)
inds <- c(5,6,1,2,3,7,8,10,11,9,4,12,13,14) 



seqlevels(tx_trans) <- seqlevels(tx_trans)[inds]
isCircular(tx_trans) <- rep(F, length(isCircular(tx_trans)))


seqinfo(tx_trans)

saveRDS(tx_trans, "../Input_KZ/toxo_cdc/rds_ME49_59/ME49_tx_trans_granges.rds")

## Now read scATAC data from cellranger 
counts <- Read10X_h5(filename = "../Input/toxo_scATAC_MJ_ME49_59/filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "../Input/toxo_scATAC_MJ_ME49_59/singlecell.csv",
  header = TRUE,
  row.names = 1
)

metadata.filt <- metadata
metadata.filt$Sample <- rownames(metadata.filt)
metadata.filt <- metadata.filt[metadata.filt$Sample %in% colnames(counts), ]
peak_anno <- read_tsv("../Input_YR/toxo_scATAC_MJ_ME49_59/filtered_peak_bc_matrix/peaks.bed", col_names = c('Chr', 'strt', 'stp'))

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = seqinfo(tx_trans),
  fragments = '../Input/toxo_scATAC_MJ_ME49_59/fragments.tsv.gz',
  min.cells = 5,
  min.features = 100
)

Tg_ATAC <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata.filt
)


Tg_ATAC[['peaks']]
granges(Tg_ATAC)
annotations <- tx_trans

Annotation(Tg_ATAC) <- annotations
Tg_ATAC <- NucleosomeSignal(object = Tg_ATAC)
Tg_ATAC <- TSSEnrichment(object = Tg_ATAC, fast = FALSE)


# add blacklist ratio and fraction of reads in peaks
Tg_ATAC$pct_reads_in_peaks <- Tg_ATAC$peak_region_fragments / Tg_ATAC$passed_filters * 100
Tg_ATAC$blacklist_ratio <- Tg_ATAC$blacklist_region_fragments / Tg_ATAC$peak_region_fragments

Tg_ATAC$high.tss <- ifelse(Tg_ATAC$TSS.enrichment > 3, 'High', 'Low')
TSSPlot(Tg_ATAC, group.by = 'high.tss') + NoLegend()

Tg_ATAC$nucleosome_group <- ifelse(Tg_ATAC$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = Tg_ATAC, group.by = 'nucleosome_group') # takes a while

VlnPlot(
  object = Tg_ATAC,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)


Tg_ATAC <- subset(
  x = Tg_ATAC,
  subset = peak_region_fragments > 200 &
    peak_region_fragments < 6000 &
    pct_reads_in_peaks > 40 &
    #blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)

# stats
dim(Tg_ATAC@assays$peaks)



Tg_ATAC <- RunTFIDF(Tg_ATAC)
Tg_ATAC <- FindTopFeatures(Tg_ATAC, min.cutoff = 'q0')
Tg_ATAC <- RunSVD(Tg_ATAC)

DepthCor(Tg_ATAC)

## Must remove highly correlating components
Tg_ATAC <- RunUMAP(object = Tg_ATAC, reduction = 'lsi', dims = seq(1:30)[-c(1,3)])
Tg_ATAC <- FindNeighbors(object = Tg_ATAC, reduction = 'lsi', dims = seq(1:30)[-c(1,3)])
Tg_ATAC <- FindClusters(object = Tg_ATAC, verbose = FALSE, algorithm = 3)


DimPlot(object = Tg_ATAC, label = TRUE, reduction = 'umap') + NoLegend()

## estimate RNA levels from atac-seq 
## this is needed when we integrate with scRNA - WT which is our reference

DefaultAssay(Tg_ATAC) <- "peaks"
gene.activities <- GeneActivity(Tg_ATAC, extend.upstream = 600,
                                extend.downstream = 200)

S.O.ATAC <- CreateSeuratObject(counts = gene.activities)
S.O.ATAC$orig.ident <- 'scATAC'

S.O.ATAC <- NormalizeData(
  object = S.O.ATAC,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(S.O.ATAC$nCount_RNA)
)

saveRDS(S.O.ATAC, '../Input_KZ/toxo_cdc/rds_ME49_59/S.O_ATAC_not_integrated_not_down_samples.rds')
S.O.ATAC <- readRDS( '../Input_KZ/toxo_cdc/rds_ME49_59/S.O_ATAC_not_integrated_not_down_samples.rds')


## integrate with scRNA

DefaultAssay(S.O.ATAC) <- 'RNA'
S.O.ATAC <- FindVariableFeatures(S.O.ATAC, selection.method = "vst", nfeatures = 6000)
S.O.list <- list(RNA = S.O.ME49, ATAC = S.O.ATAC)
features <- SelectIntegrationFeatures(object.list = S.O.list, nfeatures = 6000)
reference_dataset <- 1
anchors <- FindIntegrationAnchors(object.list = S.O.list, 
                                  anchor.features = features, reference = reference_dataset)
S.O.integrated <- IntegrateData(anchorset = anchors)
# switch to integrated assay. Make sure to set to RNA for Differential Expression
DefaultAssay(S.O.integrated) <- "integrated"
S.O.integrated <- ScaleData(object = S.O.integrated, verbose = FALSE)
S.O.integrated <- RunPCA(S.O.integrated, features = VariableFeatures(object = S.O.integrated))
S.O.integrated <- FindVariableFeatures(S.O.integrated, nfeatures = 6000)
S.O.integrated <- FindNeighbors(S.O.integrated, dims = 1:10, reduction = 'pca')
S.O.integrated <- FindClusters(S.O.integrated, resolution = 0.2)
S.O.integrated <- RunUMAP(S.O.integrated, dims = 1:13)


## Transfer labels to scATAC
Idents(S.O.integrated) <- 'orig.ident'

atac_sub <- subset(S.O.integrated, ident = 'scATAC')
rna_sub <- subset(S.O.integrated, ident = 'intra')

anchors <- FindTransferAnchors(reference = rna_sub, query = atac_sub, dims = 1:30)
predictions <- TransferData(anchorset = anchors, refdata = rna_sub@meta.data$phase,dims = 1:30)
atac_sub <- AddMetaData(object = atac_sub, metadata = predictions)
atac_sub@meta.data$phase <- atac_sub@meta.data$predicted.id
Idents(atac_sub) <- 'phase'
DimPlot(atac_sub, reduction = 'pca')

ind1 <- S.O.integrated@meta.data$orig.ident == 'scATAC'
ind2 <- match(rownames(S.O.integrated@meta.data)[ind1], rownames(atac_sub@meta.data))
S.O.integrated@meta.data$phase[ind1] <- atac_sub@meta.data$phase[ind2]

ind <- S.O.integrated$orig.ident == 'intra'
S.O.integrated$orig.ident[ind] <- 'scRNA'
Idents(S.O.integrated) <- 'phase'
S.O.integrated$phase <- factor(S.O.integrated$phase, levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))
S.O.integrated@meta.data$phase <- factor(S.O.integrated@meta.data$phase, levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))
S.O.integrated@reductions[["pca"]]@cell.embeddings[,2] <- -1 * S.O.integrated@reductions[["pca"]]@cell.embeddings[,2]

p <- DimPlot(S.O.integrated, reduction = "pca", 
             #group.by = "cell", 
             split.by = 'orig.ident',
             pt.size = 1,
             #shape.by='spp',
             label = T, label.size = 5) + NoLegend() + 
  theme(panel.spacing = unit(0.5, "lines")) + 
  theme(axis.text.x = element_text(face="bold", size=12, angle=0)) +
  theme(axis.text.y = element_text(face="bold", size=12, angle=0)) +
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )


plot(p)



saveRDS(S.O.integrated, '../Input_KZ/toxo_cdc/rds_ME49_59/S.O.intra_atac_integrated.rds')

## new
Tg_ATAC[['RNA']] <- CreateAssayObject(counts = gene.activities)

Tg_ATAC <- NormalizeData(
  object = Tg_ATAC,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(Tg_ATAC$nCount_RNA)
)

DefaultAssay(Tg_ATAC) <- 'RNA'

## For PCA
DefaultAssay(Tg_ATAC) <- 'RNA'
Tg_ATAC <- FindVariableFeatures(Tg_ATAC, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(Tg_ATAC)
Tg_ATAC <- ScaleData(Tg_ATAC, features = all.genes)
Tg_ATAC <- RunPCA(Tg_ATAC, features = VariableFeatures(object = Tg_ATAC))
Tg_ATAC <- FindNeighbors(Tg_ATAC, dims = 1:10, reduction = 'pca')
Tg_ATAC <- FindClusters(Tg_ATAC, resolution = 0.2)
Tg_ATAC <- RunTSNE(object = Tg_ATAC,features = VariableFeatures(object = Tg_ATAC) )
DimPlot(object = Tg_ATAC, reduction = "tsne", label = TRUE) + NoLegend()


## Coverage Browser

Tg_ATAC <- AddMetaData(Tg_ATAC, atac_sub@meta.data)
#levels(Tg_ATAC) <- c("CD4 Naive","CD4 Memory","CD8 Naive","CD8 Effector","DN T","NK CD56bright","NK CD56Dim","pre-B",'pro-B',"pDC","DC","CD14 Mono",'CD16 Mono')
Idents(Tg_ATAC) <- 'phase'

##Find Markers
DefaultAssay(Tg_ATAC) <- 'peaks'
da_peaks <- FindAllMarkers(
  object = Tg_ATAC,
  only.pos = T,
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)

head(da_peaks)

plot1 <- VlnPlot(
  object = Tg_ATAC,
  features = rownames(da_peaks)[1],
  pt.size = 0.4,
  idents = c("G1.a","G1.b", 'S', 'M', 'C')
)

plot2 <- FeaturePlot(
  object = Tg_ATAC,
  features = rownames(da_peaks)[1],
  reduction = 'pca',
  pt.size = 0.4
)

plot1 | plot2

head(da_peaks)

top.da <- da_peaks %>% group_by(cluster) %>% slice_max(n = 1, order_by = avg_log2FC)
#region <- StringToGRanges(regions = gsub('TGME49-', 'TGME49_', top.da$gene[4]),  sep = c("-", "-"))
#region <- StringToGRanges(regions = gsub('TGME49-', 'TGME49_', rownames(da_peaks)[5]),  sep = c("-", "-"))

xx <- findOverlaps(region, tx_trans)
tx_trans[xx@to]$gene_id

my.gene <- tx_trans[xx@to]$gene_id[1]

DefaultAssay(Tg_ATAC) <- 'RNA'

DefaultAssay(atac_sub) <- "RNA"
p1 <- FeaturePlot(
  object = atac_sub,
  features = gsub('_', '-', my.gene),
  pt.size = 0.4,
  max.cutoff = 'q0',
  ncol = 1,
  reduction = 'pca'
)

plot(p1)


saveRDS(Tg_ATAC, '../Input/toxo_cdc/rds_ME49_59/S.O_ATAC_peak.rds')

##################################
##################################
##################################

Tg_ATAC <- readRDS('../Input/toxo_cdc/rds_ME49_59/S.O_ATAC_peak.rds')
Tg_ATAC@meta.data$phase2 <- gsub("\\.", "", Tg_ATAC@meta.data$phase) 

#### region_gene assignments
## Generating region/gene peak  counts
##Find Markers

peak.regions <- rownames(Tg_ATAC@assays$peaks@data)

regions <- lapply(peak.regions, function(rr){
  region <- StringToGRanges(regions = gsub('TGME49-', 'TGME49_', rr),  sep = c("-", "-"))
  #xx <- findOverlaps(region, tx_trans, maxgap = 100)
  xx <- findOverlaps(region, tx_trans, ignore.strand=TRUE)
  tx_trans[xx@to]$gene_id
  L = list(region = region, gene_id = tx_trans[xx@to]$gene_id)
})

regions <- lapply(regions, function(rr){
  if(length(rr$gene_id) == 0){
    rr$gene_id <- NA
  }
  data.frame(rr$gene_id, as.data.frame(rr$region))
})

regions <- bind_rows(regions)
colnames(regions) <- c('GeneID', 'chr', 'strt', 'stp', 'width', 'strand')

region_gene_assignment <- regions %>% 
  transmute(chr = chr, strt = strt, stp = stp, width = width, strand = strand, GeneID = GeneID)

region_gene_assignment <- region_gene_assignment %>% na.omit()



gtf.file <- "../Input/toxo_genomics/genome/ToxoDB-59_TgondiiME49_filter.gtf"
gtf <- read.table(gtf.file, header = F, sep = '\t', quote = NULL)
gtf.filt <- gtf %>% dplyr::filter(grepl('TGME49*', V1))
## filter gtf for transcripts only 
gtf.filt.trn <- gtf.filt %>% filter(V3 == "transcript")
gtf.filt.trn$gene_name <- gsub("\\;.*", "", gsub("transcript_id ", "", gsub("-t.*", "", gtf.filt.trn$V9)))
gtf.filt.trn$gene_name <- gsub("\"", "", gtf.filt.trn$gene_name)

region_gene_assignment$strand <- gtf.filt.trn$V7[match(region_gene_assignment$GeneID, gtf.filt.trn$gene_name)]
saveRDS(region_gene_assignment, '../Input/toxo_cdc/rds_ME49_59/sc_atac_regions_gene_assignment.rds')



###### Pileup Tracks
## gene.activity matrix created with other approaches can be passed here
peak.genes.bed.merged.bed <- read.xlsx("../Input/toxo_scATAC_MJ_ME49_59/peak_gene_assigned_final.xlsx")
colnames(peak.genes.bed.merged.bed) <- gsub("V1.x", "chr", colnames(peak.genes.bed.merged.bed))

peak.genes.bed.merged.bed$peak <- paste(peak.genes.bed.merged.bed$chr, 
                                         peak.genes.bed.merged.bed$start_peak,
                                         peak.genes.bed.merged.bed$end_peak, sep = "-")
gene <- "TGME49_300100"
gene <- "TGME49_204530"
gene <- "TGME49_309820"

peak <- peak.genes.bed.merged.bed$peak[peak.genes.bed.merged.bed$gene_name %in% gene]
region <- StringToGRanges(regions = gsub('TGME49-', 'TGME49_',peak),  sep = c("-", "-"))
#region <- StringToGRanges(regions = gsub('TGME49-', 'TGME49_', rownames(da_peaks)[4]),  sep = c("-", "-"))


Idents(Tg_ATAC) <- 'phase2'
Tg_ATAC@active.ident <- factor(Tg_ATAC@active.ident, levels = c('G1a', 'G1b', 'S', 'M', 'C'))
DefaultAssay(Tg_ATAC) <- 'peaks'
p2 <- CoveragePlot(
  object = Tg_ATAC,
  sep = c("-", "-"),
  region = region,
  extend.upstream = 3000,
  extend.downstream = 1000) 


p2 <- p2 & scale_fill_manual(values = c("#b6232a",'#ed7202', '#caae05', '#6f883a','#b138ee'))
p2

ggsave(filename="../Output/toxo_cdc/ME49_59/figures_paper/RPL34_track_pileup_by_phase.pdf",
       plot=p2,
       width = 5, height = 4,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

# p3 <- DimPlot(object = Tg_ATAC, label = TRUE, reduction = 'tsne') + NoLegend()
# 
# ggsave(filename="../Output/scClockFigs/tsne_clusters.pdf",
#        plot=p3,
#        width = 5, height = 5,
#        units = "in", # other options are "in", "cm", "mm"
#        dpi = 300
# )



##





singleGRange <- GRanges(as.data.frame(GRangesList(regions)))
saveRDS(as.data.frame(singleGRange), '../Input/toxo_cdc/rds/atac_regions.rds')

gtf.file <- "../Input/toxo_genomics/genome/ToxoDB-52_TgondiiME49_filter.gtf"
gtf <- read.table(gtf.file, header = F, sep = '\t', quote = NULL)
gtf.filt <- gtf %>% dplyr::filter(grepl('TGME49*', V1))
## Remove the first Exon from transcripts.
gtf.exon <- gtf.filt %>% dplyr::filter(V3 == 'exon')
gtf.exon.sort <- gtf.exon %>% arrange(V1, V4, V5)
parse.str <- strsplit(gtf.exon$V9, split = ' ')
inds <- unlist(lapply(parse.str , function(x) which(grepl("gene_id", x)) + 1))
gtf.exon$gene_name <- gsub(";", "", unlist(lapply(1:length(inds), function(i) parse.str[[i]][inds[[i]]])))
gtf.exon$gene_name <- gsub("\"", "", gtf.exon$gene_name)
gtf.exon <- gtf.exon %>% group_by(V9) %>% mutate(exon.ord = ifelse(V7 == '+', 1:n(), seq(n(), 1, by = -1)),
                                                 multiple.exon = ifelse(n() > 1, T, F))
## Remove the exon1, but keep the Intron 1
gtf.exon.2Ton <- gtf.exon %>% mutate(V10 = ifelse(multiple.exon & V7 == '-', min(V4), 
                                                  ifelse(multiple.exon & V7 == '+', min(V5), V4)),
                                     V11 = ifelse(multiple.exon & V7 == '-', max(V4), 
                                                  ifelse(multiple.exon & V7 == '+', max(V5), V5))) %>%
  mutate(V4 = V10, V5 = V11) %>% 
  dplyr::select(-c(exon.ord,multiple.exon, V10, V11) ) %>% distinct()


## filter gtf for transcripts only 
gtf.filt.trn <- gtf.filt %>% filter(V3 == "transcript")
gtf.filt.trn$gene_name <- gsub("\\;.*", "", gsub("transcript_id ", "", gsub("-t.*", "", gtf.filt.trn$V9)))
gtf.filt.trn$gene_name <- gsub("\"", "", gtf.filt.trn$gene_name)

## Filter for first exon coordinates
tmp.neg <- gtf.exon %>% filter(V7 == "-") %>% group_by(V9) %>%  dplyr::slice(which.max(V5))
tmp.pos <- gtf.exon %>% filter(V7 == "+") %>% group_by(V9) %>%  dplyr::slice(which.min(V5))
gtf.exon1 <- bind_rows(tmp.pos, tmp.neg)
gtf.exon1.sort <- gtf.exon1 %>% arrange(V1, V4, V5)

library(bedtoolsr)
ranges.tab <- as.data.frame(singleGRange)
ranges.tab <- ranges.tab %>% dplyr::filter(grepl('TGME49*', seqnames)) %>% arrange(seqnames, start, end)

peaks.genes.dist <- bedtoolsr::bt.closest(a = as.data.frame(ranges.tab), 
                                          b = as.data.frame(gtf.exon1.sort), D = "b", k = 5)

peaks.genes.dist$gene_name <- gsub("\\;.*", "", gsub("transcript_id ", "", gsub("-t.*", "", peaks.genes.dist$V13)))
#gtf.filt.trn$gene_name <- gsub("\"", "", gtf.filt.trn$gene_name)

peaks.genes.dist.trns <- left_join(peaks.genes.dist, gtf.filt.trn, by = "gene_name")

bedtoolsr::bt.closest()
region_to_gene <- findOverlaps(singleGRange, tx_trans)
subsetByOverlaps(singleGRange, tx_trans)

ov_ab <- as(findOverlaps(singleGRange, tx_trans), "List")
tx_trans[xx@to]$gene_id

my.gene <- tx_trans[xx@to]$gene_id[1]

DefaultAssay(Tg_ATAC) <- 'RNA'

p1 <- FeaturePlot(
  object = Tg_ATAC,
  features = gsub('_', '-', my.gene),
  pt.size = 0.4,
  max.cutoff = 'q0',
  ncol = 1,
  reduction = 'tsne'
)

plot(p1)

# ggsave(filename="../Output/scClockFigs/RON2_ATAC_activity.pdf",
#        plot=p1,
#        width = 5, height = 5,
#        units = "in", # other options are "in", "cm", "mm"
#        dpi = 300
# )


Idents(Tg_ATAC) <- 'phase'
Tg_ATAC@active.ident <- factor(Tg_ATAC@active.ident, levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))
DefaultAssay(Tg_ATAC) <- 'peaks'
p2 <- CoveragePlot(
  object = Tg_ATAC,
  sep = c("-", "-"),
  #region = gsub('TGME49-', 'TGME49_', rownames(Tg_ATAC)[1:3]),
  region = region,
  #region = gsub('TGME49-', 'TGME49_', VariableFeatures(Tg_ATAC)[20]),
  extend.upstream = 4000,
  extend.downstream = 2000
)

plot(p2)
ggsave(filename="../Output/scClockFigs/RON2_track_pileup_by_phase.pdf",
       plot=p2,
       width = 5, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


# p3 <- DimPlot(object = Tg_ATAC, label = TRUE, reduction = 'tsne') + NoLegend()
# 
# ggsave(filename="../Output/scClockFigs/tsne_clusters.pdf",
#        plot=p3,
#        width = 5, height = 5,
#        units = "in", # other options are "in", "cm", "mm"
#        dpi = 300
# )


saveRDS(Tg_ATAC, '../Input/scClock/Tg_ATAC.RData')
saveRDS(S.O.integrated, '../Input/scClock/S_O_scRNA_scATAC_integrated.RData')

