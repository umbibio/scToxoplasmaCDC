library(Seurat)
library(openxlsx)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(matrixStats)
library(tidyverse)
library(RColorBrewer)
library(parallel)
library(openxlsx)
library(plotly)


source('./util_funcs.R')

num.cores <- detectCores(all.tests = FALSE, logical = TRUE)

## Count files AP2XII-8 KD
intra.file.csv <- "../Input_sub/toxo_scRNA_AP2XII8_KD_230327/AP2XII8_KD.expr.csv"

## IDs
prod.desc  <- read.xlsx('../input_sub/toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../input_sub/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')


# expression mtx with TGGT1 IDs (needed later for transfering lables from public data- Bootroyed)
getExpr <- function(in.file, TGGT1_ME49){
  file.counts <- read.csv(in.file)
  genes <- file.counts$X
  ind <- which(genes %in% TGGT1_ME49$TGME49)
  file.counts <- file.counts[ind, ]
  genes <- genes[ind]
  genes <- TGGT1_ME49$TGGT1[match(genes, TGGT1_ME49$TGME49)]
  expr <- file.counts[,-1]
  rownames(expr) <- genes
  
  return(expr)
}

intra.counts <- getExpr(intra.file.csv, TGGT1_ME49)
dim(intra.counts)


## individual Seurat objects
feats <- c("nFeature_RNA","nCount_RNA")

# Intra
S.O.intra <- CreateSeuratObject(counts = intra.counts)
S.O.intra$orig.ident <- 'intra'
VlnPlot(S.O.intra, features = feats, pt.size = 0.1,ncol = 2) + NoLegend()
FeatureScatter(S.O.intra, "nCount_RNA", "nFeature_RNA",  pt.size = 0.5)
selected_c <- WhichCells(S.O.intra, expression = nFeature_RNA > 100 & nFeature_RNA < 1100)
selected_f <- rownames(S.O.intra)[ Matrix::rowSums(S.O.intra) > 5]
S.O.intra  <- subset(S.O.intra, features=selected_f, cells=selected_c)
dim(S.O.intra@assays$RNA@data)

S.O.list <- list(intra = S.O.intra)

## Downsample to 8000 cells
set.seed(100)
S.O.list <- mclapply(S.O.list, function(S.O){
  S.O <- subset(x = S.O, downsample = 8000)
}, mc.cores = num.cores)



## transfer labels from Bootroyed

S.O.tg.boothroyd <- readRDS('../input_sub/toxo_cdc/rds_ME49_59/S.O.tg_RH_boothroyd.rds')

## split the data, process each, transfer the lables (here we have only one data)
S.Os <- mclapply(S.O.list, function(S.O){
  S.O <- prep_S.O(S.O, res = 0.4)
  anchors <- FindTransferAnchors(reference = S.O.tg.boothroyd, query = S.O, dims = 1:30)
  predictions <- TransferData(anchorset = anchors, refdata = S.O.tg.boothroyd@meta.data$phase,dims = 1:30)
  predictions$phase <- predictions$predicted.id
  #predictions$phase[which(predictions$prediction.score.max < 0.7)] <- 'NA'
  S.O <- AddMetaData(object = S.O, metadata = predictions)
  return(S.O)
}, mc.cores = num.cores)

spps <- names(S.Os)

S.Os <- lapply(1:length(S.Os), function(i){
  S.Os[[i]]@meta.data$spp <- spps[i]
  S.Os[[i]]
})

#  plots for testing

DimPlot(S.Os[[1]], reduction = "pca")
DimPlot(S.Os[[1]], reduction = "umap")


saveRDS(S.Os[[1]], '../input_sub/toxo_cdc/rds_ME49_59/S.O.rna.AP2XII8.KD.new_transferred_lables_bootroyed_TGGT1.rds')


# convert TGGT1 IDs to ME49
S.O <- readRDS('../input_sub/toxo_cdc/rds_ME49_59/S.O.rna.AP2XII8.KD.new_transferred_lables_bootroyed_TGGT1.rds')
prod.desc  <- read.xlsx('../input_sub/toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../input_sub/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')

counts = S.O@assays$RNA@counts
rownames(counts) <- gsub('_', '-', TGGT1_ME49$TGME49[match(gsub('-', '_', rownames(counts)), TGGT1_ME49$TGGT1)])

S.O.KD <- CreateSeuratObject(counts = counts)
S.O.KD <- AddMetaData(S.O.KD, S.O@meta.data)
S.O.KD@meta.data$phase <- factor(S.O.KD@meta.data$phase, levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))
S.O.KD <- prep_S.O(S.O.KD,  res = 0.4)
ind1 <- S.O.KD@meta.data$orig.ident == 'intra'
S.O.KD$orig.ident[ind1] <- 'scRNA.KD'

Idents(S.O.KD) <- 'phase'
Idents(S.O.KD) <- 'seurat_clusters'
DimPlot(S.O.KD, reduction = 'umap')
DimPlot(S.O.KD, reduction = 'pca')

S.O.KD@reductions$pca@cell.embeddings[,2] <- -1 * S.O.KD@reductions$pca@cell.embeddings[,2]
S.O.KD@reductions$umap@cell.embeddings[,2] <- -1 * S.O.KD@reductions$umap@cell.embeddings[,2]

saveRDS(S.O.KD, '../input_sub/toxo_cdc/rds_ME49_59/S.O.rna.AP2XII8.KD.new_transferred_lables_bootroyed.rds')

pca = S.O.KD[["pca"]]
eigValues = (pca@stdev)^2  ## EigenValues
varExplained = eigValues / sum(eigValues)
sum(varExplained)


## WT - convert TGGT1 IDs to ME49
## rna-WT with already transferred lables from bootroyed will be used as reference 
## for integration of WT and KD

S.O <- readRDS('../input_sub/toxo_cdc/rds_ME49_59/S.O.intra_lables.rds')
prod.desc  <- read.xlsx('../input_sub/toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../input_sub/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')

counts = S.O@assays$RNA@counts
rownames(counts) <- gsub('_', '-', TGGT1_ME49$TGME49[match(gsub('-', '_', rownames(counts)), TGGT1_ME49$TGGT1)])

S.O.ref <- CreateSeuratObject(counts = counts)
S.O.ref$orig.ident <- 'scRNA'
S.O.ref <- AddMetaData(S.O.ref, S.O@meta.data)

S.O.ref@meta.data$phase <- factor(S.O.ref@meta.data$phase, levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))
S.O.ref <- prep_S.O(S.O.ref, res = 0.4)
Idents(S.O.ref) <- 'phase'
ind1 <- S.O.ref@meta.data$orig.ident == 'intra'
S.O.ref$orig.ident[ind1] <- 'scRNA'

Idents(S.O.ref) <- "phase"
DimPlot(S.O.ref, reduction = 'umap',
        cols = c("G1.a" = "#b6232a","G1.b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee'))

#S.O.ref@reductions[["umap"]]@cell.embeddings[,2] <- -S.O.ref@reductions[["umap"]]@cell.embeddings[,2]


saveRDS(S.O.ref, "../input_sub/toxo_cdc/rds_ME49_59/S.O.rna.WT_labels.rds")


