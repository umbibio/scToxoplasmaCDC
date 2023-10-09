
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

## WT 

S.O.ref <- readRDS("../Input_sub/toxo_cdc/rds_ME49_59/S.O.rna.WT_labels.rds")

##  A2XII8 KD 

S.O.rna.KD <- readRDS('../Input_sub/toxo_cdc/rds_ME49_59/S.O.rna.AP2XII8.KD.new_transferred_lables_bootroyed.rds')


## Integration of WT and KD scRNA-seq data sets

S.Os <- list(rna.WT = S.O.rna.WT, rna_KD = S.O.rna.KD)

features <- SelectIntegrationFeatures(object.list = S.Os, nfeatures = 6000)
reference_dataset <- 1
anchors <- FindIntegrationAnchors(object.list = S.Os, 
                                  anchor.features = features, reference = reference_dataset)
S.O.integrated <- IntegrateData(anchorset = anchors)
# switch to integrated assay. Make sure to set to RNA for Differential Expression
DefaultAssay(S.O.integrated) <- "integrated"
S.O.integrated <- ScaleData(object = S.O.integrated, verbose = FALSE)
S.O.integrated <- RunPCA(S.O.integrated, features = VariableFeatures(object = S.O.integrated))
S.O.integrated <- FindVariableFeatures(S.O.integrated, nfeatures = 6000)
S.O.integrated <- FindNeighbors(S.O.integrated, dims = 1:10, reduction = 'pca')
S.O.integrated <- FindClusters(S.O.integrated, resolution = 0.2)
S.O.integrated <- RunUMAP(S.O.integrated, dims = 1:13, n.components = 3)

DimPlot(S.O.integrated, split.by = "orig.ident", reduction = "umap", dims = c(1,2))


S.O.integrated@meta.data$new.spp <- gsub("\\.KD", "", S.O.integrated@meta.data$orig.ident)
S.O.integrated@meta.data$phase <- gsub("\\.", "", S.O.integrated@meta.data$phase)
S.O.integrated@meta.data$orig.ident <- gsub("\\.", "_", S.O.integrated@meta.data$orig.ident)
S.O.integrated@meta.data$phase.spp <- paste(S.O.integrated@meta.data$orig.ident, S.O.integrated@meta.data$phase, sep = ":")


S.O.integrated@reductions$pca@cell.embeddings[,1] <- -1 * S.O.integrated@reductions$pca@cell.embeddings[,1]
S.O.integrated@reductions$umap@cell.embeddings[,2] <- -1 * S.O.integrated@reductions$umap@cell.embeddings[,2]
S.O.integrated$orig.ident <- factor(S.O.integrated$orig.ident, levels = c("scRNA", "scRNA_KD"))
S.O.integrated$inferred.phase <- factor(S.O.integrated$phase, levels = c("G1a", "G1b", "S", "M", "C"))


## Add automatic Seurat clusters 

S.O.rna.WT <- readRDS('../Input_sub/toxo_cdc/rds_ME49_59/S.O.rna.WT_labels.rds')
rna.wt.met <- data.frame(cells = rownames(S.O.rna.WT@meta.data),
                         seurat_clausters_indiv = S.O.rna.WT@meta.data$seurat_clusters)
rownames(rna.wt.met) <- paste(rownames(S.O.rna.WT@meta.data), "_1", sep = "")


S.O.rna.KD <- readRDS('../Input_sub/toxo_cdc/rds_ME49_59/S.O.rna.AP2XII8.KD.new_transferred_lables_bootroyed.rds')
rna.KD.met <- data.frame(cells = rownames(S.O.rna.KD@meta.data),
                         seurat_clausters_indiv = S.O.rna.KD@meta.data$seurat_clusters)
rownames(rna.KD.met) <- paste(rownames(S.O.rna.KD@meta.data), "_2", sep = "")


meta <- rbind(rna.wt.met, rna.KD.met)
S.O.integrated <- AddMetaData(S.O.integrated, meta)
S.O.integrated@meta.data$phase.seurat.indiv <-paste(S.O.integrated@meta.data$phase , S.O.integrated@meta.data$seurat_clausters_indiv, sep = "_")



