
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

## prep rna-WT with already transferred lables from bootroyed to use as reference for integration
## you dont need to run this just read the rds file
S.O <- readRDS('../Input_sub/toxo_cdc/rds_ME49_59/S.O.intra_lables.rds')
prod.desc  <- read.xlsx('../Input_sub/toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input_sub/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')

counts = S.O@assays$RNA@counts
rownames(counts) <- gsub('_', '-', TGGT1_ME49$TGME49[match(gsub('-', '_', rownames(counts)), TGGT1_ME49$TGGT1)])

S.O.ref <- CreateSeuratObject(counts = counts)
S.O.ref$orig.ident <- 'scRNA'
S.O.ref <- AddMetaData(S.O.ref, S.O@meta.data)
Idents(S.O.ref) <- 'phase'

S.O.ref@meta.data$phase <- factor(S.O.ref@meta.data$phase, levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))
S.O.ref <- prep_S.O(S.O.ref)
Idents(S.O.ref) <- 'phase'
ind1 <- S.O.ref@meta.data$orig.ident == 'intra'
S.O.ref$orig.ident[ind1] <- 'scRNA'

DimPlot(S.O.ref, reduction = 'pca')

saveRDS(S.O.ref, "../Input_sub/toxo_cdc/rds_ME49_59/S.O.rna.WT_labels.rds")



## individually processed rna WT
S.O.rna.WT <- readRDS("../Input_sub/toxo_cdc/rds_ME49_59/S.O.rna.WT_labels.rds")
DimPlot(S.O.rna.WT, reduction = 'umap')


umap.WT <- FetchData(object = S.O.rna.WT, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "phase"))
umap.WT$phase <- factor(umap.WT$phase)
plot_ly(umap.WT, 
        x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3,
        color =  umap.WT$phase, 
        size = 0.2)

pca.WT <- FetchData(object = S.O.rna.WT, vars = c("PC_1", "PC_2", "PC_3", "phase"))
pca.WT$phase <- factor(pca.WT$phase)
plot_ly(pca.WT, 
        x = ~PC_1, y = ~PC_2, z = ~PC_3,
        color =  pca.WT$phase, 
        size = 0.4)



##  A2XII8 KD with transferred labels and TGME49 IDs

S.O.rna.KD <- readRDS('../Input_sub/toxo_cdc/rds_ME49_59/S.O.rna.AP2XII8.KD.new_transferred_lables_bootroyed.rds')
Idents(S.O.rna.KD) <- "phase"
DimPlot(S.O.rna.KD, reduction = "pca")


umap.KD <- FetchData(object = S.O.rna.KD, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "phase", "seurat_clusters"))
umap.KD$phase<- factor(umap.KD$phase)
plot_ly(umap.KD, 
        x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3,
        color =  umap.KD$phase, 
        size = 0.2)

umap.KD$seurat_clusters<- factor(umap.KD$seurat_clusters)
plot_ly(umap.KD, 
        x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3,
        color =  umap.KD$seurat_clusters, 
        size = 0.2)


pca.KD <- FetchData(object = S.O.rna.KD, vars = c("PC_1", "PC_2", "PC_3", "phase"))
pca.KD$seurat_clusters <- factor(pca.KD$phase)
plot_ly(pca.KD, 
        x = ~PC_1, y = ~PC_2, z = ~PC_3,
        color =  pca.KD$phase, 
        size = 0.4)

Idents(S.O.rna.KD) <- "phase"
p1 <- DimPlot(S.O.rna.KD, reduction = "pca", 
              cols = c("G1.a" = "#b6232a","G1.b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee'))
p1 <- p1 + ggtitle("scRNA_KD")

p2 <- DimPlot(S.O.rna.WT, reduction = "pca",
              cols = c("G1.a" = "#b6232a","G1.b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee'))
p2 <- p2 + ggtitle("scRNA")

p <- p2|p1
p
ggsave("../Output/toxo_cdc/ME49_59/figures_paper/rna_WT_AP2XII8_new_KD_transferred_lables_bootroyed_No_integration_PCA.pdf", 
       plot = p, width = 8, height = 4, dpi = 300)

Idents(S.O.rna.KD) <- "phase"
p1 <- DimPlot(S.O.rna.KD, reduction = "umap", 
              cols = c("G1.a" = "#b6232a","G1.b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee'))
p1 <- p1 + ggtitle("scRNA_KD")

p2 <- DimPlot(S.O.rna.WT, reduction = "umap",
              cols = c("G1.a" = "#b6232a","G1.b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee'))
p2 <- p2 + ggtitle("scRNA_KD")

p <- p2|p1
p
ggsave("../Output/toxo_cdc/ME49_59/figures_paper/rna_WT_AP2XII8_new_KD_transferred_lables_bootroyed_No_integration_UMAP.pdf", 
       plot = p, width = 8, height = 4, dpi = 300)

#####################################################
#################### Integration ########################
#####################################################

# this is done just to have one seurat object for DE analysis,  we should set the assay to RNA for DE 

S.O.rna.WT <- readRDS("../Input_sub/toxo_cdc/rds_ME49_59/S.O.rna.WT_labels.rds")
S.O.rna.KD <- readRDS('../Input_sub/toxo_cdc/rds_ME49_59/S.O.rna.AP2XII8.KD.new_transferred_lables_bootroyed.rds')


S.Os <- list(rna.WT = S.O.rna.WT, rna_KD = S.O.rna.KD)

## Integration
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



## better naming of data sets and phases

S.O.integrated@meta.data$new.spp <- gsub("\\.KD", "", S.O.integrated@meta.data$orig.ident)
S.O.integrated@meta.data$phase <- gsub("\\.", "", S.O.integrated@meta.data$phase)
S.O.integrated@meta.data$orig.ident <- gsub("\\.", "_", S.O.integrated@meta.data$orig.ident)
S.O.integrated@meta.data$phase.spp <- paste(S.O.integrated@meta.data$orig.ident, S.O.integrated@meta.data$phase, sep = ":")


S.O.integrated@reductions$pca@cell.embeddings[,1] <- -1 * S.O.integrated@reductions$pca@cell.embeddings[,1]
S.O.integrated@reductions$umap@cell.embeddings[,2] <- -1 * S.O.integrated@reductions$umap@cell.embeddings[,2]
S.O.integrated$orig.ident <- factor(S.O.integrated$orig.ident, levels = c("scRNA", "scRNA_KD"))
S.O.integrated$inferred.phase <- factor(S.O.integrated$phase, levels = c("G1a", "G1b", "S", "M", "C"))


## add seurat clusters obtained from individual processing of data sets 

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


saveRDS(S.O.integrated, "../Input_sub/toxo_cdc/rds_ME49_59/S.O.integrated_rna_WT_AP2XII8KD_reference_rna_WT_transferred_lables_from_boot.rds")

# test plot
S.O.integrated <- readRDS("../Input_sub/toxo_cdc/rds_ME49_59/S.O.integrated_rna_WT_AP2XII8KD_reference_rna_WT_transferred_lables_from_boot.rds")

## plot 
Idents(S.O.integrated) <- "seurat_clausters_indiv"
Idents(S.O.integrated) <- "phase"
p <- DimPlot(S.O.integrated, reduction = "pca", dims = c(1,2), 
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
    axis.title.y = element_text(size=14, face="bold"), legend.title = element_text(size = 8)
  ) 

p


