
source('./util_funcs.R')
source('./loadlb.R')

num.cores <- detectCores(all.tests = FALSE, logical = TRUE)


# Count files
intra.file.csv <- "../Input_sub/toxo_scRNA_MJ_ME49_59/RH.intra.expr.csv"

## IDs
prod.desc  <- read.xlsx('../Input_sub/toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input_sub/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')


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


# create Seurat objects
feats <- c("nFeature_RNA","nCount_RNA")

S.O.intra <- CreateSeuratObject(counts = intra.counts)
S.O.intra$orig.ident <- 'intra'
VlnPlot(S.O.intra, features = feats, pt.size = 0.1,ncol = 2) + NoLegend()
FeatureScatter(S.O.intra, "nCount_RNA", "nFeature_RNA",  pt.size = 0.5)
selected_c <- WhichCells(S.O.intra, expression = nFeature_RNA > 100 & nFeature_RNA < 950)
selected_f <- rownames(S.O.intra)[ Matrix::rowSums(S.O.intra) > 5]
S.O.intra  <- subset(S.O.intra, features=selected_f, cells=selected_c)
dim(S.O.intra@assays$RNA@data)

 
set.seed(100)
S.O <- subset(x = S.O.intra, downsample = 8000)
S.O <- prep_S.O(S.O, res = 0.4)
DimPlot(S.O)


# Transfer cell cycle phase lables from available data https://elifesciences.org/articles/54129

S.O.list <- list(intra = S.O.intra)

S.O.tg.boothroyd <- readRDS('../Input_sub/toxo_cdc/rds_ME49_59/S.O.tg_RH_boothroyd.rds')

S.Os <- mclapply(S.O.list, function(S.O){
  S.O <- prep_S.O(S.O, res = 0.4)
  anchors <- FindTransferAnchors(reference = S.O.tg.boothroyd, query = S.O, dims = 1:30)
  predictions <- TransferData(anchorset = anchors, refdata = S.O.tg.boothroyd@meta.data$phase,dims = 1:30)
  predictions$phase <- predictions$predicted.id
  S.O <- AddMetaData(object = S.O, metadata = predictions)
  return(S.O)
}, mc.cores = num.cores)

spps <- names(S.Os)

S.Os <- lapply(1:length(S.Os), function(i){
  S.Os[[i]]@meta.data$spp <- spps[i]
  S.Os[[i]]
})


saveRDS(S.Os[[1]], '../Input_sub/toxo_cdc/rds_ME49_59/S.O.intra_lables.rds')

