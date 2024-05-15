
source('./util_funcs.R')
source('./loadlb.R')

num.cores <- detectCores(all.tests = FALSE, logical = TRUE)


prod.desc <- read.xlsx("../Input_sub/rds_ME49_59/Toxo_genomics/genes/MJ_annotation_07_27_2023.xlsx")
atac_sub <- readRDS('../Input_sub/rds_ME49_59/S.O_intra_atac_lables_pt.rds')
rna_sub <- readRDS('../Input_sub/rds_ME49_59/S.O_intra_lables_pt.rds')

# Differential gene expression (inferred cell cycle phase)

Idents(rna_sub) <- 'phase'
DefaultAssay(rna_sub) <- 'RNA'
Intra.markers <- FindAllMarkers(object = rna_sub, only.pos = T, min.pct = 0)

Intra.markers$GeneID <- gsub('-', '_', Intra.markers$gene)
Intra.markers.top <- Intra.markers %>% group_by(cluster) %>% top_n(2, avg_log2FC)

Intra.markers.sig <- Intra.markers %>% dplyr::filter(avg_log2FC > log2(1.5) & p_val_adj < 0.05)
Intra.markers.sig <- left_join(Intra.markers.sig, prod.desc, by = c("gene" = "TGME49") )

ss <- Intra.markers.sig %>% group_by(cluster) %>% summarise(num.DEG = n())


# Differential gene expression (inferred transition point)

prod.desc <- read.xlsx("../Input_sub/rds_ME49_59/Toxo_genomics/genes/MJ_annotation_07_27_2023.xlsx")
rna_sub <- readRDS('../Input_sub/rds_ME49_59/S.O.intra_rna_atac_trnasition_v2.rds')
atac_sub <- readRDS('../Input_sub/rds_ME49_59/S.O.intra_atac_atac_trnasition_v2.rds')


Idents(rna_sub) <- 'transition.rna'
rna.markers.rna.trans <- FindAllMarkers(object = rna_sub, only.pos = T, min.pct = 0.0)
rna.markers.rna.trans$GeneID <- gsub('-', '_', rna.markers.rna.trans$gene)

rna.sig.markers.rna.trans <- rna.markers.rna.trans %>% dplyr::filter(avg_log2FC > log2(1.5) & p_val_adj < 0.05)

rna.sig.markers.rna.trans <- left_join(rna.sig.markers.rna.trans, prod.desc, by = c("GeneID" = "TGME49") )
rna.sig.markers.rna.trans <- rna.sig.markers.rna.trans %>% 
  mutate(new.prod.desc = ifelse(ProductDescription == "hypothetical protein" & !is.na(new.name), new.name, ProductDescription))

ss <- rna.sig.markers.rna.trans %>% group_by(cluster) %>% summarise(num.DEG = n())






