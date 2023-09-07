library(Seurat)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(gam)
library(princurve)
library(parallel)
library(tidyverse)
library(MyEllipsefit)
library(sctransform)
library(openxlsx)
library(doParallel)
library(tidytext)
library(ggrepel)
library(dtwclust)
library(geomtextpath)


source('./util_funcs.R')
num.cores <- detectCores(all.tests = FALSE, logical = TRUE)

## Splines
sc.rna.spline.fits <- readRDS('../Input_sub/toxo_cdc/rds_ME49_59/sc_rna_spline_fits_all_genes.rds')
sc.atac.spline.fits <- readRDS('../Input_sub/toxo_cdc/rds_ME49_59/sc_atac_spline_fits_all_genes.rds')

# cyclic.genes <- readRDS('../Input_sub/toxo_cdc/rds_ME49_59/all_genes_cyclic_timing.rds')
# cyclic.genes <- cyclic.genes %>% dplyr::filter(rna.cyclic == 1 & atac.cyclic == 1)


marker.genes <- readRDS('../Input_sub/toxo_cdc/rds_ME49_59/Intra_markers_sig.rds')
marker.genes.phase <- marker.genes %>% transmute(GeneID = gene, phase = cluster) %>% distinct()


## Turn the data into wide format (time by gene) and center & scale each gene
sc.rna.dtw.wide <- sc.rna.spline.fits %>% 
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  mutate_at(vars(matches("TGME")), scale) %>%
  as.data.frame()

sc.atac.dtw.wide <- sc.atac.spline.fits %>% 
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  mutate_at(vars(matches("TGME")), scale) %>%
  as.data.frame()


sc.rna.mu.scale <- sc.rna.dtw.wide %>% 
  pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')

sc.atac.mu.scale <- sc.atac.dtw.wide %>% 
  pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')


## Filter to include markers only
sc.rna.mu.scale <- sc.rna.mu.scale %>% dplyr::filter(GeneID %in% marker.genes$gene)
sc.atac.mu.scale <- sc.atac.mu.scale %>% dplyr::filter(GeneID %in% marker.genes$gene)

# ## Filter to include cyclic
# sc.rna.mu.scale <- sc.rna.mu.scale %>% dplyr::filter(GeneID %in% cyclic.genes$GeneID)
# sc.atac.mu.scale <- sc.atac.mu.scale %>% dplyr::filter(GeneID %in% cyclic.genes$GeneID)
# 

genes <- unique(sc.atac.mu.scale$GeneID)

cc.dat <- mclapply(1:length(genes), function(i){
  atac.g <- sc.atac.mu.scale %>% dplyr::filter(GeneID == genes[i]) %>% 
    transmute(t = x, y = expr) %>% arrange(t)
  rna.g <- sc.rna.mu.scale %>% dplyr::filter(GeneID == genes[i]) %>% 
    transmute(t = x, y = expr) %>% arrange(t)
  tmp <- ccf(c(atac.g$y), c(rna.g$y), plot = F)
  L <- list(Lag = tmp$lag[which.max(tmp$acf)], cc = max(tmp$acf))
  return(L)
}, mc.cores = num.cores)

cc.dat <- data.frame(GeneID = genes, Lags = unlist(lapply(cc.dat, `[[`, 1)), ccs = unlist(lapply(cc.dat, `[[`, 2)))

#saveRDS(cc.dat, '../Input_KZ//toxo_cdc/rds_ME49_59/sc_rna_sc_atac_cross_cor_lag.rds')
p <- ggplot(cc.dat, aes(x = ccs)) + 
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "white") +
  geom_density(lwd = 1,
               linetype = 1,
               colour = 2)+
  theme_bw() + 
  theme(plot.title = element_text(face = "bold.italic", size = 18),
        axis.title = element_text(face = "bold", size = 14)) +
  ggtitle("scRNA & scATAC cross-correlation")

p

sum(cc.dat$ccs > 0.6) / nrow(cc.dat)

## Correlation to gene families

## IDs
prod.desc  <- read.xlsx('../Input_sub//toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input_sub/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')

## New
gene.fam <- read.xlsx("../Input_sub/Toxo_genomics/gene_families/gene_fam_KZ.xlsx")
gene.fam <- left_join(gene.fam, TGGT1_ME49, by = c('GeneID' = 'TGGT1'))
gene.fam$cc <- cc.dat$ccs[match(gsub('_', '-', gene.fam$TGME49), cc.dat$GeneID)]

gene.fam <- gene.fam %>% na.omit()

fam.ord <- gene.fam %>% group_by(Family) %>% summarise(m = median(cc))
gene.fam$Family <- factor(gene.fam$Family, levels = fam.ord$Family[sort(fam.ord$m, index.return = T)$ix])
p <- ggplot(data = gene.fam, aes(x = Family, y = cc, color = Family)) + 
  geom_boxplot() + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) + 
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  theme_bw(base_size = 14) +
  theme(legend.position = "right") +
  ylab('curve cross correlation') + xlab('Family') +
  ylim(c(0, 1.01)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  theme(strip.text = element_text(size = 14, face="bold", angle = 0)) + 
  theme(
    plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=14, face="bold", hjust = 1),
    axis.title.y = element_text(size=14, face="bold")
  ) + 
  theme(#legend.position = c(0.15, 0.85),
    legend.position = 'none',
    legend.title = element_text(colour="black", size=12, 
                                face="bold"),
    legend.text = element_text(colour="black", size=12, 
                               face="bold"))


p


