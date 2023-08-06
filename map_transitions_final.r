library(tidyverse)
library(openxlsx)
library(Seurat)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(tidytext)
library(parallel)


## For parallel calculations
num.cores <- detectCores(all.tests = FALSE, logical = TRUE)

getCurvePeakLoc <- function(t, y, prob = 0.8){
  
  ## Fitting the estimated kernel with smooth splines
  spline.fit <- smooth.spline(x = t, y = y)
  
  
  ## Compute the derivatives of the fitted splines
  s.0 <- predict(spline.fit, spline.fit$x, deriv=0)
  s.1 <- predict(spline.fit, spline.fit$x, deriv=1)
  s.derv <- data.frame(s0=s.0$y, s1=s.1$y)
  
  ## Get the location of the extrema
  locs <- rle(den.sign <- sign(s.derv$s1))
  
  
  ## Maxima
  inc.ind <- which(locs$values == 1)
  if(length(inc.ind) > 1){
    maxima.ind = {}
    for(i in inc.ind){
      maxima.ind = c(maxima.ind, sum(locs$lengths[1:i]))
    }
    ## Interpolate a point between the location where derivative changes sign
    maxima = (spline.fit$x[maxima.ind] + spline.fit$x[(maxima.ind + 1)]) / 2
    maxima = maxima[!is.na(maxima)]
    ## Get the maximum values
    maxval = predict(spline.fit, maxima)
    
    ## Get the outliers
    maxima.outliers = which(maxval$y >= quantile(maxval$y, prob = prob))
    
    ## Peaks for entities of interest
    entity.x = maxval$x[maxima.outliers]
    entity.y = maxval$y[maxima.outliers]
  }else{
    entity.x <- spline.fit$x[which.max(spline.fit$y)]
    entity.y <- spline.fit$y[which.max(spline.fit$y)]
  }
  
  return(entity.x)
}

## Read in the data.
## IDs
prod.desc  <- read.xlsx('../Input_sub/toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input_sub/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')
prod.desc <- left_join(prod.desc, TGGT1_ME49, by = c('GeneID' = 'TGGT1'))

marker.genes <- readRDS('../Input_sub/toxo_cdc/rds_ME49_59/Intra_markers_sig.rds')
marker.genes.phase <- marker.genes %>% transmute(GeneID = gene, phase = cluster) %>% distinct()


sc.rna.spline.fits <- readRDS('../Input_sub/toxo_cdc/rds_ME49_59/sc_rna_spline_fits_all_genes.rds')
sc.atac.spline.fits <- readRDS('../Input_sub/toxo_cdc/rds_ME49_59/sc_atac_spline_fits_all_genes.rds')

## Filter to include markers only
sc.rna.spline.fits <- sc.rna.spline.fits %>% dplyr::filter(GeneID %in% marker.genes$gene)
sc.atac.spline.fits <- sc.atac.spline.fits %>% dplyr::filter(GeneID %in% marker.genes$gene)


rna_sub <- readRDS('../Input_sub/toxo_cdc/rds_ME49_59/S.O_intra_lables_pt.rds')
atac_sub <- readRDS('../Input_sub/toxo_cdc/rds_ME49_59/S.O_intra_atac_lables_pt.rds')


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


sc.rna.peak.order <- sc.rna.mu.scale %>% group_by(GeneID) %>% summarise(peak.ord = getCurvePeakLoc(x, expr))
sc.rna.mu.scale <- left_join(sc.rna.mu.scale, sc.rna.peak.order, by = 'GeneID')

sc.rna.mu.scale$GeneID <- factor(sc.rna.mu.scale$GeneID, 
                                 levels = unique(sc.rna.mu.scale$GeneID[order(-sc.rna.mu.scale$peak.ord)]))


sc.atac.peak.order <- sc.atac.mu.scale %>% group_by(GeneID) %>% summarise(peak.ord = getCurvePeakLoc(x, expr))
sc.atac.mu.scale <- left_join(sc.atac.mu.scale, sc.atac.peak.order, by = 'GeneID')

sc.atac.mu.scale$GeneID <- factor(sc.atac.mu.scale$GeneID, 
                                  levels = unique(sc.atac.mu.scale$GeneID[order(-sc.atac.mu.scale$peak.ord)]))

## Filter to include markers only
sc.rna.mu.scale <- sc.rna.mu.scale %>% dplyr::filter(GeneID %in% marker.genes$gene)
sc.atac.mu.scale <- sc.atac.mu.scale %>% dplyr::filter(GeneID %in% marker.genes$gene)

sc.rna.mu.scale <- sc.rna.mu.scale %>% mutate(t.t = case_when(peak.ord >= 0 & peak.ord < 3 ~ 'G1',
                                                              peak.ord >= 3 & peak.ord < 4.7 ~ 'S',
                                                              peak.ord >= 4.7 & peak.ord < 5 ~ 'M',
                                                              peak.ord >= 5 ~ 'C'))

sc.rna.mu.scale$t.t <- factor(sc.rna.mu.scale$t.t, levels = c('G1', 'S', 'M', 'C'))

sc.atac.mu.scale <- sc.atac.mu.scale %>% mutate(t.t = case_when(peak.ord >= 0 & peak.ord < 3 ~ 'G1',
                                                                peak.ord >= 3 & peak.ord < 4.7 ~ 'S',
                                                                peak.ord >= 4.7 & peak.ord < 5 ~ 'M',
                                                                peak.ord >= 5 ~ 'C'))

sc.atac.mu.scale$t.t <- factor(sc.atac.mu.scale$t.t, levels = c('G1', 'S', 'M', 'C'))


tans.points <- function(mu.scale, lam = 0.1, prob = 0.1){
  peak.locs <- mu.scale %>% dplyr::select(GeneID, peak.ord) %>% distinct()
  spline.fit.peaks <- smooth.spline(x = 1:nrow(peak.locs), y = sort(peak.locs$peak.ord), lambda = lam)
  
  s.0 <- predict(spline.fit.peaks, seq(1, nrow(peak.locs), by = 0.1), deriv=0)
  s.1 <- predict(spline.fit.peaks, seq(1, nrow(peak.locs), by = 0.1), deriv=1)
  s.2 <- predict(spline.fit.peaks, seq(1, nrow(peak.locs), by = 0.1), deriv=2)
  
  spline.fit.peaks.smooth <- data.frame(x = seq(1, nrow(peak.locs), by = 0.1), 
                                        s0=s.0$y, s1=s.1$y, s2 = s.2$y)
  
  ## Get locations of peaks and valies
  max.loc <- getCurvePeakLoc(spline.fit.peaks.smooth$x, spline.fit.peaks.smooth$s1, prob = prob)
  min.loc <- getCurvePeakLoc(spline.fit.peaks.smooth$x, -spline.fit.peaks.smooth$s1, prob = prob)
  
  transition.points <- predict(spline.fit.peaks, sort(c(max.loc, min.loc)))
  
  L <- list(peak.locs = peak.locs, spline.fit.peaks.smooth = spline.fit.peaks.smooth, 
            transition.points = transition.points)
  
  return(L)
  
}



## ATAC transition points
L.trans.atac <- tans.points(sc.atac.mu.scale, lam = 0.005, prob = 0.0)

## set transition to be between 0 and 6 hrs
L.trans.atac$transition.points$y[L.trans.atac$transition.points$y < 0] <- 0
L.trans.atac$transition.points$y[L.trans.atac$transition.points$y > 6] <- 6


## removing nosy transition (manual)
L.trans.atac$transition.points$y <- L.trans.atac$transition.points$y[-c(3,5,7)]
L.trans.atac$transition.points$x <- L.trans.atac$transition.points$x[-c(3,5,7)]

saveRDS(L.trans.atac, '../Input_sub/toxo_cdc/rds_ME49_59/atac_based_transition_points_v2.rds')

## plot the spline fitted to peak time of expression (only marker genes)
## add the transition points as vertical lines
## fig 3A

atac_peaks.dat <- L.trans.atac$spline.fit.peaks.smooth %>% 
  transmute(g = x, y = s0, yp = s1) %>% pivot_longer(-g, names_to = 'drivs', values_to = 'value')
atac_peaks.dat$value[atac_peaks.dat$value < 0] <- 0
atac_peaks.dat$value[atac_peaks.dat$value > 6] <- 6
atac_peaks.dat$drivs <- factor(atac_peaks.dat$drivs, levels = c('y', 'yp'))


p1  <- ggplot(atac_peaks.dat, aes(x= g,y=value)) +
  geom_path(aes(color = drivs),alpha = 0.8, size = 1)+ 
  theme_bw(base_size = 14) +
  geom_vline(xintercept=L.trans.atac$transition.points$x, linetype=2, color = 'orange', size = 0.6) + 
  ylab('peak time') + xlab('genes') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  theme(strip.text = element_text(size = 14, face="bold", angle = 0)) + 
  
  
  facet_grid(drivs~., scales = 'free') +
  
  #ggtitle(titles[i]) +
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


plot(p1)



## RNA transition points
L.trans.rna <- tans.points(sc.rna.mu.scale, lam = 0.005, prob = 0.0)
L.trans.rna$transition.points$y[L.trans.rna$transition.points$y < 0] <- 0
L.trans.rna$transition.points$y[L.trans.rna$transition.points$y > 6] <- 6

## removing nosy transition (manual)
L.trans.rna$transition.points$y <- L.trans.rna$transition.points$y[-2]
L.trans.rna$transition.points$x <- L.trans.rna$transition.points$x[-2]

saveRDS(L.trans.rna, '../Input_sub/toxo_cdc/rds_ME49_59/rna_based_transition_points_v2.rds')

## plot the spline fitted to peak time of accessibility (only marker genes)
## add the transition points as vertical lines

rna_peaks.dat <- L.trans.rna$spline.fit.peaks.smooth %>% 
  transmute(g = x, y = s0, yp = s1) %>% pivot_longer(-g, names_to = 'drivs', values_to = 'value')
rna_peaks.dat$value[rna_peaks.dat$value < 0] <- 0
rna_peaks.dat$value[rna_peaks.dat$value > 6] <- 6


rna_peaks.dat$drivs <- factor(rna_peaks.dat$drivs, levels = c('y', 'yp'))
p2  <- ggplot(rna_peaks.dat, aes(x= g,y=value)) +
  geom_path(aes(color = drivs),alpha = 0.8, size = 1)+ 
  theme_bw(base_size = 14) +
  geom_vline(xintercept=L.trans.rna$transition.points$x, linetype=2, color = 'orange', size = 0.6) + 
  ylab('peak time') + xlab('genes') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  theme(strip.text = element_text(size = 14, face="bold", angle = 0)) + 
  
  
  facet_grid(drivs~., scales = 'free') +
  
  #ggtitle(titles[i]) +
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


plot(p2)


## Add transition time  to the meta data (according to transition time and cell cycle phase)

mapTransToPCA <- function(x, y, p){
  z = rep('T0', length(x))
  z[x < y[1]] <- paste('T', length(y), sep = '')
  for(i in 1:(length(y) - 1)){
    z[x >= y[i] & x < y[(i+1)]] <- paste('T', i, sep = '')
  }
  
  z[x >= y[length(y)]] <- paste('T', length(y), sep = '')
  z[p == 'C' & z == 'T1'] <- paste('T', length(y), sep = '')
  
  return(z)
  
}

## 0 and 6 are start of cell cycle 
## here I am using 0 as the start of cell cycle for rna transitions

L.trans.rna$transition.points$y <- c(0, L.trans.rna$transition.points$y[-4])

t.atac.over.rna <- mapTransToPCA(rna_sub@meta.data$pt.shifted.scaled, L.trans.atac$transition.points$y, rna_sub@meta.data$phase)
rna_sub@meta.data$transition.atac <- factor(t.atac.over.rna, levels = sort(unique(t.atac.over.rna)))

t.rna.over.rna <- mapTransToPCA(rna_sub@meta.data$pt.shifted.scaled, L.trans.rna$transition.points$y, rna_sub@meta.data$phase)
rna_sub@meta.data$transition.rna <- factor(t.rna.over.rna, levels = sort(unique(t.rna.over.rna)))

t.atac.over.atac <- mapTransToPCA(atac_sub@meta.data$pt.shifted.scaled, L.trans.atac$transition.points$y, atac_sub@meta.data$phase)
atac_sub@meta.data$transition.atac <- factor(t.atac.over.atac, levels = sort(unique(t.atac.over.atac)))



saveRDS(rna_sub, '../Input_sub/toxo_cdc/rds_ME49_59/S.O.intra_rna_atac_trnasition_v2.rds')
saveRDS(atac_sub, '../Input/toxo_cdc/rds_ME49_59/S.O.intra_atac_atac_trnasition_v2.rds')


## plot for test
rna_sub@reductions[["pca"]]@cell.embeddings[,2] <- -1 * rna_sub@reductions[["pca"]]@cell.embeddings[,2]

Idents(rna_sub) <- 'phase'
p1 <- DimPlot(rna_sub, reduction = 'pca', label = T) + NoLegend()
p1

Idents(rna_sub) <- 'transition.atac'
p2 <- DimPlot(rna_sub, reduction = 'pca', label = T) + NoLegend()
p2

Idents(rna_sub) <- 'transition.rna'
p3 <- DimPlot(rna_sub, reduction = 'pca', label = T) + NoLegend()
p3

pp <- p1|p2|p3
pp


## rna
Idents(rna_sub) <- 'transition.rna'
p4 <- DimPlot(rna_sub, reduction = 'pca', dims = c(2,3), label = T) + NoLegend()
p5 <- p1|p2|p3|p4


p <- p1|p3

## atac
atac_sub@reductions[["pca"]]@cell.embeddings[,2] <- -1 * atac_sub@reductions[["pca"]]@cell.embeddings[,2]

Idents(atac_sub) <- 'phase'
p1 <- DimPlot(atac_sub, reduction = 'pca', label = T) + NoLegend()
p1

Idents(atac_sub) <- 'transition.atac'
p2 <- DimPlot(atac_sub, reduction = 'pca', label = T) + NoLegend()
p2

p1 | p2 



