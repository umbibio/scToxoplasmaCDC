library(openxlsx)
library(tidyverse)
library(splines)
library(parallel)
library(ggplot2)
library(tidytext)
library(ggrepel)
library(geomtextpath)
library(ggVennDiagram)


source('./util_funcs.R')

num.cores <- detectCores(all.tests = FALSE, logical = TRUE)


getCurvePeakLoc <- function(t, y, prob = 0.8){
  
  ## Fitting the estimated kernel with smooth splines
  spline.fit <- smooth.spline(x = t, y = y)
  
  
  ## Compute the derivatives of the fitted splines
  s.0 <- predict(spline.fit, spline.fit$x, deriv=0)
  s.1 <- predict(spline.fit, spline.fit$x, deriv=1)
  s.derv <- data.frame(s0=round(s.0$y, digits = 3), s1=round(s.1$y, digits = 3))
  
  ## Get the location of the extrema
  locs <- rle(den.sign <- sign(s.derv$s1))
  
  
  ## Maxima
  inc.ind <- which(locs$values == 1)
  if(length(inc.ind) >= 1 & any(locs$values== -1)){
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
  if(length(entity.x) == 0){
    entity.x <- spline.fit$x[which.max(spline.fit$y)]
    entity.y <- spline.fit$y[which.max(spline.fit$y)]
  }
  
  return(entity.x)
}

prod.disc  <- read.xlsx('../Input_sub/rds_ME49_59/toxo_genomics/genes/ProductDescription_GT1.xlsx')
GT1.ME49 <- read.xlsx('../Input_sub/rds_ME49_59/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')
prod.disc <- left_join(prod.disc, GT1.ME49, by = c("GeneID" = "TGGT1")) %>% na.omit()
prod.disc$TGME49 <- gsub('_', '-', prod.disc$TGME49)


# sc.rna.genes.expr.pt <- readRDS('../Input_sub/toxo_cdc/rds_ME49_59/sc_rna_spline_fits_all_genes.rds')
# sc.atac.genes.expr.pt <- readRDS('../Input_sub/toxo_cdc/rds_ME49_59/sc_atac_spline_fits_all_genes.rds')

sc.rna.genes.expr.pt <- readRDS('../Input_sub/rds_ME49_59/sc_rna_spline_fits_all_genes_1.1.rds')
sc.atac.genes.expr.pt <- readRDS('../Input_sub/rds_ME49_59/sc_atac_spline_fits_all_genes_1.1.rds')

## Phase-based DEGs
Intra.markers.sig <- readRDS('../Input_sub/rds_ME49_59/Intra_markers_sig.rds')



genes <- unique(sc.rna.genes.expr.pt$GeneID)
fft.stas <- mclapply(1:length(genes), function(i){
  
  rna.sp <- sc.rna.genes.expr.pt %>% dplyr::filter(GeneID == genes[i]) %>%
    transmute(GeneID = GeneID, x = x, y = y) %>% arrange(x)

  atac.sp <- sc.atac.genes.expr.pt %>% dplyr::filter(GeneID == genes[i]) %>%
    transmute(GeneID = GeneID, x = x, y = y) %>% arrange(x)
  
  
  rna.upper.expr <- mean(rna.sp$y[rna.sp$y > quantile(rna.sp$y, p = 0.75)])
  
  rna.transform = fft(rna.sp$y)/length(rna.sp$y)
  rna.magTransform = abs(rna.transform)
  
  rna.amp <- max(rna.magTransform[-1]) * 2
  rna.freq <- which.max(rna.magTransform[-1])
  rna.peak.time <- getCurvePeakLoc(rna.sp$x, rna.sp$y, prob = 0.8)
  
  atac.upper.expr <- mean(atac.sp$y[atac.sp$y > quantile(atac.sp$y, p = 0.75)])
  
  atac.transform = fft(atac.sp$y)/length(atac.sp$y)
  atac.magTransform = abs(atac.transform)
  
  atac.amp <- max(atac.magTransform[-1]) * 2
  atac.freq <- which.max(atac.magTransform[-1])
  atac.peak.time <- getCurvePeakLoc(atac.sp$x, atac.sp$y, prob = 0.8)
  
  L <- list(rna.upper.expr, rna.amp, rna.freq, rna.peak.time,
            atac.upper.expr, atac.amp, atac.freq, atac.peak.time)
}, mc.cores = num.cores)

stats <- data.frame(GeneID = genes, 
                    rna.upper.expr = unlist(lapply(fft.stas, `[[`, 1)),
                    rna.amp = unlist(lapply(fft.stas, `[[`, 2)),
                    rna.freq = unlist(lapply(fft.stas, `[[`, 3)),
                    rna.peak.time = unlist(lapply(fft.stas, `[[`, 4)),
                    atac.upper.expr = unlist(lapply(fft.stas, `[[`, 5)),
                    atac.amp = unlist(lapply(fft.stas, `[[`, 6)),
                    atac.freq = unlist(lapply(fft.stas, `[[`, 7)),
                    atac.peak.time = unlist(lapply(fft.stas, `[[`, 8)))


stats <- left_join(stats, prod.disc, by = c('GeneID' = 'TGME49'))

rna.expr.cutoff <- quantile(stats$rna.upper.expr, prob = 0.01, na.rm = T)
stats$rna.expressed <- ifelse(stats$rna.upper.expr >= rna.expr.cutoff, 1, 0)

atac.expr.cutoff <- quantile(stats$atac.upper.expr, prob = 0.01, na.rm = T)
stats$atac.expressed <- ifelse(stats$atac.upper.expr >= atac.expr.cutoff, 1, 0)

rna.amp.cutoff <- quantile(stats$rna.amp, prob = 0.5)
stats$rna.cyclic<- ifelse(stats$rna.amp >= rna.amp.cutoff , 1, 0)

atac.amp.cutoff <- quantile(stats$atac.amp, prob = 0.5)
stats$atac.cyclic<- ifelse(stats$atac.amp >= atac.amp.cutoff , 1, 0)

stats.expressed <- stats %>% dplyr::filter(rna.expressed == 1)
nrow(stats.expressed) / nrow(stats)
nrow(stats.expressed)

stats.cyclic.rna <- stats %>% dplyr::filter(rna.expressed == 1, rna.cyclic == 1)
nrow(stats.cyclic.rna) / nrow(stats)
nrow(stats.cyclic.rna) 

stats.cyclic.atac <- stats %>% dplyr::filter(rna.expressed == 1, atac.cyclic == 1)
nrow(stats.cyclic.atac)

stats.cyclic.both <- stats %>% dplyr::filter(rna.expressed == 1, rna.cyclic == 1 & atac.cyclic == 1)
nrow(stats.cyclic.both) / nrow(stats)
nrow(stats.cyclic.both)

sum(stats$rna.cyclic) / nrow(stats)





