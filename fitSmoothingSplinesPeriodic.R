library(tidyverse)
library(openxlsx)
library(doParallel)
library(npreg)

## For parallel calculations
num.cores <- detectCores(all.tests = FALSE, logical = TRUE)


## Read in the data.

sc.rna.genes.expr.pt <- readRDS('../input_sub/toxo_cdc/rds_ME49_59/sc_rna_genes_expr_pt.rds')
sc.atac.genes.expr.pt <- readRDS('../input_sub/toxo_cdc/rds_ME49_59/sc_atac_genes_expr_pt.rds')


## Common genes between the data sets
comm.genes <- intersect(unique(sc.rna.genes.expr.pt$GeneID), unique(sc.atac.genes.expr.pt$GeneID))



## Expression
lbx <- 1.1
sc.rna.spline.fits <- mclapply(1:length(comm.genes), function(i){
  tmp <- sc.rna.genes.expr.pt %>% dplyr::filter(GeneID == comm.genes[i]) %>%
    transmute(GeneID = GeneID, x = pt.shifted.scaled, y = log2.expr)
 
  
  y <- tmp$y
  t <- tmp$x
  
  
  ### Enforcing preiodicity
  ntimes <- length(t)
  y.ext <- c(y, y, y)
  t.ext <- c(t, t+6, t+12) 
  w.ext <- rep(1, length(y.ext))
  w.ext[which(y.ext == 0)] <- 1/3
  sc.rna.sp.ext <- smooth.spline(t.ext, y.ext, spar = lbx, w = w.ext)
  #sc.rna.sp.ext <- smooth.spline(t.ext, y.ext, cv = T, w = w.ext)
  #sparx <- sc.rna.sp.ext$spar
  #if(sparx < lbx){
  #  sc.rna.sp.ext <- smooth.spline(t.ext, y.ext, spar = lbx, w = w.ext)
  #  
  #}
  sc.rna.sp.ext <- predict(sc.rna.sp.ext, seq(6, 12, by = 1/3)) 
  sc.rna.sp <- sc.rna.sp.ext
  sc.rna.sp$x <- round(sc.rna.sp$x - 6, 4)
  #plot(t, y)
  #points(sc.rna.sp$x, sc.rna.sp$y, type = 'l', col = 'red')
  
  ######
 
  mu <- data.frame(x = sc.rna.sp$x, y = sc.rna.sp$y) 
  mu <- data.frame(GeneID = rep(tmp$GeneID[1], length(mu[,1])), x = mu[,1], y = mu[,2])
  return(mu)
}, mc.cores = num.cores)

sc.rna.spline.fits <- bind_rows(sc.rna.spline.fits)

saveRDS(sc.rna.spline.fits, '../Input_sub/toxo_cdc/rds_ME49_59/sc_rna_spline_fits_all_genes_1.1.rds')

## Access
lbx <- 1.1
sc.atac.spline.fits <- mclapply(1:length(comm.genes), function(i){
  tmp <- sc.atac.genes.expr.pt %>% dplyr::filter(GeneID == comm.genes[i]) %>%
    transmute(GeneID = GeneID, x = pt.shifted.scaled, y = log2.expr)

  
  y <- tmp$y
  t <- tmp$x
  
  ### Enforcing preiodicity
  ntimes <- length(t)
  y.ext <- c(y, y, y)
  t.ext <- c(t, t+6, t+12) 
  w.ext <- rep(1, length(y.ext))
  w.ext[which(y.ext == 0)] <- 1/3
  sc.atac.sp.ext <- smooth.spline(t.ext, y.ext, spar = lbx, w = w.ext)
  #sc.atac.sp.ext <- smooth.spline(t.ext, y.ext, cv = T, w = w.ext)
  #sparx <- sc.atac.sp.ext$spar
  #if(sparx < lbx){
  #  sc.atac.sp.ext <- smooth.spline(t.ext, y.ext, spar = lbx, w = w.ext)
  #  
  #}
  sc.atac.sp.ext <- predict(sc.atac.sp.ext, seq(6, 12, by = 1/3)) 
  sc.atac.sp <- sc.atac.sp.ext
  sc.atac.sp$x <- round(sc.atac.sp$x - 6, 4)
  ######
 
  mu <- data.frame(x = sc.atac.sp$x, y = sc.atac.sp$y)
  mu <- data.frame(GeneID = rep(tmp$GeneID[1], length(mu[,1])), x = mu[,1], y = mu[,2])
  return(mu)
}, mc.cores = num.cores)

sc.atac.spline.fits <- bind_rows(sc.atac.spline.fits)

saveRDS(sc.atac.spline.fits, '../input_sub/toxo_cdc/rds_ME49_59/sc_atac_spline_fits_all_genes_1.1.rds')

