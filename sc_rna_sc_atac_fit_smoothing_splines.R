library(tidyverse)
library(openxlsx)
library(doParallel)
library(npreg)


num.cores <- detectCores(all.tests = FALSE, logical = TRUE)


# Read in the data.

sc.rna.genes.expr.pt <- readRDS('../Input_sub/toxo_cdc/rds_ME49_59/sc_rna_genes_expr_pt.rds')
sc.atac.genes.expr.pt <- readRDS('../Input_sub/toxo_cdc/rds_ME49_59/sc_atac_genes_expr_pt.rds')


## Common genes between the data sets

comm.genes <- intersect(unique(sc.rna.genes.expr.pt$GeneID), unique(sc.atac.genes.expr.pt$GeneID))


do.filt <- F

# Fit smoothing splines and sample at regular time-points

# Expression

sc.rna.spline.fits <- mclapply(1:length(comm.genes), function(i){
  tmp <- sc.rna.genes.expr.pt %>% dplyr::filter(GeneID == comm.genes[i]) %>%
    transmute(GeneID = GeneID, x = pt.shifted.scaled, y = log2.expr)
 
  if(do.filt){
    ind <- smoothFilter(tmp$y)
    tmp2 <- tmp[-ind,]
    y <- tmp2$y
    t <- tmp2$x
  }else{
    y <- tmp$y
    t <- tmp$x
    
  }
  #sc.rna.sp <- ss(t, y, periodic = T, lambda = 1e-4)
  w <- rep(1, length(y))
  w[which(y == 0)] <- 1/2
  sc.rna.sp <- smooth.spline(t, y, lambda = 0.1, w = w)
  
  sc.rna.sp <- predict(sc.rna.sp, seq(0, 6, by = 1/3)) 
  #sc.rna.pp <- fitPsplines(t, y)
  #plot(tmp$x, tmp$y)
  #points(t, y, col = 'red')
  #points(sc.rna.sp$x, sc.rna.sp$y, type = 'l', col = 'red')
  #points(sc.rna.sp$x, sc.rna.sp$y, type = 'l', col = 'blue')
  #points(sc.rna.pp$x, sc.rna.sp$y, type = 'l', col = 'green')
  mu <- data.frame(x = sc.rna.sp$x, y = sc.rna.sp$y) 
  mu <- data.frame(GeneID = rep(tmp$GeneID[1], length(mu[,1])), x = mu[,1], y = mu[,2])
  return(mu)
}, mc.cores = num.cores)

sc.rna.spline.fits <- bind_rows(sc.rna.spline.fits)



# Access

sc.atac.spline.fits <- mclapply(1:length(comm.genes), function(i){
  tmp <- sc.atac.genes.expr.pt %>% dplyr::filter(GeneID == comm.genes[i]) %>%
    transmute(GeneID = GeneID, x = pt.shifted.scaled, y = log2.expr)
  
  if(do.filt){
    ind <- smoothFilter(tmp$y)
    tmp2 <- tmp[-ind,]
    y <- tmp2$y
    t <- tmp2$x
  }else{
    y <- tmp$y
    t <- tmp$x
    
  }
  #sc.atac.sp <- ss(t, y, periodic = T, lambda = 1e-4)
  w <- rep(1, length(y))
  w[which(y == 0)] <- 1/2
  sc.atac.sp <- smooth.spline(t, y, lambda = 0.1, w = w)
  sc.atac.sp <- predict(sc.atac.sp, seq(0, 6, by = 1/3)) 
  
  mu <- data.frame(x = sc.atac.sp$x, y = sc.atac.sp$y)
  mu <- data.frame(GeneID = rep(tmp$GeneID[1], length(mu[,1])), x = mu[,1], y = mu[,2])
  return(mu)
}, mc.cores = num.cores)

sc.atac.spline.fits <- bind_rows(sc.atac.spline.fits)

