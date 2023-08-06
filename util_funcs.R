## functions
##
library(doParallel)
num.cores <- detectCores(all.tests = FALSE, logical = TRUE)

prep_S.O <- function(S.O, res = 0.1, var.features = F, down.sample = F){
  set.seed(100)
  S.O <- NormalizeData(S.O, normalization.method = "LogNormalize", scale.factor = 10000)
  S.O <- FindVariableFeatures(S.O, selection.method = "vst", nfeatures = 6000)
  if(var.features){
    ## Work on variable features only
    S.O <- subset(S.O, features = VariableFeatures(S.O))
  }
  all.genes <- rownames(S.O)
  S.O <- ScaleData(S.O, features = all.genes)
  S.O <- RunPCA(S.O, features = VariableFeatures(object = S.O))
  S.O <- FindNeighbors(S.O, dims = 1:10, reduction = 'pca')
  S.O <- FindClusters(S.O, resolution = res)
  if(down.sample){
    S.O <- subset(x = S.O, downsample = 800)
    S.O <- FindNeighbors(S.O, dims = 1:13)
    S.O <- FindClusters(S.O, resolution = res)
  }
  S.O <- RunUMAP(S.O, dims = 1:13, n.components = 3L)
  return(S.O)
}


smooth.S.O <- function(S.O, network = T){
  S.O@meta.data$smooth.id <- rownames(S.O@meta.data)
  if(network){
    cat('network smoothing \n')
    AdjacencyMat <- as.matrix(S.O@graphs$RNA_nn)
  
    con <- colSums(AdjacencyMat)
    con <- unlist(lapply(con, function(x) ifelse(x == 0, 0, 1/x)))
    
    ## Scaling adjacancy
    for(i in 1:ncol(AdjacencyMat)){
      AdjacencyMat[,i] <- con[i] * AdjacencyMat[,i]
    }
    
    ## Smoothing for a few iterations
    max.smoothing <- 3
    alpha <- 0.3
    scDat <- as.matrix(S.O[["RNA"]]@data)
    Ft <- scDat
    for(i in 1:max.smoothing){
      cat(paste('smoothing iteration', i, '/', max.smoothing))
      cat('\n')
      Ft <- alpha * Ft %*% AdjacencyMat + (1 - alpha) * scDat
    }
    exprs <- Ft
  }else{
    Idents(S.O) <- 'smooth.id'
    exprs <- AverageExpression(S.O, features = all.genes)
    exprs <- exprs$RNA
  }
  smooth_assay <- CreateAssayObject(counts = exprs)
  S.O[["smooth"]] <- smooth_assay
  #Idents(S.O) <- 'seurat_clusters'
  #DefaultAssay(S.O) <- "smooth"
  
  return(S.O)
}


scale.S.O <- function(S.O){
  S.O@meta.data$scale.id <- rownames(S.O@meta.data)
  counts <- S.O[["RNA"]]@data
  c.n <- colnames(counts)
  r.n <- rownames(counts)
  counts <- t(as.matrix(counts)) %>% as_tibble() %>% 
    mutate_all(~ (.x  - min(.x))/(max(.x) - min(.x))) %>% as.matrix() %>% t()
  colnames(counts) <- c.n
  rownames(counts) <- r.n
  
  scale_assay <- CreateAssayObject(counts = counts)
  S.O[["scale"]] <- scale_assay
  
  return(S.O)
}

fitPseudoTime <- function(S.O, reverset.time = F){

  num.cores <- detectCores(all.tests = FALSE, logical = TRUE)
  ## Fit a pseudo-time curve and align using sync data
  pc <- getPCA(S.O)
  sds.data <- getPrinCurve(pc)
  pc.sds <- left_join(pc, sds.data, by = "Sample")

  Y <- log2(S.O@assays$RNA@data + 1)
  var.genes <- names(sort(apply(Y, 1, var),decreasing = TRUE))#[1:1000] 
  Y <- Y[var.genes, ]
  
  pt <- sds.data$pt
  
  if(reverset.time){
    pt <- max(pt) - pt
  }
  ## Map the pseudo-time to 0-12:20 hours 
  t <- (12 + 1/3) * ((as.numeric(pt) - min(as.numeric(pt)))/(max(as.numeric(pt)) - min(as.numeric(pt))))
  sds.data$t <- t
  
  ## time-index cells in 20 min intervals and identify cells in each partition
  ## They will be considered as replicates
  time.breaks <- seq(1/3, 12 + 1/3, by = 1/3) 
  time.idx <- rep(0, nrow(sds.data))
  
  ind <- which(sds.data$t <= time.breaks[1])
  time.idx[ind] <- 0
  
  for(i in 2:length(time.breaks)){
    ind <- which(sds.data$t > time.breaks[(i-1)] & sds.data$t <= time.breaks[i])
    time.idx[ind] <- i - 1
  }
  
  sds.data$time.idx <- time.idx
  
  ## Update the time to 20 min increments
  sds.data$t <- (time.idx) * (1/3)
  
  sds.data <- sds.data %>%  
    group_by(time.idx) %>% mutate(rep = seq(1:n()))
  
  
  ## Run a GAM regression of expression on the pseudo-time
  ## Use parallel computation to speed things up. 16 cores
  cat('Fitting the GAM model\n')
  gam.pval <- mclapply(1:nrow(Y), function(z){
    d <- data.frame(z=as.numeric(Y[z,]), t=as.numeric(pt))
    tmp <- gam(z ~ lo(t), data=d)
    # p <- summary(tmp)[4][[1]][1,5] ## Linear Effects
    p <- summary(tmp)$anova$`Pr(F)`[2] ## nonlinear effects
    p
  }, mc.cores = num.cores)
  
  gam.pval <- unlist(gam.pval)
  names(gam.pval) <- rownames(Y)
  ## Remove the NA's and get the best fits
  if(any(is.na(gam.pval))){
    gam.pval <- gam.pval[-which(is.na(gam.pval))]
  }
  
  gam.pval.adj <- p.adjust(gam.pval, method = 'fdr', n = length(gam.pval))
  gam.pval.sig <- gam.pval[gam.pval.adj < 0.01] 
  print(length(gam.pval.sig)) ## number of correlating genes
  
  ## Sort the cells on the pt
  cell.ord <- sds.data$cell.ord
  
  topgenes <- names(sort(gam.pval.sig, decreasing = FALSE))  
  cell.cycle.genes.expr <- as.matrix(S.O@assays$RNA@data[topgenes, cell.ord])
  #cell.cycle.genes.expr <- as.matrix(S.O.bd.filt@assays$smooth@data[topgenes, cell.ord]) ## smoothed version
  
  
  cell.cycle.genes.df <- data.frame(GeneID = rownames(cell.cycle.genes.expr),
                                    cell.cycle.genes.expr) %>% 
    pivot_longer(-c(GeneID), names_to = 'Sample', values_to = 'log2.expr')
  
  
  cell.cycle.genes.df$GeneID <- gsub('-', '_', cell.cycle.genes.df$GeneID)
  cell.cycle.genes.df <- left_join(cell.cycle.genes.df, sds.data, by = 'Sample')
  cell.cycle.genes.df$cluster <- S.O@meta.data$seurat_clusters[match(cell.cycle.genes.df$Sample, 
                                                                        rownames(S.O@meta.data))]
 
  ## Filtering to include genes that fit well with pseudo time
  ##cat('Network smoothing on GAM genes\n')
  ##S.O.gam <- subset(S.O, features = names(gam.pval.sig))
  ##S.O.gam <- prep_S.O(S.O.gam)
  ##S.O.gam.smooth <- smooth.S.O(S.O.gam)
  

  L <- list(cell.cycle.genes.df = cell.cycle.genes.df, 
            sds.data = sds.data,
            gam.genes = names(gam.pval.sig))
  
  return(L)
}

alignWithBdiv <- function(S.O.b,  b.cell.cycle.genes.df, bd.cell.cycle.genes.df, b.sds.data, 
                          bd.sds.data, bd.lag.time = NA){
  
  num.cores <- detectCores(all.tests = FALSE, logical = TRUE)
  
  ## As a first pass, fit smoothing splines to both sync and sc and align the splines
  b.sc.tc.df <- b.cell.cycle.genes.df %>% 
    transmute(y = log2.expr, tme = t, ind = rep, variable = GeneID)
  
  bd.sc.tc.df <- bd.cell.cycle.genes.df %>% 
    transmute(y = log2.expr, tme = t, ind = rep, variable = GeneID)

  ## Get the common genes
  comm.genes <- unique(b.sc.tc.df$variable)[which(unique(b.sc.tc.df$variable) %in%
                                                    unique(bd.sc.tc.df$variable))]
  ## Fit smoothing splines to both and sample at regular time points (every 20 min from 0 - 12h)
  b.mu.sc.com <- mclapply(comm.genes, function(v){
    xx <- b.sc.tc.df[b.sc.tc.df$variable == v, c("y","tme","ind")]
    mu <-  smooth.spline(x = xx$tme, y = xx$y)
    mu
  }, mc.cores = num.cores)
  
  b.mu.sc.com.grid <- lapply(b.mu.sc.com, function(mu) predict(mu, seq(0, 12, by = 1/3)))
  
  bd.mu.sc.com <- mclapply(comm.genes, function(v){
    xx <- bd.sc.tc.df[bd.sc.tc.df$variable == v, c("y","tme","ind")]
    mu <-  smooth.spline(x = xx$tme, y = xx$y)
    mu
  }, mc.cores = num.cores)
  
  bd.mu.sc.com.grid <- lapply(bd.mu.sc.com, function(mu) predict(mu, seq(0, 12, by = 1/3)))
  
  
  ## Calculate the cross-correlation between the fitted smoothing splines
  cc.sc.genes <- mclapply(c(1:length(comm.genes)), function(i){
    ccc <- rep(0, length(b.mu.sc.com.grid[[i]]$y))
    for (tau in 0:(length(b.mu.sc.com.grid[[i]]$y) - 1)){
      circ.ind <- (0:(length(bd.mu.sc.com.grid[[i]]$y) - 1) + tau) %% length(bd.mu.sc.com.grid[[i]]$y) + 1
      ccc[tau+1] <- sum(b.mu.sc.com.grid[[i]]$y * bd.mu.sc.com.grid[[i]]$y[circ.ind])
    }
    ll <- which.max(ccc)
    ##ll <- ccf(mu.sc.com.grid[[i]]$y, mu.sync.com.grid[[i]]$y, plot = F, lag.max = length(mu.sc.com.grid[[i]]$y))
    ##ccc <- rep(0, length(mu.sc.com.grid[[i]]$y))
    ##pos.lag <- which(ll$lag > 0)
    ##ll <- ll$lag[pos.lag][which.max(ll$acf[pos.lag])]
    #ll <- ll$lag[which.max(ll$acf)]
    ll
  }, mc.cores = num.cores)
  
  
  # 
  # ## Calculate the cross-correlation between the fitted smoothing splines
  # cc.sc.genes <- mclapply(c(1:length(comm.genes)), function(i){
  #   ll <- ccf(b.mu.sc.com.grid[[i]]$y, bd.mu.sc.com.grid[[i]]$y, plot = F, lag.max = length(b.mu.sc.com.grid[[i]]$y))
  #   ll <- ll$lag[which.max(ll$acf)]
  #   ll
  # }, mc.cores = num.cores)
  # 
  
  # Histogram with density plot
  
  dd <- data.frame(lag = unlist(cc.sc.genes))
  
  ## calculate the optimal lag time
  den <- density(dd$lag)
  
  lag.time <- (ceiling(den$x[which.max(den$y)]) + bd.lag.time - 1) %% length(bd.mu.sc.com.grid[[i]]$y) + 1
  
  adjusted.time <- (b.sds.data$time.idx * 1/3) -  sort(unique(b.sds.data$time.idx) * 1/3)[lag.time]
  neg.ind <- ifelse(adjusted.time < 0, T, F)
  adjusted.time[neg.ind] <- adjusted.time[neg.ind] + (12 + 1/3)
  b.sds.data$adj.time <-  adjusted.time
  
  ## create fine-resolution 20 min cluster of cells
  clusters <- paste('C', 1:length(unique(b.sds.data$adj.time)), sep = '')
  b.sds.data$cluster <- clusters[as.integer((b.sds.data$adj.time) * 3 + 1)]
  
  
  ## Generate shifted curves
  time.breaks <- seq(1/3, 12 + 1/3, by = 1/3) 
  time.idx <- rep(0, nrow(b.sds.data))
  
  ind <- which(b.sds.data$adj.time <= time.breaks[1])
  time.idx[ind] <- 0
  
  for(i in 2:length(time.breaks)){
    ind <- which(b.sds.data$adj.time > time.breaks[(i-1)] & b.sds.data$adj.time <= time.breaks[i])
    time.idx[ind] <- i - 1
  }
  
  b.sds.data$adj.time.idx <- time.idx
  
  
  b.sds.data <- b.sds.data %>%  ungroup() %>%
    group_by(adj.time.idx) %>% mutate(rep = seq(1:n()))
  
  b.sds.data <- as.data.frame(b.sds.data)
  rownames(b.sds.data) <- b.sds.data$Sample
  
  ## Add the new clusters as meta-data
  S.O.b <- AddMetaData(S.O.b, b.sds.data)
  
  pc <- S.O.b[['pca']]@cell.embeddings
  pc <- data.frame(pc) %>% dplyr::mutate(Sample = rownames(pc))
  pc$cluster <- S.O.b$cluster
  
  
  
  pc.sds.adj <- left_join(pc, b.sds.data, by = "Sample")
  
  lvs <- paste('C', unique(sort(as.numeric(gsub('C', '', pc.sds.adj$cluster.y)))), sep = '')
  pc.sds.adj$cluster.y <- factor(pc.sds.adj$cluster.y, levels = lvs)
  
  
  b.cell.cycle.genes.df.adj <- left_join(b.cell.cycle.genes.df, b.sds.data[,c('Sample', 'adj.time', 
                                                                              'adj.time.idx', 'rep', 'cluster')], 
                                          by = 'Sample')
  b.sc.tc.df.adj <- b.cell.cycle.genes.df.adj %>% 
    transmute(y = log2.expr, tme = adj.time, ind = rep.y, variable = GeneID)
  
  
  L <- list(S.O.bd.update = S.O.b,
            mu.sync.com.grid = NULL,
            mu.sc.com.grid = b.mu.sc.com.grid,
            pc.sds.adj = pc.sds.adj,
            dd = dd,
            den = den,
            lag.time = lag.time,
            cell.cycle.genes.df.adj = b.cell.cycle.genes.df.adj,
            sync.tc.df = NULL,
            sc.tc.df.adj = b.sc.tc.df.adj)
  
  return(L)
}

alignBdivPseudoTimeWithBulkSync <- function(S.O.bd, bd.cell.cycle.genes.df, sds.data, bd.tc.logCPM, lag.time = NA){
  
  num.cores <- detectCores(all.tests = FALSE, logical = TRUE)
  
  ## As a first pass, fit smoothing splines to both sync and sc and align the splines
  sc.tc.df <- bd.cell.cycle.genes.df %>% 
    transmute(y = log2.expr, tme = t, ind = rep, variable = GeneID)
  
  sync.tc.df <- bd.tc.logCPM %>% 
    transmute(y = as.numeric(expr), 
              tme = as.numeric(time_point), 
              ind = rep, variable = GeneID)
  
  ## Get the common genes
  comm.genes <- unique(sync.tc.df$variable)[which(unique(sync.tc.df$variable) %in%
                                                    unique(sc.tc.df$variable))]
  ## Fit smoothing splines to both and sample at regular time points (every 20 min from 0 - 12h)
  mu.sync.com <- mclapply(comm.genes, function(v){
    xx <- sync.tc.df[sync.tc.df$variable == v, c("y","tme","ind")]
    mu <-  smooth.spline(x = xx$tme, y = xx$y)
    mu
  }, mc.cores = num.cores)
  
  mu.sync.com.grid <- lapply(mu.sync.com, function(mu) predict(mu, seq(0, 12, by = 1/3)))
  
  mu.sc.com <- mclapply(comm.genes, function(v){
    xx <- sc.tc.df[sc.tc.df$variable == v, c("y","tme","ind")]
    mu <-  smooth.spline(x = xx$tme, y = xx$y)
    mu
  }, mc.cores = num.cores)
  
  mu.sc.com.grid <- lapply(mu.sc.com, function(mu) predict(mu, seq(0, 12, by = 1/3)))
  
  ## Calculate the cross-correlation between the fitted smoothing splines
  cc.sc.sync.genes <- mclapply(c(1:length(comm.genes)), function(i){
    # ccc <- rep(0, length(mu.sc.com.grid[[i]]$y))
    # for (tau in 0:(length(mu.sc.com.grid[[i]]$y) - 1)){
    #   circ.ind <- (0:(length(mu.sync.com.grid[[i]]$y) - 1) + tau) %% length(mu.sync.com.grid[[i]]$y) + 1
    #   ccc[tau+1] <- sum(mu.sc.com.grid[[i]]$y * mu.sync.com.grid[[i]]$y[circ.ind])
    # }
    # ll <- which.max(ccc)
    ll <- ccf(mu.sc.com.grid[[i]]$y, mu.sync.com.grid[[i]]$y, plot = F, lag.max = length(mu.sc.com.grid[[i]]$y))
    pos.ind <- which(ll$lag >= 0)
    ll <- ll$lag[pos.ind][which.max(ll$acf[pos.ind])]
    ll
  }, mc.cores = num.cores)
  
  
  # Histogram with density plot
  
  dd <- data.frame(lag = unlist(cc.sc.sync.genes))
  
  ## calculate the optimal lag time
  den <- density(dd$lag)
  
  
  ## If already passed as argument, do not calculate
  if(is.na(lag.time)){
    lag.time <- ceiling(den$x[which.max(den$y)])
  }

  adjusted.time <- (sds.data$time.idx * 1/3) -  sort(unique(sds.data$time.idx) * 1/3)[lag.time]
  neg.ind <- ifelse(adjusted.time < 0, T, F)
  adjusted.time[neg.ind] <- adjusted.time[neg.ind] + (12 + 1/3)
  sds.data$adj.time <-  adjusted.time
  
  ## create fine-resolution 20 min cluster of cells
  clusters <- paste('C', 1:length(unique(sds.data$adj.time)), sep = '')
  sds.data$cluster <- clusters[as.integer((sds.data$adj.time) * 3 + 1)]
  
  
  ## Generate shifted curves
  time.breaks <- seq(1/3, 12 + 1/3, by = 1/3) 
  time.idx <- rep(0, nrow(sds.data))
  
  ind <- which(sds.data$adj.time <= time.breaks[1])
  time.idx[ind] <- 0
  
  for(i in 2:length(time.breaks)){
    ind <- which(sds.data$adj.time > time.breaks[(i-1)] & sds.data$adj.time <= time.breaks[i])
    time.idx[ind] <- i - 1
  }
  
  sds.data$adj.time.idx <- time.idx
  
  
  sds.data <- sds.data %>%  ungroup() %>%
    group_by(adj.time.idx) %>% mutate(rep = seq(1:n()))
  
  sds.data <- as.data.frame(sds.data)
  rownames(sds.data) <- sds.data$Sample
  
  ## Add the new clusters as meta-data
  S.O.bd <- AddMetaData(S.O.bd, sds.data)
  
  pc <- S.O.bd[['pca']]@cell.embeddings
  pc <- data.frame(pc) %>% dplyr::mutate(Sample = rownames(pc))
  pc$cluster <- S.O.bd$cluster
  
  
  
  pc.sds.adj <- left_join(pc, sds.data, by = "Sample")
  
  lvs <- paste('C', unique(sort(as.numeric(gsub('C', '', pc.sds.adj$cluster.y)))), sep = '')
  pc.sds.adj$cluster.y <- factor(pc.sds.adj$cluster.y, levels = lvs)
  

  bd.cell.cycle.genes.df.adj <- left_join(bd.cell.cycle.genes.df, sds.data[,c('Sample', 'adj.time', 
                                                                              'adj.time.idx', 'rep', 'cluster')], 
                                          by = 'Sample')
  
  bd.cell.cycle.genes.df.adj <- left_join(bd.cell.cycle.genes.df.adj, pc, by = "Sample")
  sc.tc.df.adj <- bd.cell.cycle.genes.df.adj %>% 
    transmute(y = log2.expr, tme = adj.time, ind = rep.y, variable = GeneID)

  L <- list(S.O.bd.update = S.O.bd,
            mu.sync.com.grid = mu.sync.com.grid,
            mu.sc.com.grid = mu.sc.com.grid,
            pc.sds.adj = pc.sds.adj,
            dd = dd,
            den = den,
            lag.time = lag.time,
            cell.cycle.genes.df.adj = bd.cell.cycle.genes.df.adj,
            sync.tc.df = sync.tc.df,
            sc.tc.df.adj = sc.tc.df.adj)
  
  return(L)
}

getNormExpr <- function(S.O){
  expr.norm <- as.matrix(S.O[["RNA"]]@data)
  rownames(expr.norm) <- gsub('-', '_',rownames(expr.norm))
  expr.norm <- as.data.frame(expr.norm) 
  expr.norm <- expr.norm %>% dplyr::mutate(GeneID = rownames(expr.norm))
  expr.norm <- expr.norm %>% gather(key = Sample, value = expr, -GeneID)
  expr.norm <- expr.norm  %>% group_by(GeneID) %>% 
    dplyr::mutate(quantile = ifelse(expr <= quantile(expr)[2], 1, 
                                    ifelse(expr <= quantile(expr)[3], 2, 
                                           ifelse(expr <= quantile(expr)[4], 3,4)))) %>% ungroup()
  return(expr.norm)
}


getUmap <- function(S.O){
  umap <- S.O[['umap']]@cell.embeddings
  umap <- data.frame(umap) %>% dplyr::mutate(Sample = rownames(umap)) 
  #umap <- left_join(umap, S.O@meta.data, by = 'Sample')
  umap$cluster <- S.O$seurat_clusters
  
  return(umap)
}

getPCA <- function(S.O){
  pc <- S.O[['pca']]@cell.embeddings
  pc <- data.frame(pc) %>% dplyr::mutate(Sample = rownames(pc)) 
  #pc <- left_join(pc, S.O@meta.data, by = 'Sample')
  pc$cluster <- S.O$seurat_clusters
  
  return(pc)
}

getPrinCurve <- function(pc.dat){
  ## Initial elliptical fit
  ell <- Ellipsefit(pc.dat,PC_1, PC_2, coords = TRUE, bbox = TRUE)
  coord <- ell$Coord
  
  x <- as.matrix(pc.dat[,c(1,2)])
  #fit <- principal_curve(x, start = as.matrix(coord))
  fit <- principal_curve(x, start = as.matrix(coord), smoother = "periodic_lowess", maxit = 0)
  pt <- fit$lambda
  
  ## reversing the order of time
  pt <- max(pt) - pt
  pt <- (pt - min(pt))/(max(pt) - min(pt))
  cell.ord <- fit$ord[seq(length(fit$ord), 1, by = -1)]
  
  
  s.data <- data.frame(Sample = rownames(fit$s), 
                       cell.ord = cell.ord,
                       pt = pt,
                       sc1 = fit$s[,1], 
                       sc2 = fit$s[,2])
  
  return(s.data)
}

wiskerPlot <- function(S.O){
  pc.dat <- getPCA(S.O)
  ## Initial elliptical fit
  ell <- Ellipsefit(pc.dat,PC_1, PC_2, coords = TRUE, bbox = TRUE)
  coord <- ell$Coord
  
  x <- as.matrix(pc.dat[,c(1,2)])
  #fit <- principal_curve(x, start = as.matrix(coord))
  fit <- principal_curve(x, start = as.matrix(coord), smoother = "periodic_lowess", maxit = 0)

  return(list(pc = pc.dat, fit = fit))
}

getSlingShot <- function(S.O, method = 'pca'){
  sds <- slingshot(Embeddings(S.O, method), 
                   clusterLabels = S.O$seurat_clusters, 
                   start.clus = 0, end.clus = 3, stretch = 2)
  
  pt <- slingPseudotime(sds)
  
  s.data <- data.frame(cell.ord = sds@curves$curve1$ord, pt = pt, sds@curves$curve1$s[,1:2])
  s.data$Sample <- rownames(s.data)
  s.data <- s.data %>% transmute(Sample = Sample, cell.ord = cell.ord, 
                                 pt = curve1, sc1 = PC_1, sc2 = PC_2)
  
  #s.data <- left_join(s.data, S.O.combined@meta.data, by = 'Sample')
  
  return(s.data)
}


cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}


fitSingleSme <-function(my.tc, v){
  tryCatch(
    expr = {
      fit <- sme(my.tc[my.tc$variable==v,c("y","tme","ind")], criteria = 'AIC')
      return(fit)
    },
    error = function(v){
      message(paste('error:', v))
    },
    warning = function(w){
      message(paste('error:', v))
    },
    finally = {
      fit <- sme(my.tc[my.tc$variable==v,c("y","tme","ind")])
      return(fit)
    }
  )    
}

## spline the fitted values to get the means
splineSmeFits <- function(fits, variables){
  mus <- lapply(fits, function(x) spline(x = as.numeric(colnames(x$coefficients)),
                                         y = x$coefficients[1,], 
                                         n = as.numeric(colnames(x$coefficients))[ncol(x$coefficients)] - 
                                           as.numeric(colnames(x$coefficients))[1] + 1, 
                                         method = "natural")) 

  mus <- lapply(fits, function(x) spline(x = as.numeric(colnames(x$coefficients)),
                                         y = x$coefficients[1,], 
                                         method = "natural")) 
  
  mus.y <- unlist(lapply(mus, `[[`, 2))
  mus.x <- unlist(lapply(mus, `[[`, 1))
  lens  <- unlist(lapply(lapply(mus, `[[`, 1), length))
  
  mus <- data.frame(variable = rep(variables, times = lens),
                    tme = mus.x,
                    y = mus.y)
  
  return(mus)
  
}

## spline the fitted values to get the means
splinefunctionSmeFits <- function(fits, variables, extend = F){
  musfun <- lapply(fits, function(x) splinefun(x = as.numeric(colnames(x$coefficients)),
                                         y = x$coefficients[1,], 
                                         method = "natural")) 
  if(extend){
    mus <- lapply(musfun, function(ff){
      x <- seq(0, 12, by = 1/3)
      y <- ff(x)
      mu = list(x = x, y = y)
    })
  }else{
    mus <- lapply(musfun, function(ff){
      x <- as.numeric(colnames(x$coefficients))
      y <- ff(x)
      mu = list(x = x, y = y)
    })
  }
  
  mus.y <- unlist(lapply(mus, `[[`, 2))
  mus.x <- unlist(lapply(mus, `[[`, 1))
  lens  <- unlist(lapply(lapply(mus, `[[`, 1), length))
  
  mus <- data.frame(variable = rep(variables, times = lens),
                    tme = mus.x,
                    y = mus.y)
  
  return(mus)
  
}


## spline the fitted values to get the means
smoothSplineSmeFits <- function(fits, variables, extend = F){
  ## Fitting the estimated kernel with smooth splines
  
  
  mus <- lapply(fits, function(x) 
    smooth.spline(x = as.numeric(colnames(x$coefficients)), 
                  y = x$coefficients[1,])) 
  if(extend){
    mus <- lapply(mus, function(x)
      predict(x, seq(0, 12, by = 1/3)))
  }
  
  mus.y <- unlist(lapply(mus, `[[`, 2))
  mus.x <- unlist(lapply(mus, `[[`, 1))
  lens  <- unlist(lapply(lapply(mus, `[[`, 1), length))
  
  mus <- data.frame(variable = rep(variables, times = lens),
                    tme = mus.x,
                    y = mus.y)
  
  return(mus)
  
}


## spline the fitted values to get the means
smoothSplineSmeFitsConfBand <- function(fits, variables){
  ## Fitting the estimated kernel with smooth splines
  
  mus <- lapply(fits, function(x) 
    smooth.spline(x = as.numeric(colnames(x$coefficients)), 
                  y = x$coefficients[1,])) 
  # if(extend){
  #   mus <- lapply(mus, function(x)
  #     predict(x, seq(0, 12, by = 1/3)))
  # }
  
  mus.y <- unlist(lapply(mus, `[[`, 2))
  mus.x <- unlist(lapply(mus, `[[`, 1))
  lens  <- unlist(lapply(lapply(mus, `[[`, 1), length))
  
  
  ## Fitting confidence bands
  mu.variances <- mclapply(fits, function(fit){
    tryCatch(
      expr = {
        mu.variance <- diag(vcov(fit))
        mu.variance[which(mu.variance < 0)] <- 0
        upper.band <- smooth.spline(x = as.numeric(colnames(fit$coefficients)), 
                                    y = fit$coefficients[1, ] + 1.96 * sqrt(mu.variance)) 
        lower.band <- smooth.spline(x = as.numeric(colnames(fit$coefficients)),
                                    y = fit$coefficients[1, ] - 1.96 * sqrt(mu.variance))
        return(list(lb = lower.band, ub = upper.band))
      },
      error = function(v){
        message(paste('error:', v))
      }
    ) 
    
  }, mc.cores = num.cores) 
  #err.ind <- which(unlist(mclapply(mu.variances, function(x) is.null(x$lb),  mc.cores = num.cores)))
  
  lbs <- unlist(lapply(lapply(mu.variances, `[[`, 1), `[[`, 2))
  ubs <- unlist(lapply(lapply(mu.variances, `[[`, 2), `[[`, 2))
  
  mus <- data.frame(variable = rep(variables, times = lens),
                    tme = mus.x,
                    y = mus.y,
                    lb = lbs,
                    ub = ubs)
  
   
  return(mus)
  
}

plot.sme <-function(fit, v, conf = T){
  mu <- spline(x = as.numeric(colnames(fit$coefficients)), 
               y = fit$coefficients[1,], n = 500, 
               method = "natural")
  fs <- lapply(2:nrow(fit$coefficients), function(i) {
    spline(x = as.numeric(colnames(fit$coefficients)), 
           y = fit$coefficients[1,] + fit$coefficients[i, ], 
           method = "natural", 
           n = 500)})
  
  
  ylim <- range(fit$data$y, mu$y, sapply(fs, "[[", "y"))
  xlim <- range(as.numeric(colnames(fit$coefficients)))
  mu.variance <- diag(vcov(fit))
  if(conf){
    upper.band <- spline(x = as.numeric(colnames(fit$coefficients)), 
                         y = fit$coefficients[1, ] + 1.96 * sqrt(mu.variance), 
                         method = "natural", n = 500)
    lower.band <- spline(x = as.numeric(colnames(fit$coefficients)), 
                         y = fit$coefficients[1, ] - 1.96 * sqrt(mu.variance), 
                         method = "natural", n = 500)
    ylim <- range(ylim, upper.band$y, lower.band$y)
  }
  
  plot(x = fit$data$tme, y = fit$data$y, ylim = ylim, xlim = xlim, xaxt="none",
       xlab = '', ylab = '', col = 'black', cex = 1.2, main = '',
       cex.lab = 1.2, font = 2)
  grid()

  for (i in 1:length(fs)) {
    lines(fs[[i]], lty = "dashed", col = 'black', lwd = 0.8)
  }
  
  lines(mu, lwd = 2, col = 'red')
  
  col.meanCurve.rgb <- col2rgb('red')
  if(conf){
    polygon(x = c(upper.band$x, rev(lower.band$x)), 
            y = c(upper.band$y,rev(lower.band$y)), 
            col = rgb(col.meanCurve.rgb[1], 
                      col.meanCurve.rgb[2], col.meanCurve.rgb[3], alpha = 125, 
                      maxColorValue = 255), border = NA)
  }
 
  axis(1, seq(min(fit$data$tme), max(fit$data$tme), length = 7),
       labels = seq(0, 6), font=2)
  mtext(side=1, line=2, "Time (h)", col="black", font=2,cex=1.1)
  mtext(side=2, line=2, "log2(expr)", col="black", font=2,cex=1.1)
  title(main = v , cex.lab = 1.2, line = 0.5)
}



dtwClustCurves <- function(tc.mus, nclust = 6L){
  ## Calculate clusters in parallel
  num.cores <- detectCores(all.tests = FALSE, logical = TRUE)
  
  workers <- makeCluster(num.cores)
  invisible(clusterEvalQ(workers, library("dtwclust")))
  registerDoParallel(workers)
  tc.mus <- lapply(1:ncol(tc.mus), function(i) c(as.numeric(tc.mus[,i])))
  hc_dtw <- tsclust(tc.mus, 
                    type = "h", 
                    k = nclust, 
                    distance = "dtw", 
                    control = hierarchical_control(method = "complete"),
                    centroid = shape_extraction, 
                    #preproc = NULL, 
                    preproc = zscore,
                    trace = T,
                    args = tsclust_args(dist = list(window.size = 4L))
  )
  
  stopCluster(workers)
  registerDoSEQ()
  
  return(hc_dtw)
}

withinCalssReOrder <- function(tc.mus){
  clust.ord <- tc.mus %>% dplyr::select(GeneID, cluster) %>% distinct() %>% 
    group_by(cluster) %>% summarise(GeneSet = list(GeneID)) 
  
  clusters <- unique(tc.mus$cluster)
  hc_eucledian.df <- {}
  for(i in 1:length(clusters)){
    class.i <- tc.mus %>% dplyr::filter(cluster == clusters[i]) %>% dplyr::select(GeneID, y, t) %>%
      pivot_wider(names_from = GeneID, values_from = y) %>% as.data.frame()
    
    
    ## Hierarchical clustering with Eucledian distance
    hc_eucledian <- hclust(dist(t(as.matrix(class.i[,-1] ))), method = "ward.D")
    hc_eucledian.df <- rbind(hc_eucledian.df, 
                             data.frame(GeneID = colnames(class.i[,-1]), 
                                        hc_eucledian.order = hc_eucledian$order,
                                        hc_eucledian.cluster = cutree(hc_eucledian,k = 10)))
  }
  
  return(hc_eucledian.df)
}

## Calculate th overlap of clusters with marker genes and re-order clusters accordingly
matchClustersToPhase <- function(hc_dtw.df, markers.sig){
  join.hc_dtw.df <- inner_join(hc_dtw.df, markers.sig, by = 'GeneID')
  totals <- join.hc_dtw.df %>% group_by(cluster.x) %>% summarise(totals = n())
  overlap <- join.hc_dtw.df %>% group_by(cluster.x, cluster.y) %>% summarise(overlap = n())
  overlap <- left_join(overlap, totals, by = 'cluster.x') %>% mutate(precent = overlap/totals)
  
  colnames(overlap) <- c('cluster', 'markers', 'counts', 'totals', 'percent')
  overlap <- overlap %>% dplyr::select(c('cluster','markers', 'percent')) %>% 
    pivot_wider(names_from = 'markers', values_from = 'percent')
  overlap[is.na(overlap)] <- 0
  overlap <- overlap %>% pivot_longer(-c('cluster'), names_to = 'markers', values_to = 'percent')
  
  overlap$markers[gsub('.*_', '', overlap$markers) == 0] <- 'Ph0'
  overlap$markers[gsub('.*_', '', overlap$markers) == 3] <- 'Ph1'
  overlap$markers[gsub('.*_', '', overlap$markers)== 2] <- 'Ph3'
  overlap$markers[gsub('.*_', '', overlap$markers) == 1] <- 'Ph2'
  overlap$markers[gsub('.*_', '', overlap$markers)== 4] <- 'Ph4'
  overlap$markers[gsub('.*_', '', overlap$markers)== 5] <- 'Ph5'
  
  
  overlap$cluster <- factor(overlap$cluster, levels = unique(sort(overlap$cluster)))
  overlap$markers <- factor(overlap$markers, levels = unique(sort(overlap$markers)))
  
  return(overlap)
  
}

crossCompareMarkers <- function(marker1, marker2, cond1, cond2){
  marker1 <- marker1 %>% 
    select(GeneID, cluster)
  
  marker2 <- marker2 %>% 
    select(GeneID, cluster)
  
  marker1.lists <- marker1 %>% 
    group_by(cluster) %>% summarise(genes = list(GeneID))
  
  marker2.lists <- marker2 %>% 
    group_by(cluster) %>% summarise(genes = list(GeneID))
  
  
  
  marker1.lists$dummy <- 1
  marker2.lists$dummy <- 1
  
  
  
  # RH vs pb 10x
  marker1.marker2.lists <- full_join(marker1.lists, marker2.lists, by = 'dummy')
  
  marker1.marker2.lists <- marker1.marker2.lists %>% rowwise() %>%
    mutate(common.genes = list(intersect(unlist(genes.x), unlist(genes.y))), 
           num.common = length(intersect(unlist(genes.x), unlist(genes.y))))
  
  cluster.sim.marker1.marker2 <- marker1.marker2.lists  %>% 
    select(cluster.x, cluster.y, common.genes, num.common) %>% unnest(common.genes)
  
  colnames(cluster.sim.marker1.marker2) = c(cond1, cond2, 'GeneID', 'num.comm.genes')
  
  cluster.sim.marker1.marker2 <- cluster.sim.marker1.marker2 %>% 
    select(all_of(cond1), all_of(cond2) ,'num.comm.genes') %>% distinct()
  
  cluster.sim.marker1.marker2 <- cluster.sim.marker1.marker2 %>% 
    pivot_wider(names_from = !!cond2, values_from = 'num.comm.genes')
  
  cluster.sim.marker1.marker2[is.na(cluster.sim.marker1.marker2)] <- 0
  
  cluster.sim.marker1.marker2 <- cluster.sim.marker1.marker2 %>% 
    pivot_longer(-cond1, names_to = paste0(cond2), values_to = 'num.comm.genes')
  
  return(cluster.sim.marker1.marker2)
}


processCount <- function(input.dir, filename, tt, rr, down.sample = T){
  file.counts <- read.csv(paste(input.dir, filename, sep = ''))
  genes <- file.counts$X
  expr <- file.counts[,-1]
  rownames(expr) <- genes
  
  
  S.O <- CreateSeuratObject(counts = expr, min.cells = 10, min.features = 100)
  
  #VlnPlot(S.O, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
  #FeatureScatter(S.O, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  
  cutoffs <- quantile(S.O$nCount_RNA, probs = c(0.01, 0.9))
  print(cutoffs)
  S.O <- subset(S.O, subset = nFeature_RNA > cutoffs[1] & nFeature_RNA < cutoffs[2] )
  
  if(down.sample){
    S.O <- subset(x = S.O, downsample = 2000)
  }
  
  pheno <- data.frame(Sample = names(S.O$orig.ident))
  spp <- paste('BDiv', tt, rr, sep = '')
  pheno$spp <- spp
  pheno$time <- tt
  pheno$reactivate <- rr
  
  pheno$NAME <- paste(pheno$spp, pheno$Sample, sep = '_')
  
  
  L <- list(pheno = pheno, S.O = S.O)
  
  return(L)
}


mergeS.O <- function(L){
  num.objs <- length(L)
  
  phenos <- lapply(L, `[[`, 1)
  S.Os <-  lapply(L, `[[`, 2)
  
  all.samples <- bind_rows(phenos)
  
  S.O.merge <- merge(S.Os[[1]], y = S.Os[2:num.objs ], add.cell.ids = unique(all.samples$spp))
  
  rownames(all.samples) <- all.samples$NAME
  S.O.merge <- AddMetaData(S.O.merge, metadata = all.samples)  
  
  return(S.O.merge)
}


processeMergedS.O <- function(S.O.list, orthologs, data.ind = NA, ref.ind = NA, res = 0.2, SC = FALSE){
  ## Merge the S.O, add metadata, and re-split by spp and update the S.O.list
  if(is.na(data.ind[1])){
    data.ind = 1:length(S.O.list)
  }
  
  S.O.merge <- mergeS.O(S.O.list[data.ind])
  S.O.list <- SplitObject(S.O.merge, split.by = "spp")

  
  if(SC){
    S.O.list <- lapply(X = S.O.list, FUN = SCTransform)
    features <- SelectIntegrationFeatures(object.list = S.O.list, nfeatures = 3000)
    S.O.list <- PrepSCTIntegration(object.list = S.O.list, anchor.features = features)
    
    ## Reference-based
    if(!is.na(ref.ind[1])){
      print('using reference index')
      reference_dataset <- ref.ind
      anchors <- FindIntegrationAnchors(object.list = S.O.list, normalization.method = "SCT", 
                                             anchor.features = features, reference = reference_dataset)
    }else{
      
      anchors <- FindIntegrationAnchors(object.list = S.O.list, normalization.method = "SCT", 
                                        anchor.features = features)
    }
    S.O.integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
  }else{
    S.O.list <- lapply(X = S.O.list, FUN = function(x) {
      ## Extract the count data
      
      ## extract the count data from each as.matrix(S.O.list[[1]][["RNA"]]@data)
      ## Replace genes with Bdiv orthologous when needed
      ## recreate the new Seurat object.
      x <- NormalizeData(x)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
    })
    features <- SelectIntegrationFeatures(object.list = S.O.list)

    ## Reference-based
    if(!is.na(ref.ind[1])){
      print('using reference index')
      reference_dataset <- ref.ind
      anchors <- FindIntegrationAnchors(object.list = S.O.list, 
                                        anchor.features = features, reference = reference_dataset)
    }else{
      anchors <- FindIntegrationAnchors(object.list = S.O.list, anchor.features = features)
    }
    S.O.integrated <- IntegrateData(anchorset = anchors)
  }
  
  
  # switch to integrated assay. Make sure to set to RNA for Differential Expression
  DefaultAssay(S.O.integrated) <- "integrated"
  
  # Run the standard workflow for visualization and clustering
  S.O.integrated <- ScaleData(S.O.integrated, verbose = FALSE)
  S.O.integrated <- RunPCA(S.O.integrated, npcs = 30, verbose = FALSE)
  S.O.integrated <- RunUMAP(S.O.integrated, reduction = "pca", dims = 1:30)
  S.O.integrated <- FindNeighbors(S.O.integrated, reduction = "pca", dims = 1:30)
  S.O.integrated <- FindClusters(S.O.integrated, resolution = res)
  
  S.O.integrated$phase.cond <- paste(S.O.integrated@meta.data$spp, 
                                     Idents(S.O.integrated), sep = "_")
  Idents(S.O.integrated) <- "phase.cond"
  
  return(S.O.integrated)
}


getCellCyclePhaseMarkers <- function(all.spp.list){
  all.markers.list <- mclapply(all.spp.list, function(x) FindAllMarkers(object = x, only.pos = TRUE))
  all.markers.list <- lapply(all.markers.list, function(x) {
    x$GeneID = gsub('-', '_', x$gene)
    x$glob.clust <- gsub('.*_', '', x$cluster)
    return(x)
  })
  
  
  all.markers.list.sig <- lapply(all.markers.list, function(x) {
    sig.marker = x %>% dplyr::filter(avg_log2FC > log2(1.5) & p_val_adj < 0.05)
    return(sig.marker)
  })
  
  L <- list(all.markers.list = all.markers.list, 
            all.markers.list.sig = all.markers.list.sig)
  return(L)
}

crossCompareMarkers <- function(marker1, marker2, cond1, cond2){
  marker1 <- marker1 %>% 
    select(GeneID, glob.clust)
  
  marker2 <- marker2 %>% 
    select(GeneID, glob.clust)
  
  marker1.lists <- marker1 %>% 
    group_by(glob.clust) %>% summarise(genes = list(GeneID))
  
  marker2.lists <- marker2 %>% 
    group_by(glob.clust) %>% summarise(genes = list(GeneID))
  
  
  
  marker1.lists$dummy <- 1
  marker2.lists$dummy <- 1
  
  
  marker1.marker2.lists <- full_join(marker1.lists, marker2.lists, by = 'dummy')
  
  marker1.marker2.lists <- marker1.marker2.lists %>% rowwise() %>%
    mutate(common.genes = list(intersect(unlist(genes.x), unlist(genes.y))), 
           num.common = length(intersect(unlist(genes.x), unlist(genes.y))))
  
  cluster.sim.marker1.marker2 <- marker1.marker2.lists  %>% 
    select(glob.clust.x, glob.clust.y, common.genes, num.common) %>% unnest(common.genes)
  
  colnames(cluster.sim.marker1.marker2) = c(cond1, cond2, 'GeneID', 'num.comm.genes')
  
  cluster.sim.marker1.marker2 <- cluster.sim.marker1.marker2 %>% 
    select(all_of(cond1), all_of(cond2),'num.comm.genes') %>% distinct()
  
  cluster.sim.marker1.marker2 <- cluster.sim.marker1.marker2 %>% 
    pivot_wider(names_from = !!cond2, values_from = 'num.comm.genes')
  
  cluster.sim.marker1.marker2[is.na(cluster.sim.marker1.marker2)] <- 0
  
  cluster.sim.marker1.marker2 <- cluster.sim.marker1.marker2 %>% 
    pivot_longer(-cond1, names_to = paste0(cond2), values_to = 'num.comm.genes')
  
  return(cluster.sim.marker1.marker2)
}


markerContrasts <- function(marker1, marker2, cond1, cond2){
  marker1 <- marker1 %>% 
    select(GeneID, glob.clust)
  
  marker2 <- marker2 %>% 
    select(GeneID, glob.clust)
  
  marker1.lists <- marker1 %>% 
    group_by(glob.clust) %>% summarise(genes = list(GeneID))
  
  marker2.lists <- marker2 %>% 
    group_by(glob.clust) %>% summarise(genes = list(GeneID))
  
  
  marker1.lists$dummy <- 1
  marker2.lists$dummy <- 1
  
  marker1.marker2.lists <- full_join(marker1.lists, marker2.lists, by = 'dummy')
  marker1.marker2.lists <- marker1.marker2.lists %>% dplyr::filter(glob.clust.x == glob.clust.y)
  
  marker1.marker2.lists <- marker1.marker2.lists %>% rowwise() %>%
    mutate(common.genes = list(intersect(unlist(genes.x), unlist(genes.y))), 
           num.common = length(intersect(unlist(genes.x), unlist(genes.y))),
           cond1.only.genes = list(setdiff(unlist(genes.x), unlist(genes.y))),
           num.cond1.only.genes = length(setdiff(unlist(genes.x), unlist(genes.y))),
           cond2.only.genes = list(setdiff(unlist(genes.y), unlist(genes.x))),
           num.cond2.only.genes = length(setdiff(unlist(genes.y), unlist(genes.x))))
  
  
  marker1.marker2.lists <- marker1.marker2.lists %>% select(-dummy)
  colnames(marker1.marker2.lists) <- c(paste(cond1, 'cluster', sep = '.'), 
                                       paste(cond1, 'genes', sep = '.'),
                                       paste(cond2, 'cluster', sep = '.'), 
                                       paste(cond2, 'genes', sep = '.'),
                                       'common.genes',
                                       'num.common.genes',
                                       paste(cond1, 'only.genes', sep = '.'),
                                       paste('num', cond1, 'only.genes', sep = '.'),
                                       paste(cond2, 'only.genes', sep = '.'),
                                       paste('num', cond2, 'only.genes', sep = '.'))
  
  
  return(marker1.marker2.lists)
}

makeGlobalHeatMap <- function(global.mat, prod.desc, cut.level = 4){
  mat <- global.mat
  row_dend = hclust(dist(mat))
  base_mean = rowMeans(mat[,2:ncol(mat)])
  mat_scaled = t(apply(mat, 1, scale))
  Sample = factor(colnames(mat) , levels = c("BDiv0hrN", "BDiv4hrN", "BDiv12hrN", "BDiv36hrN", "BDiv7dN"))
  ha = HeatmapAnnotation(
    Sample = Sample, 
    col = list(
      Sample = c("BDiv0hrN" = 'red', "BDiv4hrN" = 'purple', "BDiv12hrN" = 'green', "BDiv36hrN" = 'gold', "BDiv7dN" = 'blue')
      
    ),
    annotation_name_side = "left"
  )
  ht_list = Heatmap(mat_scaled, name = "expression", cluster_columns = F, 
                    col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                    top_annotation = ha, 
                    row_split = cut.level,
                    show_column_names = FALSE, row_title = NULL, show_row_names = F, show_row_dend = FALSE) 
  
  tmp <- cutree(row_dend, cut.level)
  heat.clust <- data.frame(GeneID = gsub('-', '_', names(tmp)), cluster = c(tmp))
  heat.clust <- left_join(heat.clust, prod.desc, by = 'GeneID')
  
  L <- list(heat.clust = heat.clust, ht_list = ht_list)
  return(L)
  
}

makeMatchedHeatMap <- function(matched.mat, prod.desc, cut.level = 4){
  mat <- matched.mat
  row_dend = hclust(dist(mat))
  base_mean = rowMeans(mat[,2:ncol(mat)])
  mat_scaled = t(apply(mat, 1, scale))
  
  
  sample.ord <- as.numeric(gsub('hrN', '', gsub('BDiv', '', gsub('_.*', '', gsub('7d', '168hr', colnames(mat))))))
  clust.ord <- as.numeric(gsub('.*_', '', colnames(mat)))
  tmp <- data.frame(sample.ord = sample.ord, clust.ord = clust.ord)
  col.orders <- with(tmp, order(clust.ord, sample.ord))
  #re.ord <- rep(c(0:(length(unique(sample.ord)) - 1)), rle(sample.ord[order(sample.ord)])$lengths)
  #re.ord <- re.ord[order(sample.ord)]
  #col.orders <- clust.ord + re.ord * length(unique(clust.ord)) + 1
  #col.orders <- order(re.ord +  clust.ord * length(unique(re.ord)) + 1)
  
  
  Sample = factor(gsub('_.*', '', colnames(mat)) , levels = c("BDiv0hrN", "BDiv4hrN", "BDiv12hrN", "BDiv36hrN", "BDiv7dN"))
  Cluster = factor(gsub('.*_', '', colnames(mat)), level=unique(sort(gsub('.*_', '', colnames(mat)))))
  
  ha = HeatmapAnnotation(
    Cluster = Cluster,
    Sample = Sample, 
    col = list(
      Cluster = c("0" = 'yellow', "1" = 'cyan', "2" = 'brown', "3" = 'orange'),
      Sample = c("BDiv0hrN" = 'red', "BDiv4hrN" = 'purple', "BDiv12hrN" = 'green', "BDiv36hrN" = 'gold', "BDiv7dN" = 'blue')
    ),
    annotation_name_side = "left"
  )
  
  ht_list = Heatmap(mat_scaled, name = "expression", cluster_columns = F, column_order = col.orders,
                    col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                    top_annotation = ha, 
                    row_split = cut.level,
                    show_column_names = FALSE, row_title = NULL, show_row_names = F, show_row_dend = FALSE) 
  
  tmp <- cutree(row_dend, cut.level)
  heat.clust <- data.frame(GeneID = gsub('-', '_', names(tmp)), cluster = c(tmp))
  heat.clust <- left_join(heat.clust, prod.desc, by = 'GeneID')
  
  L <- list(heat.clust = heat.clust, ht_list = ht_list)
  return(L)
}

plot_rna_atac_trends.ord <- function(sc.rna.sc.atac.joint.long.sub){
  p  <- ggplot(sc.rna.sc.atac.joint.long.sub, aes(x= time,y=normExpr)) +
    geom_path(aes(color = GeneID),alpha = 0.8, size = 0.8)+ 
    theme_bw() +
    #theme_bw(base_size = 16) +
    theme(legend.position = "right") +
    ylab('normExpr') + xlab('Time') +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 22, face="bold", colour = "black")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 18, face="bold", colour = "black")) +
    theme(strip.background = element_rect(colour="black", fill="white",size=0.5, linetype="solid")) +
    theme(strip.text = element_text(size = 22, face="bold", angle = 0)) + 
    
    coord_cartesian(xlim = c(0,6.5)) + 
    facet_grid(cluster.RNA.ordered ~ data, scales = 'free', space = 'free') +
    theme(
      plot.title = element_text(size=20, face = "bold.italic", color = 'red'),
      axis.title.x = element_text(size=22, face="bold", hjust = 1),
      axis.title.y = element_text(size=22, face="bold")
    ) + 
    theme(#legend.position = c(0.15, 0.85),
      legend.position = 'none',
      legend.title = element_text(colour="black", size=12, 
                                  face="bold"),
      legend.text = element_text(colour="black", size=12, 
                                 face="bold"))
  
  
  return(p)
  
}


plot_atac_sub_clust_ordered <- function(tab){
  p  <- ggplot(tab, aes(x= time,y=normExpr)) +
    geom_path(aes(color = gene_name),alpha = 0.8, size = 0.8)+ 
    theme_bw() +
    theme(legend.position = "right") +
    ylab('normExpr') + xlab('Time') +
    # theme(axis.text.x = element_text(angle = 0, size = 20, face="bold", colour = "black")) +
    # theme(axis.text.y = element_text(angle = 0, size = 20, face="bold", colour = "black")) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(), 
          axis.title = element_blank()) + 
    theme(strip.background = element_rect(colour="black", fill="white",
                                          size=0.5, linetype="solid")) +
    theme(strip.text = element_text(size = 20, face="bold", angle = 0)) + 
    
    coord_cartesian(xlim = c(0,6.5)) + 
    facet_grid(atac.clust.ord ~ ., scales = 'free', space = 'free') +
    theme(
      plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
      # axis.title.x = element_text(size=20, face="bold", hjust = 1),
      # axis.title.y = element_text(size=20, face="bold")
    ) + 
    theme(#legend.position = c(0.15, 0.85),
      legend.position = 'none',
      legend.title = element_text(colour="black", size=12, 
                                  face="bold"),
      legend.text = element_text(colour="black", size=12, 
                                 face="bold"))
  return(p)
  
} 

