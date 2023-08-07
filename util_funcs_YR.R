
## this function plots rna and atac profile of genes in a table of interest
plot_rna_atac <- function(sc.rna.sc.atac.joint.long.sub){
  p  <- ggplot(sc.rna.sc.atac.joint.long.sub, aes(x= time,y=normExpr)) +
    geom_path(aes(color = GeneID),alpha = 0.8, size = 0.8)+ 
    theme_bw() +
    theme(legend.position = "right") +
    ylab('normExpr') + xlab('Time') +
    theme(axis.text.x = element_text(angle = 0, size = 20, face="bold", colour = "black")) +
    theme(axis.text.y = element_text(angle = 0, size = 20, face="bold", colour = "black")) +
    theme(strip.background = element_rect(colour="black", fill="white",
                                          size=0.5, linetype="solid")) +
    theme(strip.text = element_text(size = 20, face="bold", angle = 0)) + 
    
    coord_cartesian(xlim = c(0,6.5)) + 
    facet_grid(. ~ data, scales = 'free', space = 'free') +
    theme(
      plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
      axis.title.x = element_text(size=20, face="bold", hjust = 1),
      axis.title.y = element_text(size=20, face="bold")
    ) + 
    theme(#legend.position = c(0.15, 0.85),
      legend.position = 'none',
      legend.title = element_text(colour="black", size=12, 
                                  face="bold"),
      legend.text = element_text(colour="black", size=12, 
                                 face="bold"))
  
  return(p)
  
}

## this function extracts the expression and accessibility profile of genes of interest
## need to give rna and atac splines and genes of interest as input
## set the scale T/F
get_rna_atac_profile <- function(rna.splines, atac.splines, genes.tab, scale = T) {
  
  
  sc.rna.dtw.wide <- rna.splines %>% 
    pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
    mutate_at(vars(matches("TGME")), ~scale(., center = T, scale = scale)) %>%
    as.data.frame()
  
  sc.atac.dtw.wide <- atac.splines %>% 
    pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
    mutate_at(vars(matches("TGME")), ~scale(., center = T, scale = scale)) %>%
    as.data.frame()
  
  
  sc.rna.mu.scale <- sc.rna.dtw.wide %>% 
    pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')
  
  sc.atac.mu.scale <- sc.atac.dtw.wide %>% 
    pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')
  
  
  ## table of genes to plot their expression 
  tab.genes <- data.frame(TGME49 = gsub("_", "-", genes.tab$gene_name), 
                          Name = genes.tab$ProductDescription)
  
  
  sc.rna.long <- inner_join(sc.rna.mu.scale, tab.genes, by = c('GeneID' = 'TGME49')) 
  sc.rna.long <- sc.rna.long %>% 
    transmute(time = x, GeneID = GeneID, normExpr = expr, Name = Name)  %>% distinct()
  
  sc.atac.long <- inner_join(sc.atac.mu.scale, tab.genes, by = c('GeneID' = 'TGME49')) 
  sc.atac.long <- sc.atac.long %>% 
    transmute(time = x, GeneID = GeneID, normExpr = expr, Name = Name) %>% distinct()
  
  
  sc.rna.sc.atac.joint <- inner_join(sc.rna.long, sc.atac.long, 
                                     by = c("time", "GeneID", "Name"))
  colnames(sc.rna.sc.atac.joint) <- c("time", "GeneID", "scRNA", "Name", "scATAC")
  
  sc.rna.sc.atac.joint.long <- sc.rna.sc.atac.joint %>% 
    pivot_longer(-c('time', "GeneID", "Name"), 
                 names_to = 'data', values_to = 'normExpr') 
  
  return(sc.rna.sc.atac.joint.long)
  
}

## This function performs dynamic time warping clustering for both rna and atac profiles
## it gets a table of genes, 
## the table should include the product description (gene ID and gene Name) 
## you can specify the number of clusters you are looking for (cannot be less than 2)

clust.df <- function(tab , num.clust) {
  
  k <- num.clust
  sc.rna <- sc.rna.dtw.wide[,colnames(sc.rna.dtw.wide) %in% tab$gene_name ]
  sc.atac <- sc.atac.dtw.wide[,colnames(sc.atac.dtw.wide) %in% tab$gene_name]
  sc.atac<- sc.atac.dtw.wide[,colnames(sc.atac.dtw.wide) %in% colnames(sc.rna)]
  
  sc.rna.markers.hc_dtw <- dtwClustCurves(sc.rna, nclust = k)
  sc.atac.markers.hc_dtw <- dtwClustCurves(sc.atac, nclust = k)
  
  tab <- tab[tab$gene_name %in% colnames(sc.rna),]
  
  sc.rna.long <- inner_join(sc.rna.mu.scale, tab, by = c('GeneID' = 'gene_name'))
  sc.rna.long <- sc.rna.long %>%
    transmute(time = x, GeneID = GeneID, normExpr = expr, Name = ProductDescription)  %>% distinct()
  
  sc.atac.long <- inner_join(sc.atac.mu.scale, tab, by = c('GeneID' = 'gene_name'))
  sc.atac.long <- sc.atac.long %>%
    transmute(time = x, GeneID = GeneID, normExpr = expr, Name = ProductDescription) %>% distinct()
  
  
  sc.rna.clust.info <- data.frame(GeneID = colnames(sc.rna), cluster = cutree(sc.rna.markers.hc_dtw, k = k))
  sc.atac.clust.info <- data.frame(GeneID = colnames(sc.atac), cluster = cutree(sc.atac.markers.hc_dtw, k = k))
  
  sc.rna.long.clust <- inner_join(sc.rna.long, sc.rna.clust.info, by = 'GeneID')
  sc.atac.long.clust <- inner_join(sc.atac.long, sc.atac.clust.info, by = 'GeneID')
  
  sc.rna.sc.atac.joint <- inner_join(sc.rna.long.clust, sc.atac.long.clust,
                                     by = c("time", "GeneID", "Name"))
  colnames(sc.rna.sc.atac.joint) <- c("time", "GeneID", "scRNA", "Name", "cluster.RNA", "scATAC", "cluster.ATAC")
  
  sc.rna.sc.atac.joint.long <- sc.rna.sc.atac.joint %>%
    pivot_longer(-c('time', "GeneID", "Name", "cluster.RNA", "cluster.ATAC"),
                 names_to = 'data', values_to = 'normExpr')
  
  sc.rna.sc.atac.joint.long$cluster.RNA <- paste('C', sc.rna.sc.atac.joint.long$cluster.RNA)
  sc.rna.sc.atac.joint.long$cluster.ATAC <- paste('C', sc.rna.sc.atac.joint.long$cluster.ATAC)
  
  return(sc.rna.sc.atac.joint.long)
  
}


clust.atac.df <- function(tab, num.clust = num.clust){
  
  tab <- tab
  k <- num.clust
  
  sc.atac <- sc.atac.dtw.wide[,colnames(sc.atac.dtw.wide) %in% tab$gene_name]
  sc.atac.markers.hc_dtw <- dtwClustCurves(sc.atac, nclust = k)
  tab <- tab[tab$gene_name %in% colnames(sc.atac),]
  sc.atac.long <- inner_join(sc.atac.mu.scale, tab, by = c('GeneID' = 'gene_name'))
  sc.atac.long <- sc.atac.long %>%
    transmute(time = x, GeneID = GeneID, normExpr = expr, group = group, trans.cluster.rna) %>% 
    distinct()
  sc.atac.clust.info <- data.frame(GeneID = colnames(sc.atac), cluster = cutree(sc.atac.markers.hc_dtw, k = k))
  sc.atac.long.clust <- inner_join(sc.atac.long, sc.atac.clust.info, by = 'GeneID')
  colnames(sc.atac.long.clust) <- c("time", "GeneID", "normExpr", "group", "trans.cluster.rna", "cluster.ATAC")
  sc.atac.long.clust$cluster.ATAC <- paste('C', sc.atac.long.clust$cluster.ATAC)
  
  return(sc.atac.long.clust)
}


plot_atac_trand <- function(sc.atac.long.clust){
  
  p  <- ggplot(sc.atac.long.clust, aes(x= time,y=normExpr)) +
    geom_path(aes(color = GeneID),alpha = 0.8, size = 0.8)+ 
    theme_bw(base_size = 14) +
    theme(legend.position = "right") +
    ylab('normExpr') + xlab('Time') +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 15, face="bold", colour = "black")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 15, face="bold", colour = "black")) +
    theme(strip.background = element_rect(colour="black", fill="white",size=0.5, linetype="solid")) +
    theme(strip.text = element_text(size = 18, face="bold", angle = 0)) + 
    
    coord_cartesian(xlim = c(0,6.5)) + 
    facet_grid(cluster.ATAC ~ ., scales = 'free', space = 'free') +
    theme(
      plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
      axis.title.x = element_text(size=18, face="bold", hjust = 1),
      axis.title.y = element_text(size=18, face="bold")
    ) + 
    theme(#legend.position = c(0.15, 0.85),
      legend.position = 'none',
      legend.title = element_text(colour="black", size=12, 
                                  face="bold"),
      legend.text = element_text(colour="black", size=12, 
                                 face="bold")) +
    theme(axis.ticks =  element_blank())
  
  return(p)
}


clust.atac.df <- function(tab, num.clust = num.clust){
  
  tab <- tab
  k <- num.clust
  
  sc.atac <- sc.atac.dtw.wide[,colnames(sc.atac.dtw.wide) %in% tab$gene_name]
  sc.atac.markers.hc_dtw <- dtwClustCurves(sc.atac, nclust = k)
  tab <- tab[tab$gene_name %in% colnames(sc.atac),]
  sc.atac.long <- inner_join(sc.atac.mu.scale, tab, by = c('GeneID' = 'gene_name'))
  sc.atac.long <- sc.atac.long %>%
    transmute(time = x, GeneID = GeneID, normExpr = expr, group = group, trans.cluster.rna) %>% 
    distinct()
  sc.atac.clust.info <- data.frame(GeneID = colnames(sc.atac), cluster = cutree(sc.atac.markers.hc_dtw, k = k))
  sc.atac.long.clust <- inner_join(sc.atac.long, sc.atac.clust.info, by = 'GeneID')
  colnames(sc.atac.long.clust) <- c("time", "GeneID", "normExpr", "group", "trans.cluster.rna", "cluster.ATAC")
  sc.atac.long.clust$cluster.ATAC <- paste('C', sc.atac.long.clust$cluster.ATAC)
  
  return(sc.atac.long.clust)
}


plot_atac_trand <- function(sc.atac.long.clust){
  
  p  <- ggplot(sc.atac.long.clust, aes(x= time,y=normExpr)) +
    geom_path(aes(color = GeneID),alpha = 0.8, size = 0.8)+ 
    theme_bw(base_size = 14) +
    theme(legend.position = "right") +
    ylab('normExpr') + xlab('Time') +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 15, face="bold", colour = "black")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 15, face="bold", colour = "black")) +
    theme(strip.background = element_rect(colour="black", fill="white",size=0.5, linetype="solid")) +
    theme(strip.text = element_text(size = 18, face="bold", angle = 0)) + 
    
    coord_cartesian(xlim = c(0,6.5)) + 
    facet_grid(cluster.ATAC ~ ., scales = 'free', space = 'free') +
    theme(
      plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
      axis.title.x = element_text(size=18, face="bold", hjust = 1),
      axis.title.y = element_text(size=18, face="bold")
    ) + 
    theme(#legend.position = c(0.15, 0.85),
      legend.position = 'none',
      legend.title = element_text(colour="black", size=12, 
                                  face="bold"),
      legend.text = element_text(colour="black", size=12, 
                                 face="bold")) +
    theme(axis.ticks =  element_blank())
  
  return(p)
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


getPcaMetaData <- function(S.O){
  pc <- S.O@reductions$pca@cell.embeddings
  pc <- data.frame(pc) %>% dplyr::mutate(Sample = rownames(pc)) %>% 
    transmute(Sample = Sample, PC_1 = PC_1, PC_2 = PC_2, PC_3 = PC_3)
  umap <- S.O@reductions$umap@cell.embeddings
  umap <- data.frame(umap) %>% dplyr::mutate(Sample = rownames(umap)) %>% 
    transmute(Sample = Sample, UMAP_1 = UMAP_1, UMAP_2 = UMAP_2)
  
  meta.data <- data.frame(Sample = rownames(S.O@meta.data), 
                          spp = S.O@meta.data$spp2, phase = S.O@meta.data$phase)
  meta.data <- left_join(meta.data,
                         pc, by = 'Sample')
  meta.data <- left_join(meta.data, umap, by = 'Sample')
  return(meta.data)  
}


getPcaMetaData.trans <- function(S.O){
  pc <- S.O@reductions$pca@cell.embeddings
  pc <- data.frame(pc) %>% dplyr::mutate(Sample = rownames(pc)) %>% 
    transmute(Sample = Sample, PC_1 = PC_1, PC_2 = PC_2, PC_3 = PC_3)
  umap <- S.O@reductions$umap@cell.embeddings
  umap <- data.frame(umap) %>% dplyr::mutate(Sample = rownames(umap)) %>% 
    transmute(Sample = Sample, UMAP_1 = UMAP_1, UMAP_2 = UMAP_2)
  
  meta.data <- data.frame(Sample = rownames(S.O@meta.data), 
                          spp = S.O@meta.data$spp2, phase = S.O@meta.data$phase,
                          transition.rna = S.O@meta.data$transition.rna, 
                          transition.atac = S.O@meta.data$transition.atac)
  meta.data <- left_join(meta.data,
                         pc, by = 'Sample')
  meta.data <- left_join(meta.data, umap, by = 'Sample')
  return(meta.data)  
}

getPcaMetaData.atac.trans <- function(S.O){
  
  pc <- S.O@reductions$pca@cell.embeddings
  pc <- data.frame(pc) %>% dplyr::mutate(Sample = rownames(pc)) %>% 
    transmute(Sample = Sample, PC_1 = PC_1, PC_2 = PC_2, PC_3 = PC_3)
  umap <- S.O@reductions$umap@cell.embeddings
  umap <- data.frame(umap) %>% dplyr::mutate(Sample = rownames(umap)) %>% 
    transmute(Sample = Sample, UMAP_1 = UMAP_1, UMAP_2 = UMAP_2)
  
  meta.data <- data.frame(Sample = rownames(S.O@meta.data), 
                          spp = S.O@meta.data$spp2, 
                          phase = S.O@meta.data$phase,
                          #transition.rna = S.O@meta.data$transition.rna, 
                          transition.atac = S.O@meta.data$transition.atac)
  meta.data <- left_join(meta.data,
                         pc, by = 'Sample')
  meta.data <- left_join(meta.data, umap, by = 'Sample')
  return(meta.data) 
}
num.cores <- detectCores(all.tests = FALSE, logical = TRUE)


picewise_scale_V2 <- function(sds.data){
  
  
  # Scale the pt based on know biology: Radke et. al 2000
  G <- c(0, 3) # 3h
  S <- c(3, 4.7) # 1.7h
  M <- c(4.7, 5) # ~20 min
  C <- c(5, 6) # 1h
  
  ind.G1.a <- which(sds.data$phase == 'G1.a')
  ind.G1.b <- which(sds.data$phase == 'G1.b')
  ind.S <- which(sds.data$phase == 'S')
  ind.M <- which(sds.data$phase == 'M')
  ind.C <- which(sds.data$phase == 'C')
  
  
  #t1 <- quantile(t[ind.G1.b], prob=0.25) ## Start of G1.b
  t1 <- (quantile(sds.data$pt[ind.G1.a], prob=0.75) +
           quantile(sds.data$pt[ind.G1.b], prob=0.25))/2 ## Start of G1.b
  t2 <- (quantile(sds.data$pt[ind.G1.b], prob=0.75) + 
           quantile(sds.data$pt[ind.S], prob=0.25)) / 2 ## End of G1.b, Start of S
  t3 <- (quantile(sds.data$pt[ind.S], prob=0.75) + 
           quantile(sds.data$pt[ind.M], prob=0.25)) / 2## End of S, Start of M
  t4 <- (quantile(sds.data$pt[ind.M], prob=0.75) + 
           quantile(sds.data$pt[ind.C], prob=0.25)) / 2## End of M, Start of C
  
  t0 <- 0
  t5 <- 6
  
  # sds.data$pt.shift <- (sds.data$pt + t.shift) %% 6
  # sds.data$pt.shift <- 6 * ((sds.data$pt.shift - min(sds.data$pt.shift)) / 
  #                             (max(sds.data$pt.shift) - min(sds.data$pt.shift)))
  # plot(sds.data$phase, sds.data$pt.shift)
  
  slp.g <- (G[2] - G[1]) / (t2 - t0)
  inc.g  <- c(t0, G[1])
  
  slp.s <- (S[2] - S[1]) / (t3 - t2)
  inc.s  <- c(t2, S[1])
  
  slp.m <- (M[2] - M[1]) / (t4 - t3)
  inc.m  <- c(t3, M[1])
  
  slp.c <- (C[2] - C[1]) / (t5 - t4)
  inc.c  <- c(t4, C[1])
  
  s.t <- data.frame(pt = seq(0, 6, 0.01))  
  s.t <- s.t %>% mutate(phase = case_when(pt >= t0 & pt < t1 ~ 'G1.a',
                                          pt >= t1 & pt < t2 ~ 'G1.b',
                                          pt >= t2 & pt < t3 ~ 'S',
                                          pt >= t3 & pt < t4 ~ 'M',
                                          pt >= t4 ~ 'C'), 
                        t = case_when(pt >= t0 & pt < t2 ~ inc.g[2] + slp.g * (pt - inc.g[1]),
                                      pt >= t2 & pt < t3 ~ inc.s[2] + slp.s * (pt - inc.s[1]),
                                      pt >= t3 & pt < t4 ~ inc.m[2] + slp.m * (pt - inc.m[1]),
                                      pt >= t4 ~ inc.c[2] + slp.c * (pt - inc.c[1])))
  
}


clust.df <- function(tab , num.clust) {
  
  k <- num.clust
  sc.rna <- sc.rna.dtw.wide[,colnames(sc.rna.dtw.wide) %in% tab$gene_name ]
  sc.atac <- sc.atac.dtw.wide[,colnames(sc.atac.dtw.wide) %in% tab$gene_name]
  sc.atac<- sc.atac.dtw.wide[,colnames(sc.atac.dtw.wide) %in% colnames(sc.rna)]
  
  sc.rna.markers.hc_dtw <- dtwClustCurves(sc.rna, nclust = k)
  sc.atac.markers.hc_dtw <- dtwClustCurves(sc.atac, nclust = k)
  
  tab <- tab[tab$gene_name %in% colnames(sc.rna),]
  
  sc.rna.long <- inner_join(sc.rna.mu.scale, tab, by = c('GeneID' = 'gene_name'))
  sc.rna.long <- sc.rna.long %>%
    transmute(time = x, GeneID = GeneID, normExpr = expr, Name = ProductDescription)  %>% distinct()
  
  sc.atac.long <- inner_join(sc.atac.mu.scale, tab, by = c('GeneID' = 'gene_name'))
  sc.atac.long <- sc.atac.long %>%
    transmute(time = x, GeneID = GeneID, normExpr = expr, Name = ProductDescription) %>% distinct()
  
  
  sc.rna.clust.info <- data.frame(GeneID = colnames(sc.rna), cluster = cutree(sc.rna.markers.hc_dtw, k = k))
  sc.atac.clust.info <- data.frame(GeneID = colnames(sc.atac), cluster = cutree(sc.atac.markers.hc_dtw, k = k))
  
  sc.rna.long.clust <- inner_join(sc.rna.long, sc.rna.clust.info, by = 'GeneID')
  sc.atac.long.clust <- inner_join(sc.atac.long, sc.atac.clust.info, by = 'GeneID')
  
  sc.rna.sc.atac.joint <- inner_join(sc.rna.long.clust, sc.atac.long.clust,
                                     by = c("time", "GeneID", "Name"))
  colnames(sc.rna.sc.atac.joint) <- c("time", "GeneID", "scRNA", "Name", "cluster.RNA", "scATAC", "cluster.ATAC")
  
  sc.rna.sc.atac.joint.long <- sc.rna.sc.atac.joint %>%
    pivot_longer(-c('time', "GeneID", "Name", "cluster.RNA", "cluster.ATAC"),
                 names_to = 'data', values_to = 'normExpr')
  
  sc.rna.sc.atac.joint.long$cluster.RNA <- paste('C', sc.rna.sc.atac.joint.long$cluster.RNA)
  sc.rna.sc.atac.joint.long$cluster.ATAC <- paste('C', sc.rna.sc.atac.joint.long$cluster.ATAC)
  
  return(sc.rna.sc.atac.joint.long)
  
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

# 

## This function performs dynamic time warping clustering for both rna and atac profiles
## it gets a table of genes, 
## the table should include the product description (gene ID and gene Name) 
## you can specify the number of clusters you are looking for (cannot be less than 2)

clust.df <- function(tab , num.clust) {
  
  k <- num.clust
  sc.rna <- sc.rna.dtw.wide[,colnames(sc.rna.dtw.wide) %in% tab$gene_name ]
  sc.atac <- sc.atac.dtw.wide[,colnames(sc.atac.dtw.wide) %in% tab$gene_name]
  sc.atac<- sc.atac.dtw.wide[,colnames(sc.atac.dtw.wide) %in% colnames(sc.rna)]
  
  sc.rna.markers.hc_dtw <- dtwClustCurves(sc.rna, nclust = k)
  sc.atac.markers.hc_dtw <- dtwClustCurves(sc.atac, nclust = k)
  
  tab <- tab[tab$gene_name %in% colnames(sc.rna),]
  
  sc.rna.long <- inner_join(sc.rna.mu.scale, tab, by = c('GeneID' = 'gene_name'))
  sc.rna.long <- sc.rna.long %>%
    transmute(time = x, GeneID = GeneID, normExpr = expr, Name = ProductDescription)  %>% distinct()
  
  sc.atac.long <- inner_join(sc.atac.mu.scale, tab, by = c('GeneID' = 'gene_name'))
  sc.atac.long <- sc.atac.long %>%
    transmute(time = x, GeneID = GeneID, normExpr = expr, Name = ProductDescription) %>% distinct()
  
  
  sc.rna.clust.info <- data.frame(GeneID = colnames(sc.rna), cluster = cutree(sc.rna.markers.hc_dtw, k = k))
  sc.atac.clust.info <- data.frame(GeneID = colnames(sc.atac), cluster = cutree(sc.atac.markers.hc_dtw, k = k))
  
  sc.rna.long.clust <- inner_join(sc.rna.long, sc.rna.clust.info, by = 'GeneID')
  sc.atac.long.clust <- inner_join(sc.atac.long, sc.atac.clust.info, by = 'GeneID')
  
  sc.rna.sc.atac.joint <- inner_join(sc.rna.long.clust, sc.atac.long.clust,
                                     by = c("time", "GeneID", "Name"))
  colnames(sc.rna.sc.atac.joint) <- c("time", "GeneID", "scRNA", "Name", "cluster.RNA", "scATAC", "cluster.ATAC")
  
  sc.rna.sc.atac.joint.long <- sc.rna.sc.atac.joint %>%
    pivot_longer(-c('time', "GeneID", "Name", "cluster.RNA", "cluster.ATAC"),
                 names_to = 'data', values_to = 'normExpr')
  
  sc.rna.sc.atac.joint.long$cluster.RNA <- paste('C', sc.rna.sc.atac.joint.long$cluster.RNA)
  sc.rna.sc.atac.joint.long$cluster.ATAC <- paste('C', sc.rna.sc.atac.joint.long$cluster.ATAC)
  
  return(sc.rna.sc.atac.joint.long)
  
}


## plot the expression and accessibility of genes within each cluster
# facet by rna cluster
plot_rna_atac_trends <- function(sc.rna.sc.atac.joint.long.sub){
  p  <- ggplot(sc.rna.sc.atac.joint.long.sub, aes(x= time,y=normExpr)) +
    geom_path(aes(color = GeneID),alpha = 0.8, size = 0.8)+ 
    theme_bw(base_size = 14) +
    theme(legend.position = "right") +
    ylab('normExpr') + xlab('Time') +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 18, face="bold", colour = "black")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 18, face="bold", colour = "black")) +
    theme(strip.background = element_rect(colour="black", fill="white",size=0.5, linetype="solid")) +
    theme(strip.text = element_text(size = 18, face="bold", angle = 0)) + 
    
    coord_cartesian(xlim = c(0,6.5)) + 
    facet_grid(cluster.RNA ~ data, scales = 'free', space = 'free') +
    theme(
      plot.title = element_text(size=20, face = "bold.italic", color = 'red'),
      axis.title.x = element_text(size=18, face="bold", hjust = 1),
      axis.title.y = element_text(size=18, face="bold")
    ) + 
    theme(#legend.position = c(0.15, 0.85),
      legend.position = 'none',
      legend.title = element_text(colour="black", size=12, 
                                  face="bold"),
      legend.text = element_text(colour="black", size=12, 
                                 face="bold"))
  
  
  return(p)
  
}

## plot the expression and accessibility of genes within each cluster
#facet by atac cluster
plot_atac_trends <- function(sc.rna.sc.atac.joint.long.sub){
  p  <- ggplot(sc.rna.sc.atac.joint.long.sub, aes(x= time,y=normExpr)) +
    geom_path(aes(color = GeneID),alpha = 0.8, size = 0.8)+ 
    theme_bw(base_size = 14) +
    theme(legend.position = "right") +
    ylab('normExpr') + xlab('Time') +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 18, face="bold", colour = "black")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 18, face="bold", colour = "black")) +
    theme(strip.background = element_rect(colour="black", fill="white",size=0.5, linetype="solid")) +
    theme(strip.text = element_text(size = 18, face="bold", angle = 0)) + 
    
    coord_cartesian(xlim = c(0,6.5)) + 
    facet_grid(cluster.ATAC ~ data, scales = 'free', space = 'free') +
    theme(
      plot.title = element_text(size=20, face = "bold.italic", color = 'red'),
      axis.title.x = element_text(size=18, face="bold", hjust = 1),
      axis.title.y = element_text(size=18, face="bold")
    ) + 
    theme(#legend.position = c(0.15, 0.85),
      legend.position = 'none',
      legend.title = element_text(colour="black", size=12, 
                                  face="bold"),
      legend.text = element_text(colour="black", size=12, 
                                 face="bold"))
  
  
  return(p)
  
}

# expression and atac profile of one gene at a time
plot_trends <- function(my.GeneID, sc.rna.spline.fits,sc.atac.spline.fits ){
  
  
  ## Turn the data into wide format (time by gene) and center & scale each gene
  sc.rna.dtw.wide <- sc.rna.spline.fits %>% 
    pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
    mutate_at(vars(matches("TGME")), ~scale(., center = F, scale = F)) %>%
    as.data.frame()
  
  sc.atac.dtw.wide <- sc.atac.spline.fits %>% 
    pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
    mutate_at(vars(matches("TGME")), ~scale(., center = F, scale = F)) %>%
    as.data.frame()
  
  
  sc.rna.mu.scale <- sc.rna.dtw.wide %>% 
    pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')
  
  sc.atac.mu.scale <- sc.atac.dtw.wide %>% 
    pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')
  
  
  my.rna <- sc.rna.mu.scale %>% dplyr::filter(GeneID == my.GeneID)
  my.atac <- sc.atac.mu.scale %>% dplyr::filter(GeneID == my.GeneID)
  
  p1  <- ggplot(my.rna , aes(x= x,y=expr)) +
    geom_line(color = 'blue',alpha = 0.8, size = 0.8)+ 
    theme_bw(base_size = 14) +
    theme(legend.position = "right") +
    ylab('normExpr') + xlab('Time') +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
    theme(strip.background = element_rect(colour="black", fill="white",
                                          size=0.5, linetype="solid")) +
    ggtitle(paste('rna', my.GeneID)) + 
    theme(strip.text = element_text(size = 14, face="bold", angle = 0)) + 
    theme(
      plot.title = element_text(size=14, face = "bold.italic", color = 'black'),
      axis.title.x = element_text(size=14, face="bold", hjust = 1),
      axis.title.y = element_text(size=14, face="bold")
    ) 
  
  
  p2  <- ggplot(my.atac , aes(x= x,y=expr)) +
    geom_line(color = 'red',alpha = 0.8, size = 0.8)+ 
    theme_bw(base_size = 14) +
    theme(legend.position = "right") +
    ylab('normExpr') + xlab('Time') +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
    theme(strip.background = element_rect(colour="black", fill="white",
                                          size=0.5, linetype="solid")) +
    ggtitle(paste('atac', my.GeneID)) + 
    theme(strip.text = element_text(size = 14, face="bold", angle = 0)) + 
    theme(
      plot.title = element_text(size=14, face = "bold.italic", color = 'black'),
      axis.title.x = element_text(size=14, face="bold", hjust = 1),
      axis.title.y = element_text(size=14, face="bold")
    ) 
  
  p <- grid.arrange(p1, p2)
  
  return(p)
}


## this function plots rna and atac profile of genes in a table of interest

plot_rna_atac <- function(sc.rna.sc.atac.joint.long.sub){
  p  <- ggplot(sc.rna.sc.atac.joint.long.sub, aes(x= time,y=normExpr)) +
    geom_path(aes(color = GeneID),alpha = 0.8, size = 0.8)+ 
    theme_bw() +
    theme(legend.position = "right") +
    ylab('normExpr') + xlab('Time') +
    theme(axis.text.x = element_text(angle = 0, size = 20, face="bold", colour = "black")) +
    theme(axis.text.y = element_text(angle = 0, size = 20, face="bold", colour = "black")) +
    theme(strip.background = element_rect(colour="black", fill="white",
                                          size=0.5, linetype="solid")) +
    theme(strip.text = element_text(size = 20, face="bold", angle = 0)) + 
    
    coord_cartesian(xlim = c(0,6.5)) + 
    facet_grid(. ~ data, scales = 'free', space = 'free') +
    theme(
      plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
      axis.title.x = element_text(size=20, face="bold", hjust = 1),
      axis.title.y = element_text(size=20, face="bold")
    ) + 
    theme(#legend.position = c(0.15, 0.85),
      legend.position = 'none',
      legend.title = element_text(colour="black", size=12, 
                                  face="bold"),
      legend.text = element_text(colour="black", size=12, 
                                 face="bold"))
  
  return(p)
  
}


## this function extracts the expression and accessibility profile of genes of interest
## need to give rna and atac splines and genes of interest as input
## set the scale T/F
get_rna_atac_profile <- function(rna.splines, atac.splines, genes.tab, scale = T) {
  
  
  sc.rna.dtw.wide <- rna.splines %>% 
    pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
    mutate_at(vars(matches("TGME")), ~scale(., center = T, scale = scale)) %>%
    as.data.frame()
  
  sc.atac.dtw.wide <- atac.splines %>% 
    pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
    mutate_at(vars(matches("TGME")), ~scale(., center = T, scale = scale)) %>%
    as.data.frame()
  
  
  sc.rna.mu.scale <- sc.rna.dtw.wide %>% 
    pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')
  
  sc.atac.mu.scale <- sc.atac.dtw.wide %>% 
    pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')
  
  
  ## table of genes to plot their expression 
  tab.genes <- data.frame(TGME49 = gsub("_", "-", genes.tab$gene_name), 
                          Name = genes.tab$ProductDescription)
  
  
  sc.rna.long <- inner_join(sc.rna.mu.scale, tab.genes, by = c('GeneID' = 'TGME49')) 
  sc.rna.long <- sc.rna.long %>% 
    transmute(time = x, GeneID = GeneID, normExpr = expr, Name = Name)  %>% distinct()
  
  sc.atac.long <- inner_join(sc.atac.mu.scale, tab.genes, by = c('GeneID' = 'TGME49')) 
  sc.atac.long <- sc.atac.long %>% 
    transmute(time = x, GeneID = GeneID, normExpr = expr, Name = Name) %>% distinct()
  
  
  sc.rna.sc.atac.joint <- inner_join(sc.rna.long, sc.atac.long, 
                                     by = c("time", "GeneID", "Name"))
  colnames(sc.rna.sc.atac.joint) <- c("time", "GeneID", "scRNA", "Name", "scATAC")
  
  sc.rna.sc.atac.joint.long <- sc.rna.sc.atac.joint %>% 
    pivot_longer(-c('time', "GeneID", "Name"), 
                 names_to = 'data', values_to = 'normExpr') 
  
  return(sc.rna.sc.atac.joint.long)
  
}
