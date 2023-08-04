
source('./util_funcs.R')
source('./loadlb.R')

num.cores <- detectCores(all.tests = FALSE, logical = TRUE)

## NEW - prod desc MJ
prod.desc <- read.xlsx("../Input/Toxo_genomics/genes/MJ_annotation_07_27_2023.xlsx")
atac_sub <- readRDS('../Input/toxo_cdc/rds_ME49_59/S.O_intra_atac_lables_pt.rds')
rna_sub <- readRDS('../Input/toxo_cdc/rds_ME49_59/S.O_intra_lables_pt.rds')

### Differential gene expression
Idents(rna_sub) <- 'phase'
DefaultAssay(rna_sub) <- 'RNA'
Intra.markers <- FindAllMarkers(object = rna_sub, only.pos = T, min.pct = 0)

Intra.markers$GeneID <- gsub('-', '_', Intra.markers$gene)
Intra.markers.top <- Intra.markers %>% group_by(cluster) %>% top_n(2, avg_log2FC)

Intra.markers.sig <- Intra.markers %>% dplyr::filter(avg_log2FC > log2(1.5) & p_val_adj < 0.05)
Intra.markers.sig <- left_join(Intra.markers.sig, prod.desc, by = c("gene" = "TGME49") )

ss <- Intra.markers.sig %>% group_by(cluster) %>% summarise(num.DEG = n())

write.xlsx(Intra.markers.sig, "../OutPut/toxo_cdc/ME49_59/tables/DEGs_Canonical_cell_cycle_phases.xlsx")
saveRDS(Intra.markers.sig, '../Input/toxo_cdc/rds_ME49_59/Intra_markers_sig.rds')


###################################################
## marker analysis using inferred transition points 
###################################################

prod.desc <- read.xlsx("../Input/Toxo_genomics/genes/MJ_annotation_07_27_2023.xlsx")
rna_sub <- readRDS('../Input/toxo_cdc/rds_ME49_59/S.O.intra_rna_atac_trnasition_v2.rds')
atac_sub <- readRDS('../Input/toxo_cdc/rds_ME49_59/S.O.intra_atac_atac_trnasition_v2.rds')

## Differential gene expression usin rna transition
## rna markers using rna transition points


Idents(rna_sub) <- 'transition.rna'
rna.markers.rna.trans <- FindAllMarkers(object = rna_sub, only.pos = T, min.pct = 0.0)
rna.markers.rna.trans$GeneID <- gsub('-', '_', rna.markers.rna.trans$gene)

rna.sig.markers.rna.trans <- rna.markers.rna.trans %>% dplyr::filter(avg_log2FC > log2(1.5) & p_val_adj < 0.05)
dim(rna.sig.markers.rna.trans)


rna.sig.markers.rna.trans <- left_join(rna.sig.markers.rna.trans, prod.desc, by = c("gene" = "TGME49") )
rna.sig.markers.rna.trans <- rna.sig.markers.rna.trans %>% 
  mutate(new.prod.desc = ifelse(ProductDescription == "hypothetical protein" & !is.na(new.name), new.name, ProductDescription))

ss.rna <- rna.sig.markers.rna.trans %>% group_by(cluster) %>% summarise(num.DEG = n())
sum(ss.rna$num.DEG)

saveRDS(rna.sig.markers.rna.trans, '../Input/toxo_cdc/rds_ME49_59/rna_markers_rna_trns_sig_v2.rds')
write.xlsx(rna.sig.markers.rna.trans, '../Output/toxo_cdc/ME49_59/tables/rna_markers_sig_rna_trns_v2.xlsx')



rna.sig.markers.rna.trans <- readRDS('../Input/toxo_cdc/rds_ME49_59/rna_markers_rna_trns_sig_v2.rds')
rna.sig.markers.rna.trans <- rna.sig.markers.rna.trans %>%
  mutate(Category = ifelse(str_detect(new.prod.desc, "hypothetical protein"), "hypo.", "others"))

ss.rna <- rna.sig.markers.rna.trans %>% group_by(cluster, Category) %>% summarise(num.DEG = n())
ss.rna$Color <- c(rep("#ff9a00", 2), rep("#9ca820", 2), rep("#615FB1", 2), rep("#8f139f", 2))

ss.rna <- ss.rna %>%
  mutate( ## for this you will need to remove the whitespace 
    Category = stringr::str_trim(Category),
    newcolor = ifelse(grepl("hypo", Category), alpha(Color, .5), Color)
  ) 
sum(ss.rna$num.DEG)

saveRDS(ss.rna, "../Input/toxo_cdc/rds_ME49_59/rna_markers_rna_trns_sig_v2_sum_plt.rds")


p <- ggplot(ss.rna, aes(cluster, num.DEG)) +
  geom_col(aes(fill = I(newcolor), color = Category),
           position = position_stack(reverse = FALSE),
           ##remove the border 
           linewidth = 0) +
  geom_text(aes(label = num.DEG, group = Category),size = 8, color = "white",
            fontface = 'bold', position = position_stack(vjust = .5, reverse = TRUE)) +
  ## change the legend fill manually 
  guides(color = guide_legend(
    reverse = FALSE,
    override.aes = list(fill = c("grey25", alpha("grey25", .5))))) +
  theme(legend.position = "none") +
  theme_bw()+
  theme(
    axis.text.x = element_blank(),
    #axis.text.y = element_text(size = 10),
    axis.text.y = element_blank(),
    axis.ticks = element_blank())+
  theme(strip.background=element_rect(fill='white', color = 'black'),
        panel.spacing = unit(1.5, "lines"), 
        strip.text.x=element_text(angle=0, hjust=0.5,vjust=0.5, size = 14,face = 'bold'),
        plot.title = element_text(size=16, face = "bold.italic", color = 'black'),
        axis.title.x = element_text(size=22, face="bold"),
        axis.title.y = element_text(size=22, face="bold"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        strip.text = element_text(face = "bold", size = 24,  angle = 0), 
        strip.placement = "outside") +
  theme(legend.text = element_text(face = "bold", size = 16),
        legend.title = element_text(face = "bold", size = 18))+
  theme(panel.spacing = unit(1.5, "lines")) +
  ggtitle("scRNA markers - rna transition")+
  theme(legend.position = "none")
p


ggsave(filename="../Output/toxo_cdc/ME49_59/figures_paper/rna_numbers_rna_transition_v2_new_annotation.pdf",
       plot=p,
       width = 6, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

## table for performing GO term 
rna.sig.markers.rna.trans <- readRDS('../Input/toxo_cdc/rds_ME49_59/rna_markers_rna_trns_sig_v2.rds')

rna.sig.markers.rna.trans.sum <- rna.sig.markers.rna.trans %>%
  group_by(cluster) %>% summarise(genes = list(GeneID), total = n())
write.xlsx(rna.sig.markers.rna.trans.sum, "../Output/toxo_cdc/ME49_59/tables/rna_sig_markers_rna_trans_sum_v2.xlsx")



## Differential gene expression usin atac transition
## rna markers using atac transition points

Idents(rna_sub) <- 'transition.atac'
transition.markers <- FindAllMarkers(object = rna_sub, only.pos = T, min.pct = 0)
transition.markers$GeneID <- gsub('-', '_', transition.markers$gene)
transition.markers.sig <- transition.markers %>% dplyr::filter(avg_log2FC > log2(1.5) & p_val_adj < 0.05)
dim(transition.markers.sig)

saveRDS(transition.markers.sig, '../Input/toxo_cdc/rds_ME49_59/rna_markers_sig_atac_trans_v2.rds')
write.xlsx(transition.markers.sig, '../Output/toxo_cdc/ME49_59/tables/rna_markers_sig_atac_trns_v2.xlsx')

ss.atac <- transition.markers.sig %>% group_by(cluster) %>% summarise(num.DEG = n())
ss.atac$cluster <- factor(ss.atac$cluster, levels = c('T1', 'T2', 'T3', 'T4'))



## ATAC
## atac markers using atac transitions  (this is based on gene activity)
## we do not use this analysis 
Idents(atac_sub) <- 'transition.atac'
DefaultAssay(atac_sub) <- "RNA"

atac.markers.atac.trans <- FindAllMarkers(object = atac_sub, only.pos = T, min.pct = 0.0)
atac.markers.atac.trans$GeneID <- gsub('-', '_', atac.markers.atac.trans$gene)

atac.sig.markers.atac.trans <- atac.markers.atac.trans %>% dplyr::filter(avg_log2FC > log2(1.5) & p_val_adj < 0.05)
dim(atac.sig.markers.atac.trans)


### intersect the rna marker genes (rna and atac transitions were used to perform DEG)
## we do not use this analysis 
rna.sig.markers.rna.trans <- read.xlsx('../Output/toxo_cdc/ME49_59/tables/rna_markers_sig_rna_trns_v2.xlsx')
rna.sig.markers.atac.trans <- read.xlsx('../Output/toxo_cdc/ME49_59/tables/rna_markers_sig_atac_trns_v2.xlsx')

rna.markers.rna.trans.sum <- rna.sig.markers.rna.trans %>% 
  group_by(cluster) %>% summarise( genes = list(GeneID), total = n())

rna.markers.atac.trans.sum <- rna.sig.markers.atac.trans %>%
  group_by(cluster) %>% summarise(genes = list(GeneID), total = n())

XX <- left_join(rna.markers.rna.trans.sum, rna.markers.atac.trans.sum, by = "cluster")
XX <- XX %>% rowwise() %>%
  mutate(overlap = list(intersect(genes.x, genes.y)), 
         overlap.num = length(intersect(genes.x, genes.y)),
         diff = list(setdiff(genes.x, genes.y)), 
         diff.num = length(setdiff(genes.x, genes.y)))


## contingency 
rna.sig.markers.rna.trans <- read.xlsx('../Output/toxo_cdc/ME49_59/tables/rna_markers_sig_rna_trns_v2.xlsx')
rna.sig.markers.rna.trans$data <- "rna.transition"
#rna.sig.markers.rna.trans <- rna.sig.markers.rna.trans %>% dplyr::select(GeneID, cluster, data) 
rna.sig.markers.rna.trans$cluster.trans <- paste(rna.sig.markers.rna.trans$cluster, "rna", sep = ".")

rna.sig.markers.atac.trans <- read.xlsx('../Output/toxo_cdc/ME49_59/tables/rna_markers_sig_atac_trns_v2.xlsx')
rna.sig.markers.atac.trans$data <- "atac.transition"
#rna.sig.markers.atac.trans <- rna.sig.markers.atac.trans %>% dplyr::select(GeneID, cluster , data)
rna.sig.markers.atac.trans$cluster.trans <- paste(rna.sig.markers.atac.trans$cluster, "atac", sep = ".")

DD <- left_join(rna.sig.markers.rna.trans, rna.sig.markers.atac.trans, by = "GeneID" ,  multiple = "all") %>% na.omit()
DD.tab <- table(DD$cluster.trans.x, DD$cluster.trans.y)

DD.overlap <- DD %>% dplyr::filter(cluster.x == cluster.y)
write.xlsx(DD.overlap, "../Output/toxo_cdc/ME49_59/tables/rna_markers_sig_rna_atac_trans_ovlp_lng.xlsx")
write.xlsx(DD.tab,"../Output/toxo_cdc/ME49_59/tables/rna_markers_sig_rna_atac_trans_ovlp_contingency.xlsx" )

matched.genes <- DD.tab %>% data.frame()
p <- ggplot(matched.genes, aes(x = Var2, y = Var1, fill = Freq)) +
  geom_tile(color = "white",
            lwd = 1,
            linetype = 1) + 
  geom_text(aes(label = Freq), color = "white", size = 8, fontface = "bold") +
  theme_bw() +
  scale_fill_continuous(high = "#132B43", low = "#56B1F7") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 20, face="bold.italic", color = "black")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 20, face="bold.italic", color = "black")) +
  theme(legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 18, face = "bold")) +
  theme(
    plot.title = element_text(size=20, face = "bold.italic", color = 'black'),
    axis.title.x = element_text(size=20, face="bold", hjust = 0.5),
    axis.title.y = element_text(size=20, face="bold", hjust = 0.5)
  ) +
  ylab("") + xlab("")+
  theme(legend.position = "none")+
  guides(color = guide_legend(override.aes = list(size = 7)))

p
ggsave( "../Output/toxo_cdc/ME49_59/figures_paper/matched_DEGs_rna_atac_transitions.pdf", 
        plot = p, height = 8, width = 8, dpi = 300)

