library(data.table)
library(dplyr)
library(msigdbr)
library(ComplexHeatmap)
library(circlize)

setwd('/home/axd497/MDS_Splicing/')
expr_feat <- fread('data/expr_feat.tsv.gz', header=TRUE, sep='\t')
expr_mat <- data.matrix(fread('data/expr_mat.tsv.gz', header=FALSE, sep='\t'))
obs_ft <- fread('data/combined_annotations.tsv', header=TRUE, sep='\t')

require(limma)
# comp_ft <- obs_ft %>%
#   filter(Cluster %in% paste0('C-', c(5,9,0,1,11,3,14,12)))
#comp_ft$Group <- ifelse(comp_ft$Cluster %in% paste0('C-', c(5,9,0)), 'AML', 'MDS')
#comp_ft$GroupType <- paste(comp_ft$Group, comp_ft$SFType, sep='_')

comp_ft <- obs_ft
comp_ft$Group <- ifelse(comp_ft$Cluster %in% paste0('C-', c(5,9,0,1,11,12,3,14)), 'Progressed', 'Control')
design_mat <- model.matrix(~0+Group, data=comp_ft)
colnames(design_mat) <- gsub(pattern='Group', replacement = '', colnames(design_mat))
contrast_mat <- makeContrasts(Diff = Progressed-Control,
                              levels = design_mat)


comp_ft$GroupType <- paste(comp_ft$Group, comp_ft$type, comp_ft$SFType, sep='_')
comp_ft <- comp_ft %>%
  filter(GroupType %in% c('Control_AML_NoSF', 'Control_AML_SF', 'Control_MDS_NoSF', 'Control_MDS_SF', 'Progressed_AML_NoSF', 'Progressed_AML_SF', 'Progressed_MDS_NoSF', 'Progressed_MDS_SF'))
design_mat <- model.matrix(~0+GroupType, data=comp_ft)
colnames(design_mat) <- gsub(pattern='GroupType', replacement = '', colnames(design_mat))


contrast_mat <- makeContrasts(AML_SF_Diff = Progressed_AML_SF-Control_AML_SF,
                              AML_NoSF_Diff = Progressed_AML_NoSF-Control_AML_NoSF,
                              MDS_SF_Diff = Progressed_MDS_SF-Control_MDS_SF,
                              MDS_NoSF_Diff = Progressed_MDS_NoSF-Control_MDS_NoSF,
                              AML_MDS_SF_Diff = ((Progressed_AML_SF-Control_AML_SF)-(Progressed_MDS_SF-Control_MDS_SF)),
                              AML_MDS_NoSF_Diff = ((Progressed_AML_NoSF-Control_AML_NoSF)-(Progressed_MDS_NoSF-Control_MDS_NoSF)),
                              SF_NoSF_Diff = ((Progressed_AML_SF+Progressed_MDS_SF)-(Control_AML_SF+Control_MDS_SF))-((Progressed_AML_NoSF+Progressed_MDS_NoSF)-(Control_AML_NoSF+Control_MDS_NoSF)),
                              AML_SF_NoSF_Diff = (Progressed_AML_SF-Control_AML_SF)-(Progressed_AML_NoSF-Control_AML_NoSF),
                              MDS_SF_NoSF_Diff = (Progressed_MDS_SF-Control_MDS_SF)-(Progressed_MDS_NoSF-Control_MDS_NoSF),
                              AML_MDS_Diff = ((Progressed_AML_SF-Control_AML_SF)-(Progressed_AML_NoSF-Control_AML_NoSF)) - ((Progressed_MDS_SF-Control_MDS_SF)-(Progressed_MDS_NoSF-Control_MDS_NoSF)),
                              levels = design_mat)


local_expr <- t(expr_mat[match(comp_ft$id, table=obs_ft$id),])
local_expr <- normalizeBetweenArrays(local_expr, method = 'cyclicloess')
rownames(local_expr) <- expr_feat$gene

fit <- lmFit(local_expr, design=design_mat)
cfit <- contrasts.fit(fit, contrast_mat)
efit <- eBayes(cfit)
lres_sf <- as.data.frame(topTable(efit, coef = 8, number = Inf, adjust.method = 'fdr', sort.by = 'none'))
lres_nosf <- as.data.frame(topTable(efit, coef = 9, number = Inf, adjust.method = 'fdr', sort.by = 'none'))
lres <- as.data.frame(topTable(efit, coef = 1, number = Inf, adjust.method = 'fdr', sort.by = 'none'))


tfit <- treat(efit, fc=1.2)
tres_sf <- as.data.frame(topTreat(tfit, coef = 8, number = Inf, adjust.method = 'fdr', sort.by = 'none'))
tres_nosf <- as.data.frame(topTreat(tfit, coef = 9, number = Inf, adjust.method = 'fdr', sort.by = 'none'))

fwrite(lres, file = 'data/deg_res.tsv', col.names = TRUE, row.names = TRUE, sep='\t', append=FALSE, quote=FALSE)

# Plot Volcano #
p <- ggplot(lres, aes(x=logFC, y=-log10(adj.P.Val))) +
  xlab('logFC') +
  ylab('-log10(PValue)') +
  geom_point() +
  geom_vline(xintercept = 0.0, linetype='dashed') +
  theme_minimal() +
  theme(axis.title = element_text(face='bold', size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(p, filename = 'plots/VolcanoPlot.pdf', width = 6, height = 6)

require(msigdbr)
require(fgsea)

gene_set <- msigdbr('Homo sapiens', category='C2', subcategory = 'CP:REACTOME')
gene_set_ft <- split(gene_set$human_gene_symbol, gene_set$gs_name)
#gsea_res_nosf <- fgsea(gene_set_ft, sort(setNames(lres_nosf$logFC, rownames(lres_nosf))))
#gsea_res_sf <- fgsea(gene_set_ft, sort(setNames(lres_sf$logFC, rownames(lres_sf))))
gsea_res <- fgsea(gene_set_ft, sort(setNames(lres$logFC, rownames(lres))))
sig_res <- gsea_res %>%
  filter(padj < 0.05)
fwrite(sig_res, file='data/fgsea_res.tsv', col.names = TRUE, row.names = FALSE, sep='\t', append=FALSE, quote=FALSE)


# Plot Pathway Heatmaps #
for(path in c('REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION', 'REACTOME_ADAPTIVE_IMMUNE_SYSTEM', 'REACTOME_RESPIRATORY_ELECTRON_TRANSPORT')){
  local_sig_res <- sig_res$leadingEdge[which(sig_res$pathway == path)][[1]]
  if(length(local_sig_res) > 50){
    local_sig_res <- sample(local_sig_res, size=50, replace = FALSE)
  }
  sig_expr <- local_expr[rownames(local_expr) %in% local_sig_res,]
  sig_expr_agg <- sapply(split(1:nrow(obs_ft), obs_ft$Cluster), function(i){
    rowMeans(sig_expr[,i])
  })
  
  filt_expr <- t(scale(t(sig_expr_agg)))
  thr <- ceiling(max(abs(filt_expr)))
  col_fun <- colorRamp2(breaks=seq(-thr,thr,0.01), colors=colorRampPalette(c("dodgerblue", "white", "firebrick"))(length(seq(-thr, thr, 0.01))))
  
  h <- Heatmap(filt_expr,
               col = col_fun,
               cluster_rows = TRUE,
               cluster_columns = TRUE,
               show_column_names = TRUE,
               show_row_names = TRUE,
               show_column_dend = FALSE,
               show_row_dend = FALSE,
               na_col = 'white',
               heatmap_legend_param = list(title='Scaled Expression', 
                                           title_position='leftcenter-rot', 
                                           legend_height=unit(4,'cm')),
               show_heatmap_legend = TRUE,
               use_raster = FALSE,
               column_names_gp = gpar(fontsize=8, fontface='bold'),
               row_names_gp = gpar(fontsize=6))
  png(sprintf('plots/Heatmap_%s.png', path), width = 3, height = 4, units='in', res=300)
  draw(h)
  dev.off()
}


sig_res$pathway <- make.unique(stringr::str_trunc(gsub(pattern='^REACTOME_', replacement = '', sig_res$pathway), width = 45, side = 'center'))

sig_res$pathway <- factor(sig_res$pathway, levels=sig_res$pathway[order(sig_res$NES, decreasing = TRUE)])
p <- ggplot(sig_res, aes(x=pathway, y=NES, fill=NES)) +
  xlab('Reactome Pathways') +
  ylab('Normalized Enrichment Score') +
  geom_col() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=4),
        axis.title = element_text(face='bold', size=12),
        legend.position = 'none',
        plot.margin = margin(1.75, 0.25, 0.25, 1.75, unit='cm')) +
  scale_fill_gradient2(low='dodgerblue', mid='white', high='firebrick')
ggsave(p, filename = 'plots/FigureS7_Updated.pdf', width = 8, height = 4)


fgsea::plotEnrichment(gene_set_ft[['REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION']],
                      sort(setNames(lres$logFC, rownames(lres))))


fgsea::plotEnrichment(gene_set_ft[['REACTOME_ADAPTIVE_IMMUNE_SYSTEM']],
                      sort(setNames(lres$logFC, rownames(lres))))


fgsea::plotEnrichment(gene_set_ft[['REACTOME_RESPIRATORY_ELECTRON_TRANSPORT']],
                      sort(setNames(lres$logFC, rownames(lres))))

View(gsea_res_nosf)
View(gsea_res_sf)


aml_leading_genes <- Reduce('union', gsea_res_sf$leadingEdge[order(gsea_res_sf$padj, decreasing = FALSE)[1:2]])
mds_leading_genes <- Reduce('union', gsea_res_nosf$leadingEdge[order(gsea_res_nosf$padj, decreasing = FALSE)[1:2]])
leading_genes <- union(aml_leading_genes, mds_leading_genes)
leading_genes <- Reduce('union', gsea_res$leadingEdge[order(gsea_res$padj, decreasing = FALSE)[1:3]])

filt_expr <- t(scale(t(local_expr[match(leading_genes, table=rownames(local_expr)),])))
thr <- ceiling(max(abs(filt_expr)))
col_fun <- colorRamp2(breaks=seq(-thr,thr,0.01), colors=colorRampPalette(c("dodgerblue", "white", "firebrick"))(length(seq(-thr, thr, 0.01))))

h <- Heatmap(filt_expr,
             col = col_fun,
             cluster_rows = TRUE,
             cluster_columns = TRUE,
             show_column_names = FALSE,
             show_row_names = TRUE,
             show_column_dend = FALSE,
             show_row_dend = FALSE,
             na_col = 'white',
             heatmap_legend_param = list(title='Scaled Expression', 
                                         title_position='leftcenter-rot', 
                                         legend_height=unit(4,'cm')),
             show_heatmap_legend = TRUE,
             use_raster = FALSE,
             column_split = data.table('Disease' = comp_ft$type,
                                       'SF' = comp_ft$SFType,
                                       'Group' = comp_ft$Group),
             column_title_gp = gpar(fontsize=8, fontface='bold'),
             row_names_gp = gpar(fontsize=4),
             cluster_column_slices = FALSE,
             column_title= "%s\n%s\n%s")
draw(h)
