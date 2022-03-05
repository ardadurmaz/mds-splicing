library(data.table)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

setwd('Research/splicing/')

psi_mat <- data.matrix(fread('data/psi_mat.tsv.gz', header=FALSE, sep='\t'))
expr_mat <- data.matrix(fread('data/expr_mat.tsv.gz', header=FALSE, sep='\t'))
obs_ft <- fread('data/combined_annotations.tsv', header=TRUE, sep='\t')
embedd = data.matrix(fread('data/VAE_Embedding.tsv.gz', header=FALSE, sep='\t'))
psi_feat <- fread('data/psi_feat.tsv', header=TRUE, sep='\t')
expr_feat <- fread('data/expr_feat.tsv.gz', header=TRUE, sep='\t')

# Initial Filtering #
psi_feat$expr <- colMeans(expr_mat[,match(psi_feat$gene, table=expr_feat$gene)])
idx <- na.omit(which(psi_feat$expr > 5))
psi_feat <- psi_feat[idx,]
psi_mat <- psi_mat[,idx]

psi_cor <- cor(embedd, psi_mat)
expr_cor <- cor(embedd, expr_mat)

psi_top_cor_idx <- Reduce(union, apply(psi_cor, 1, function(x){order(abs(x), decreasing = TRUE)[1:15]}))
expr_top_cor_idx <- Reduce(union, apply(expr_cor, 1, function(x){order(abs(x), decreasing = TRUE)[1:15]}))

mean_embedd <- sapply(split(1:nrow(obs_ft), obs_ft$Cluster), function(x){
  colMeans(embedd[x,])
})
mean_embedd <- t(scale(t(mean_embedd)))
mean_psi <- sapply(split(1:nrow(obs_ft), obs_ft$Cluster), function(x){
  colMeans(psi_mat[x,psi_top_cor_idx])
})
emb_thr <- ceiling(max(abs(mean_embedd)))
mean_expr <- sapply(split(1:nrow(obs_ft), obs_ft$Cluster), function(x){
  colMeans(expr_mat[x,expr_top_cor_idx])
})


colAnnot <- HeatmapAnnotation(PSI=t(scale(t(mean_psi[,match(paste0('C-', c(2,4,10,7,8,13,6,14,3,0,11,12,1,9,5)), table=colnames(mean_psi))]))),
                              col = list(PSI=colorRamp2(breaks=seq(-3,3,0.01), colors=scico::scico(length(seq(-3,3,0.01)), palette = 'vik'))))
rowAnnot <- HeatmapAnnotation(Embedd=mean_embedd[,match(paste0('C-', c(2,4,10,7,8,13,6,14,3,0,11,12,1,9,5)), table=colnames(mean_embedd))], col=list(Embedd=colorRamp2(breaks=seq(-emb_thr,emb_thr,0.01), 
                                                                             colors=colorRampPalette(c('dodgerblue', 'white', 'firebrick'))(length(seq(-emb_thr,emb_thr,0.01))))),which='row')
col_fun = colorRamp2(seq(-1, 1, 0.01), colors = scico::scico(n=length(seq(-1,1,0.01)), palette='berlin'))

# km
hc_res <- hclust(dist(t(psi_cor[,psi_top_cor_idx])), method = 'ward.D2')
clust_res <- cutree(hc_res, k = 11)
col_split <- data.frame('Cluster' = paste0('Cluster-', clust_res))
h <- Heatmap(psi_cor[,psi_top_cor_idx],
             col = col_fun,
             cluster_rows = TRUE,
             cluster_columns = TRUE,
             show_column_dend = TRUE,
             show_row_dend = FALSE,
             show_row_names = TRUE,
             show_column_names = TRUE,
             row_labels = paste0('Dim.', 1:nrow(psi_cor)),
             column_labels = psi_feat$gene[psi_top_cor_idx],
             heatmap_legend_param = list(title = 'Latent-PSI Correlation', 
                                         title_position='leftcenter-rot', 
                                         legend_height=unit(4,'cm')),
             show_heatmap_legend = TRUE,
             left_annotation = rowAnnot,
             top_annotation = colAnnot,
             use_raster = FALSE,
             column_split = col_split,
             column_names_gp = gpar(fontsize=2),
             row_names_gp = gpar(fontsize=8),
             column_title_gp = gpar(fontsize=8))

pdf('plots/top_psi_cor.pdf', width = 16, height = 8)
draw(h)
dev.off()

temp <- psi_feat[psi_top_cor_idx,]
temp$cluster <- col_split$Cluster
fwrite(temp, file ='data/psi_top_cor.tsv', col.names = TRUE, row.names = FALSE, sep='\t', append=FALSE, quote=FALSE)

colAnnot <- HeatmapAnnotation(Expression=t(scale(t(mean_expr[,match(paste0('C-', c(2,4,10,7,8,13,6,14,3,0,11,12,1,9,5)), table=colnames(mean_expr))]))),
                              col = list(Expression=colorRamp2(breaks=seq(-3,3,0.01), colors=scico::scico(length(seq(-3,3,0.01)), palette = 'vik'))))
rowAnnot <- HeatmapAnnotation(Embedd=mean_embedd[,match(paste0('C-', c(2,4,10,7,8,13,6,14,3,0,11,12,1,9,5)), table=colnames(mean_embedd))], col=list(Embedd=colorRamp2(breaks=seq(-emb_thr,emb_thr,0.01), 
                                                                                                                                                                       colors=colorRampPalette(c('dodgerblue', 'white', 'firebrick'))(length(seq(-emb_thr,emb_thr,0.01))))),which='row')
hc_res <- hclust(dist(t(expr_cor[,expr_top_cor_idx])), method = 'ward.D2')
clust_res <- cutree(hc_res, k = 12)
col_split <- data.frame('Cluster' = paste0('Cluster-', clust_res))

col_fun = colorRamp2(seq(-1, 1, 0.01), colors = scico::scico(n=length(seq(-1,1,0.01)), palette='berlin'))
h <- Heatmap(expr_cor[,expr_top_cor_idx],
             col = col_fun,
             cluster_rows = TRUE,
             cluster_columns = TRUE,
             show_column_dend = TRUE,
             show_row_dend = FALSE,
             show_row_names = TRUE,
             show_column_names = TRUE,
             row_labels = paste0('Dim.', 1:nrow(expr_cor)),
             column_labels = expr_feat$gene[expr_top_cor_idx],
             heatmap_legend_param = list(title = 'Latent-Expression Correlation', 
                                         title_position='leftcenter-rot', 
                                         legend_height=unit(4,'cm')),
             show_heatmap_legend = TRUE,
             left_annotation = rowAnnot,
             top_annotation = colAnnot,
             use_raster = FALSE,
             column_split = col_split,
             column_names_gp = gpar(fontsize=2),
             row_names_gp = gpar(fontsize=7),
             column_title_gp = gpar(fontsize=8))


pdf('plots/top_expr_cor.pdf', width = 16, height = 8)
draw(h)
dev.off()

temp <- expr_feat[expr_top_cor_idx,]
temp$cluster <- col_split$Cluster

fwrite(temp, file ='data/expr_top_cor.tsv', col.names = TRUE, row.names = FALSE, sep='\t', append=FALSE, quote=FALSE)
