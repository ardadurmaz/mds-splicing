library(data.table)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(cluster)
library(data.table)
library(dplyr)
library(ggplot2)
library(cluster)

setwd('Research/splicing/')
psi_ft <- data.matrix(fread('data/psi_mat.tsv.gz', header=FALSE, sep='\t'))
expr_ft <- data.matrix(fread('data/expr_mat.tsv.gz', header=FALSE, sep='\t'))
mut_ft <- data.matrix(fread('data/mut_mat.tsv.gz', header=FALSE, sep='\t'))
obs_ft <- fread('data/combined_annotations.tsv', header=TRUE, sep='\t')
embedd = data.matrix(fread('data/VAE_Embedding.tsv.gz', header=FALSE, sep='\t'))
mut_feat <- fread('data/mut_feat.tsv.gz', header=TRUE, sep='\t')
expr_feat <- fread('data/expr_feat.tsv.gz', header=TRUE, sep='\t')
psi_feat <- fread('data/psi_feat.tsv.gz', header=TRUE, sep='\t')

exprSim <- cor(embedd, expr_ft)
psiSim <- cor(embedd, psi_ft)
expr_top <- Reduce('union', sapply(1:32, simplify = FALSE, function(i){
  return(c(order(exprSim[i,], decreasing = TRUE)[1:3],
           order(exprSim[i,], decreasing = FALSE)[1:3]))
}))
psi_top <- Reduce('union', sapply(1:32, simplify = FALSE, function(i){
  return(c(order(psiSim[i,], decreasing = TRUE)[1:3],
           order(psiSim[i,], decreasing = FALSE)[1:3]))
}))




# Heatmap #
man_colors <- obs_ft %>%
  distinct(type, color)

embedd_thr <- ceiling(max(abs(embedd)))

col_fun <- colorRamp2(breaks=seq(0,1,0.01), colors=viridis::inferno(length(seq(0,1,0.01))))
col_fun_ae <- colorRamp2(breaks=seq(-embedd_thr,embedd_thr,0.01), colors=colorRampPalette(c("dodgerblue", "white", "firebrick"))(length(seq(-embedd_thr, embedd_thr, 0.01))))
row_annot <- HeatmapAnnotation("Type"=obs_ft$type,
                               col=list("Type"=setNames(man_colors$color, man_colors$type)), 
                               show_annotation_name = FALSE, 
                               which = 'row',
                               show_legend = TRUE)
h <- Heatmap(embedd,
             col = col_fun_ae,
             cluster_rows = TRUE,
             cluster_columns = TRUE,
             show_column_names = TRUE,
             show_row_names = FALSE,
             show_column_dend = FALSE,
             show_row_dend = TRUE,
             column_labels = paste0('Dim.', 1:32),
             column_names_gp = gpar(fontface='bold', fontsize=4),
             na_col = 'white',
             left_annotation = row_annot,
             heatmap_legend_param = list(title='Latent Components', 
                                         title_position='leftcenter-rot', 
                                         legend_height=unit(4,'cm')),
             show_heatmap_legend = TRUE,
             use_raster = FALSE,
             row_split = data.table('Cluster' = obs_ft$Cluster),
             row_title_gp = gpar(fontsize=12, fontface='bold'),
             column_title = 'Latent Embedding')

h.ae <- Heatmap(mut_ft,
                col = col_fun,
                cluster_rows = TRUE,
                cluster_columns = TRUE,
                show_column_names = TRUE,
                show_row_names = FALSE,
                show_column_dend = FALSE,
                show_row_dend = FALSE,
                column_labels = mut_feat$id,
                column_names_gp = gpar(fontface='bold', fontsize=4),
                na_col = 'white',
                heatmap_legend_param = list(title='Variant Status', 
                                            title_position='leftcenter-rot', 
                                            legend_height=unit(4,'cm')),
                show_heatmap_legend = FALSE,
                use_raster = FALSE,
                column_title = 'Mutations')

expr_thr <- ceiling(max(abs(scale(expr_ft[, expr_top]))))
col_fun_expr <- colorRamp2(breaks=seq(-expr_thr,expr_thr,0.01), colors=scico::scico(length(seq(-expr_thr, expr_thr, 0.01)), palette = 'berlin'))
h.expr <- Heatmap(scale(expr_ft[,expr_top]),
                  col = col_fun_expr,
                  cluster_rows = TRUE,
                  cluster_columns = TRUE,
                  show_column_names = TRUE,
                  show_row_names = FALSE,
                  show_column_dend = FALSE,
                  show_row_dend = FALSE,
                  column_labels = expr_feat$gene[expr_top],
                  column_names_gp = gpar(fontface='bold', fontsize=4),
                  na_col = 'white',
                  heatmap_legend_param = list(title='Scaled Expression', 
                                              title_position='leftcenter-rot', 
                                              legend_height=unit(4,'cm')),
                  show_heatmap_legend = TRUE,
                  use_raster = FALSE,
                  column_title = 'Gene Expression')

col_fun_psi <- colorRamp2(breaks=seq(0,1,0.01), colors=scico::scico(length(seq(0, 1, 0.01)), palette = 'davos'))
h.psi <- Heatmap(psi_ft[,psi_top],
                  col = col_fun_psi,
                  cluster_rows = TRUE,
                  cluster_columns = TRUE,
                  show_column_names = TRUE,
                  show_row_names = FALSE,
                  show_column_dend = FALSE,
                  show_row_dend = FALSE,
                  column_labels = psi_feat$gene[psi_top],
                  column_names_gp = gpar(fontface='bold', fontsize=4),
                  na_col = 'white',
                  heatmap_legend_param = list(title='PSI', 
                                              title_position='leftcenter-rot', 
                                              legend_height=unit(4,'cm')),
                  show_heatmap_legend = TRUE,
                  use_raster = FALSE,
                  column_title = 'PSI')


h.combined <- h + h.ae + h.expr + h.psi

pdf('CombinedHeatmap.pdf', width = 24, height = 12)
draw(h.combined)
dev.off()
