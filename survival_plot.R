library(data.table)
library(dplyr)
library(ggplot2)
library(ggsci)
library(survival)
library(survminer)
library(ComplexHeatmap)
library(circlize)
library(viridis)

setwd('/media/ardadurmaz/HelperWin/Research/splicing/')

obs_ft <- fread('data/combined_annotations.tsv')
embedd <- data.matrix(fread('data/VAE_Embedding.tsv.gz'))
os_ft <- fread('data/combined_annotations_os.tsv.gz', header=TRUE, sep='\t')
mds_clusters <- c('C-1', 'C-2', 'C-3', 'C-7', 'C-8', 'C-11', 'C-12', 'C-13', 'C-14')
aml_clusters <- c('C-0', 'C-4', 'C-5', 'C-6', 'C-9', 'C-10')
#ipssr_data <- fread('~/ipssm_validation/AssignedScore_with_survival_added.csv')
#ipssr_data <- fread('~/ipssm_validation/AssignedScore.csv')
gender <- fread('data/2018_11_01_finalsamplelist.txt')[,1:4]
colnames(gender) <- c('id', 'gender', 'age', 'disease')
os_ft$sex <- ifelse(gender$gender == 1, 'Female', ifelse(gender$gender == 2, 'Male', 'Unknown'))[match(os_ft$id, table=gender$id)]
os_ft$age <- gender$age[match(os_ft$id, table=gender$id)]

# Treatment Info 
treat <- readxl::read_xlsx('data/VV_MergedData.xlsx')
idx <- os_ft$os > 75
os_ft$event[idx] <- 0
os_ft$os[idx] <- 75
os_ft$treat <- treat$threrapy_VV[match(os_ft$id, table=treat$id)]


os_ft <- os_ft %>%
  filter(treat %in% c('HMA', 'Induction chemotherapy/chemotherapy drugs', 'Supportive'))
surv_ft <- as.data.frame(embedd)
colnames(surv_ft) <- paste0('Dim', 1:ncol(embedd))
surv_ft <- surv_ft[match(os_ft$id, table=obs_ft$id),]
surv_ft$age <- os_ft$age
surv_ft$sex <- os_ft$sex
surv_ft$treat <- os_ft$treat
surv_ft$os <- os_ft$os
surv_ft$event <- os_ft$event

library(plyr)
spl_obs <- surv_ft %>% 
  timeSplitter(by = 5,
               event_var = "event",
               event_start_status = 0,
               time_var = "os",
               time_related_vars = c("Dim2", "Dim3"))

cph_fit <- coxph(Surv(os, event)~tt(Dim2)+., data=surv_ft)
surv_ft_tt <- survSplit(Surv(os, event)~., data=surv_ft, cut=c(35), episode='tgroup', id='id')
cph_fit_tt <- coxph(Surv(tstart, os, event)~.+tgroup:Dim2+tgroup:Dim5+tgroup:Dim8+tgroup:Dim18+tgroup:Dim26+tgroup:Dim27, data=surv_ft_tt)
cph_fit_tt <- coxph(Surv(tstart, os, event)~Dim1+Dim2+Dim3+Dim4+Dim5+Dim6+Dim7+Dim8+Dim9+Dim10+Dim11+Dim12+Dim13+Dim14+Dim15+Dim16+Dim17+Dim18+Dim19+Dim20+Dim21+Dim22+Dim23+Dim24+Dim25+Dim26+Dim27+Dim28+Dim29+Dim30+Dim31+Dim32+tgroup:Dim2+tgroup:Dim5+tgroup:Dim8+tgroup:Dim18+tgroup:Dim26+tgroup:Dim27+treat+sex+age, data=surv_ft_tt)

cox.zph(cph_fit_tt)


clust_size <- os_ft %>% group_by(Cluster) %>% summarise(Count=n())
os_ft <- os_ft %>%
  filter(Cluster %in% clust_size$Cluster[clust_size$Count > 10])


os_ft$Cluster <- factor(os_ft$Cluster, levels=c('C-12', paste0('C-', c(1, 11, 13, 14, 2, 3, 7, 8))))
cph_fit <- coxph(Surv(os, event)~tt(Dim-2)., data=surv_ft)
cph_fit <- coxph(Surv(os, event)~Cluster+sex+age+treat, data=os_ft)
temp <- summary(cph_fit_tt)
plot_ft <- data.table('Cluster' = rownames(temp$conf.int),
                      'Coef' = log(temp$conf.int[,1]),
                      'low' = log(temp$conf.int[,3]),
                      'upp' = log(temp$conf.int[,4])) %>%
  filter(!is.na(Coef))
p <- ggplot(plot_ft, aes(x=Coef, y=Cluster)) +
  geom_point(color='dodgerblue') +
  geom_segment(aes(x=low, xend=upp, y=Cluster, yend=Cluster), color='dodgerblue') +
  geom_vline(xintercept = 0.0, linetype='dashed') +
  theme_minimal() +
  xlab('log(HR)') +
  ylab('Features') +
  ggtitle('Relative to HMA/Female') +
  theme(axis.title = element_text(face='bold', size=10),
        plot.title = element_text(face='bold', size=12))
ggsave(p, filename = 'CoxPH_Embedding.pdf', width = 6, height = 8)

res_ft <- surv_ft[,33:37]
res_ft$diag <- obs_ft$type[match(os_ft$id, table=obs_ft$id)]
res_ft$id <- os_ft$id
res_ft$clust <- obs_ft$Cluster[match(os_ft$id, table=obs_ft$id)]

fwrite(res_ft, file='SurvFT.tsv', col.names = TRUE, row.names = FALSE, append = FALSE, quote = FALSE, sep='\t')

risk_ft <- readxl::read_xlsx('data/Risk.xlsx')
os_ft$risk <- risk_ft$Risk[match(os_ft$id, table=risk_ft$array_ID)]
os_ft <- na.omit(os_ft)
table(os_ft$Cluster, os_ft$risk)

os_ft <- subset(os_ft, os_ft$Cluster %in% names(table(os_ft$Cluster))[table(os_ft$Cluster)>20], drop=TRUE)
os_ft$Cluster <- factor(os_ft$Cluster, levels = unique(c('C-12', unique(os_ft$Cluster))))
cph_fit <- coxph(Surv(os, event)~Cluster+risk+sex+age, data=os_ft)
temp <- summary(cph_fit)
plot_ft <- data.table('Cluster' = rownames(temp$conf.int),
                      'Coef' = log(temp$conf.int[,1]),
                      'low' = log(temp$conf.int[,3]),
                      'upp' = log(temp$conf.int[,4]))
plot_ft$Cluster <- gsub(pattern='Cluster', replacement = '', plot_ft$Cluster)
plot_ft$Cluster <- factor(plot_ft$Cluster, levels = plot_ft$Cluster[order(plot_ft$Coef, decreasing = TRUE)])

p <- ggplot(plot_ft, aes(x=Coef, y=Cluster)) +
  geom_point(color='dodgerblue') +
  geom_segment(aes(x=low, xend=upp, y=Cluster, yend=Cluster), color='dodgerblue') +
  geom_vline(xintercept = 0.0, linetype='dashed') +
  theme_minimal() +
  xlab('log(HR)') +
  ylab('Clusters') +
  ggtitle('Relative to C-12') +
  theme(axis.title = element_text(face='bold', size=10),
        plot.title = element_text(face='bold', size=12))
ggsave(p, filename = 'Cluster_Forest_Risk.pdf', width = 4, height = 4)

# Check Latent Embedding #
require(randomForestSRC)
surv_ft <- as.data.table(embedd)
colnames(surv_ft) <- paste0('Dim.', 1:ncol(surv_ft))
surv_ft$os <- os_ft$os[match(obs_ft$id, table=os_ft$id)]
surv_ft$event <- os_ft$event[match(obs_ft$id, table=os_ft$id)]
rownames(surv_ft) <- obs_ft$id

imp_res <- sapply(1:10, function(i){
  rf_fit <- rfsrc(Surv(os, event)~., data=surv_ft)
  return(vimp(rf_fit, importance='permute')$importance)
})
colnames(imp_res) <- paste0('Run-', 1:ncol(imp_res))
imp_ft <- melt(imp_res)
med_ft <- imp_ft %>%
  group_by(Var1) %>%
  summarise(med = mean(value))
imp_ft$Var1 <- factor(imp_ft$Var1, levels=med_ft$Var1[order(med_ft$med, decreasing=TRUE)])

p <- ggplot(imp_ft, aes(x=Var1, y=value)) +
  ggtitle('Permutation Importance') +
  ylab('Decrease in Concordance') +
  xlab('Latent Features') +
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p, filename = 'plots/Importance.pdf')


rf_fit <- rfsrc(Surv(os, event)~., data=surv_ft)
p_list <- sapply(paste0('Dim.', c(21,24,1,29)), simplify = FALSE, function(s){
  r_idx <- sample(1:nrow(embedd), size=1)
  local_res <- do.call('rbind', sapply(1:11, simplify = FALSE, function(i){
    return(embedd[r_idx,])
  }))
  colnames(local_res) <- paste0('Dim.', 1:32)
  r_idx <- which(colnames(local_res) == s)
  local_res[,r_idx] <- seq(-5, 5, 1)
  local_pred <- predict(rf_fit, as.data.frame(local_res))$survival
  colnames(local_pred) <- paste0('t.', 1:ncol(local_pred))
  rownames(local_pred) <- paste0('sample.', 1:nrow(local_pred))
  temp_ft <- reshape2::melt(local_pred)
  colnames(temp_ft) <- c('id', 'time', 'surv')
  temp_ft$dim <- local_res[,r_idx][match(temp_ft$id, table=rownames(local_pred))]
  temp_ft$timeFT <- as.numeric(gsub(pattern='^t\\.', replacement = '', temp_ft$time))
  
  p <- ggplot(temp_ft, aes(x=timeFT, y=surv, color=dim, group=id)) +
    ylim(0, 1) +
    xlim(0, 150) +
    ylab('S(t)') +
    xlab('Months') +
    ggtitle(s) +
    geom_step() +
    theme_minimal() +
    scale_color_viridis_c() +
    theme(legend.title = element_blank())
  return(p)
})
require(cowplot)
p_combined <- plot_grid(plotlist = p_list, nrow = 2)
save_plot(p_combined, nrow=2, base_height = 4, base_width = 12, filename = 'plots/RF_Survival_Top.pdf')


surv_pred <- predict(rf_fit, newdata=surv_ft[,1:32])
temp <- surv_pred$survival
colnames(temp) <- paste0('t.', 1:ncol(temp))
rownames(temp) <- obs_ft$id
temp_ft <- reshape2::melt(temp)
colnames(temp_ft) <- c('id', 'time', 'surv')
temp_ft$cluster <- obs_ft$Cluster[match(temp_ft$id, table=obs_ft$id)]
agg_ft <- temp_ft %>%
  group_by(cluster, time) %>%
  summarise(survM = median(surv))
agg_ft$timeFT <- as.numeric(gsub(pattern='^t\\.', replacement = '', agg_ft$time))
agg_ft$type <- ifelse(agg_ft$cluster %in% mds_clusters, 'MDS', 'AML')
agg_ft$cluster <- factor(agg_ft$cluster, levels=paste0('C-',0:14))

man_colors <- obs_ft %>%
  distinct(Cluster, ClusterColor)

p <- ggplot(agg_ft, aes(x=timeFT, y=survM, color=cluster)) +
  ylim(0, 1) +
  xlim(0, 150) +
  ylab('S(t)') +
  xlab('Months') +
  geom_step() +
  #geom_ribbon(aes(ymin=survL, ymax=survU, fill=cluster), alpha=0.3, color=NA, show.legend = FALSE) +
  theme_minimal() +
  scale_color_manual(values=setNames(man_colors$ClusterColor, man_colors$Cluster)) +
  scale_fill_manual(values=setNames(man_colors$ClusterColor, man_colors$Cluster)) +
  facet_wrap(~type) +
  theme(strip.text.x = element_text(face='bold', size=12),
        legend.title = element_blank())
ggsave(p, filename = 'plots/RF_Survival.pdf', width = 8, height = 4)

# Predict #
pred_res <- sapply(unique(obs_ft$Cluster), simplify = FALSE, function(x){
  local_res <- embedd[which(obs_ft$Cluster == x),]
  local_sample <- MASS::mvrnorm(n=100, mu=colMeans(local_res), Sigma = cov(local_res))
  colnames(local_sample) <- paste0('Dim.', 1:32)
  surv_pred <- predict(rf_fit, newdata=as.data.frame(local_sample))
  temp <- surv_pred$survival
  colnames(temp) <- paste0('t.', 1:ncol(temp))
  rownames(temp) <- paste0('sample.', 1:nrow(temp))
  temp_ft <- reshape2::melt(temp)
  colnames(temp_ft) <- c('id', 'time', 'surv')
  temp_ft$cluster <- x
  return(temp_ft)
})
pred_res_ft <- do.call('rbind', pred_res)
pred_res_ft <- pred_res_ft %>%
  group_by(cluster, time) %>%
  summarise(survM = median(surv),
            survU = quantile(surv, p=0.975),
            survL = quantile(surv, p=0.025))

pred_res_ft$timeFT <- as.numeric(gsub(pattern='^t\\.', replacement = '', pred_res_ft$time))
pred_res_ft$col = obs_ft$color[match(pred_res_ft$cluster, table=obs_ft$Cluster)]

man_colors <- obs_ft %>%
  distinct(Cluster, ClusterColor)

pred_res_ft$type <- ifelse(pred_res_ft$cluster %in% mds_clusters, 'MDS', 'AML')

p <- ggplot(pred_res_ft, aes(x=timeFT, y=survM, color=cluster)) +
  ylim(0, 1) +
  xlim(0, 150) +
  geom_step() +
  #geom_ribbon(aes(ymin=survL, ymax=survU, fill=cluster), alpha=0.3, color=NA, show.legend = FALSE) +
  theme_minimal() +
  scale_color_manual(values=setNames(man_colors$ClusterColor, man_colors$Cluster)) +
  scale_fill_manual(values=setNames(man_colors$ClusterColor, man_colors$Cluster)) +
  facet_wrap(~type)


# KM
os_ft$ClusterGroup <- ifelse(os_ft$Cluster %in% mds_clusters, 'MDS', 'AML')
lr_km <- survfit(Surv(os, event)~Cluster, subset(os_ft, os_ft$risk == 'LR' & os_ft$Cluster %in% c('C-12', 'C-11')))
p <- ggsurvplot(lr_km, 
                palette = setNames(man_colors$ClusterColor, paste0('Cluster=', man_colors$Cluster)), 
                risk.table = TRUE,
                conf.int = TRUE,
                pval = TRUE)
pdf('plots/lr_km_comp2.pdf', width = 10, height = 10)
print(p, newpage = FALSE)
dev.off()

aml_km <- survfit(Surv(os, event)~Cluster, subset(os_ft, os_ft$ClusterGroup == 'AML' & os_ft$risk == 'LR'))
mds_km <- survfit(Surv(os, event)~Cluster, subset(os_ft, os_ft$ClusterGroup == 'MDS' & os_ft$risk == 'LR'))

p <- ggsurvplot(mds_km, palette = setNames(man_colors$ClusterColor, paste0('Cluster=', man_colors$Cluster)), risk.table = TRUE)
pdf('plots/mds_km.pdf', width = 10, height = 10)
print(p, newpage = FALSE)
dev.off()


p <- ggsurvplot(mds_km, palette = setNames(man_colors$ClusterColor, paste0('Cluster=', man_colors$Cluster)))
pdf('plots/mds_km.pdf', width = 12, height = 6)
print(p, newpage = FALSE)
dev.off()

require(MASS)


cph_fit <- coxph(Surv(os, event)~., data=surv_ft)
MASS::stepAIC(cph_fit)

cph_fit <- coxph(Surv(os, event) ~  Dim.1 + Dim.2 + Dim.4 + Dim.5 + 
                   Dim.8 + Dim.11 + Dim.12 + Dim.13 + Dim.14 + Dim.17 + Dim.19 + 
                   Dim.21 + Dim.23 + Dim.25 + Dim.29 + Dim.30, data=surv_ft)

fit <- survfit(Surv(os, event) ~ Cluster, data = subset(os_ft, os_ft$Cluster %in% mds_clusters))


ggsurvplot(fit,
           data = os_ft,
           color = 'Cluster',
           palette = setNames(man_colors$ClusterColor, man_colors$Cluster))


fit <- survfit(Surv(os, event) ~ Cluster, data = subset(os_ft, os_ft$Cluster %in% aml_clusters))
man_colors <- os_ft %>%
  filter(Cluster %in% aml_clusters) %>%
  distinct(Cluster, ClusterColor)

ggsurvplot(fit,
           data = os_ft,
           color = 'Cluster',
           palette = setNames(man_colors$ClusterColor, man_colors$Cluster))

  

comp_res <- matrix(0, ncol=15, nrow=15)
for(i in 1:14){
  for(j in (i+1):15){
    local_os <- subset(os_ft, os_ft$Cluster %in% paste0('C-', c(i-1, j-1)) & os_ft$risk == 'HR')
    if(is.null(local_os) | nrow(local_os) < 20 | min(table(local_os$Cluster)) < 20 | length(unique(local_os$Cluster)) < 2){
      comp_res[i,j] <- NA
    }else{
      fit <- survdiff(Surv(os, event) ~ Cluster, data = local_os)
      comp_res[i,j] <- -log10(1-pchisq(fit$chisq, 1))
    }
  }
}

freq_ft <- os_ft %>%
  group_by(Cluster, type) %>%
  summarise(Count=n()) %>%
  mutate(Freq = Count/sum(Count))
freq_mat <- reshape2::dcast(freq_ft, Cluster~type, value.var = 'Freq')
freq_mat <- freq_mat[match(paste0('C-', 0:14), table=freq_mat$Cluster),-1]
freq_mat[is.na(freq_mat)] <- 0
type_colors <- os_ft %>%
  distinct(type, color)

col_fun <- colorRamp2(breaks=seq(0,ceiling(max(comp_res, na.rm = TRUE)),0.001), colors=viridis::viridis(length(seq(0,ceiling(max(comp_res, na.rm = TRUE)),0.001))))
top_annot <- HeatmapAnnotation(foo = anno_barplot(freq_mat, 
                                                  gp = gpar(fill = type_colors$color[match(colnames(freq_mat), table=type_colors$type)],
                                                            col = 'white'), 
                                                  bar_width = 0.85, 
                                                  height = unit(2, "cm"),
                                                  border=FALSE),
                               show_annotation_name = FALSE)
h.ae <- Heatmap(comp_res,
                col = col_fun,
                cell_fun = function(j, i, x, y, width, height, fill) {
                  if(j > i){
                    grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "gray", fill = fill))
                    grid.text(sprintf("%.1f", comp_res[i, j]), x, y, gp = gpar(fontsize = 6))
                  }else{
                    grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "white", fill = 'white'))
                  }
                },
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                show_column_names = TRUE,
                show_row_names = TRUE,
                column_labels = paste0('C-', 0:14),
                row_labels = paste0('C-', 0:14),
                column_names_gp = gpar(fontface='bold', fontsize=8),
                row_names_gp = gpar(fontface='bold', fontsize=8),
                column_names_side = 'top',
                column_names_rot = 45,
                #top_annotation = top_annot,
                na_col = 'white',
                heatmap_legend_param = list(title='-log10(P.Value)', 
                                            title_position='leftcenter-rot', 
                                            legend_height=unit(4,'cm')),
                show_heatmap_legend = TRUE,
                use_raster = FALSE)

lgd = Legend(labels = type_colors$type, legend_gp = gpar(fill = type_colors$color))

pdf('plots/SurvivalHeatmap.pdf', width = 6, height = 4)
draw(h.ae, annotation_legend_list = list(lgd))
dev.off()

pdf('plots/Survival_HR.pdf', width = 6, height = 4)
draw(h.ae)
dev.off()