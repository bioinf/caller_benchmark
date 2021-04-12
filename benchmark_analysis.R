library(ggplot2)
library(reshape2)
library(dplyr)
library(cowplot)
library(colorRamps)
library(lattice)
library(colorspace)
library(ggsci)
library(scales)
library(stringr)

#setwd('')

mycol1 = rgb(110, 224, 255, maxColorValue = 255)
mycol2 = rgb(255, 192, 78, maxColorValue = 255)
mycol3 = rgb(255, 177, 177, maxColorValue = 255)
mycol4 = rgb(221, 221, 221, maxColorValue = 255)
mycol5 = rgb(100, 89, 89, maxColorValue = 255)
mycol6 = rgb(0, 189, 189, maxColorValue = 255)

########################Coverage analysis######################################

cov_data = read.table('GIAB_CDS_GIABhiconf.metrics', sep='\t', header=T)
head(cov_data)
cov_data$TYPE = ifelse(grepl('EXOME', cov_data$SAMPLE), 'WES', 'WGS')
ggplot(cov_data, aes(x=SAMPLE, y=MEAN_BAIT_COVERAGE, fill=TYPE)) +
  geom_bar(stat='identity', col='black', position='dodge') + 
  theme_bw() + theme(panel.grid=element_blank())

ggplot(cov_data, aes(x=SAMPLE, y=PCT_TARGET_BASES_10X, fill=TYPE)) +
  geom_bar(stat='identity', col='black', position='dodge') + 
  theme_bw() + theme(panel.grid=element_blank())

cov_scatter <- ggplot(cov_data, aes(x=MEAN_BAIT_COVERAGE, y=PCT_TARGET_BASES_10X, col=TYPE)) +
  geom_point(size=4) + theme_bw() + scale_color_npg() +
  xlab('Mean CDS coverage (all CDS)') + 
  ylab('% 10x bases')
print(cov_scatter)

var_counts = read.table('PASS_variant_counts_revamp.tsv', sep=' ', header=F)
head(var_counts)
var_counts$sample = sapply(strsplit(var_counts$V1, '_'), function(x) x[1])
var_counts$TYPE = ifelse(str_extract(var_counts$V1, '[A-Z]+OME') == 'EXOME',
                         'WES', 'WGS')
pf_counts <- ggplot(var_counts, aes(x=sample, y=V2/1000, fill=TYPE)) +
  geom_boxplot() + theme_bw() + scale_fill_npg() +
  xlab('Sample') + ylab('PASS variants (x1000)')

plot_grid(cov_scatter, pf_counts, nrow=2)


aggregate(V2~sample+TYPE, var_counts, median)
mean(aggregate(V2~sample, aggregate(V2~sample+TYPE, var_counts, median), 
               IQR)[, 'V2'])


cov_data_nods = read.table('GIAB_CDS_GIABhiconf_nodownsample.metrics', sep='\t', header=T)
head(cov_data_nods)
cov_data_nods$TYPE = ifelse(grepl('EXOME', cov_data_nods$SAMPLE), 'WES', 'WGS')
cov_data_nods = rbind(cov_data_nods, cov_data[cov_data$SAMPLE == 'HG001_EXOME_BWA', ])
cov_data_nods$SAMPLE[11] = 'HG001_DS_EXOME_BWA'
cov_data_nods$TOTAL_READS = round(cov_data_nods$TOTAL_READS/1000000, digits=0)

s1a <- ggplot(cov_data_nods, aes(x=SAMPLE, y=TOTAL_READS, fill=TYPE)) + 
  geom_bar(stat='identity', col='black', position='dodge') +
  theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1))

s1b <- ggplot(cov_data_nods, aes(x=SAMPLE, y=PCT_TARGET_BASES_20X, fill=TYPE)) +
  geom_bar(stat='identity', col='black', position='dodge') +
  theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1)) +
  scale_y_continuous(limits=c(0, 1))

s1c <- ggplot(cov_data_nods, aes(x=SAMPLE, y=MEAN_BAIT_COVERAGE, fill=TYPE)) +
  geom_bar(stat='identity', col='black', position='dodge') +
  theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1))

plot_grid(s1a, s1b, s1c, nrow=2)



###################Reading and formatting benchmarking data####################

bdata = read.table('all_stats_by_strata_rework.tsv', sep='\t', header=T)
bdata$Precision.AM = (bdata$QUERY.TP + bdata$FP.gt)/(bdata$QUERY.TP + bdata$QUERY.FP)
bdata$Recall.AM = (bdata$QUERY.TP + bdata$FP.gt)/(bdata$QUERY.TP + bdata$FP.gt + bdata$TRUTH.FN)
bdata$F1_Score.AM = 2/(1/bdata$Precision.AM + 1/bdata$Recall.AM)
bdata$F1_Score.AM[bdata$METRIC.F1_Score == 0] = 0

hist(bdata$Recall.AM)
hist(bdata$Precision.AM)
ggplot(bdata, aes(x=METRIC.F1_Score, y=F1_Score.AM)) + geom_point()

str(bdata)
unique(bdata$Type)
unique(bdata$Subset)


aliases = data.frame(src = c(unique(bdata$CallerFilter), unique(bdata$Aligner)),
            alias = c('DV', 'FB', 'G1', 'G2', 'GH', 'ST', 'BT', 'BW', 'IS', 'NO'))
bdata$Aligner = sapply(bdata$Aligner, 
        function(x) aliases[aliases$src == as.character(x), 'alias'])
bdata$CallerFilter = sapply(bdata$CallerFilter, 
        function(x) aliases[aliases$src == as.character(x), 'alias'])



generate_label_df <- function(positioning) {
  df <- data.frame(Aligner = rep(unique(bdata$Aligner), each = 6),
                   CallerFilter = rep(unique(bdata$CallerFilter), 4),
                   METRIC.Precision = rep(positioning, 24),
                   METRIC.Recall = rep(positioning, 24),
                   METRIC.F1_Score = rep(positioning, 24),
                   Label = rep(c('DV', 'FB', 'G1', 'G2', 'GH', 'ST')))
  return(df)
}
  
  
#######Overall performance comparison#################
# All exons, with padding

pad_data = bdata[bdata$Subset == "*" & bdata$Subtype == "*", ]
snp_pass_pad = pad_data[pad_data$Type == 'SNP' & pad_data$Filter == 'PASS', ]
snp_all_pad = pad_data[pad_data$Type == 'SNP' & pad_data$Filter == 'ALL', ]
indel_pass_pad = pad_data[pad_data$Type == 'INDEL' & pad_data$Filter == 'PASS', ]
indel_all_pad = pad_data[pad_data$Type == 'INDEL' & pad_data$Filter == 'ALL', ]



# Precision
ggplot(snp_all_pad, aes(x=Aligner, y=METRIC.Precision, fill=CallerFilter)) + 
  geom_boxplot() + scale_fill_simpsons()

prec_snp_pass <- ggplot(snp_pass_pad, aes(x=Aligner, y=METRIC.Precision, fill=CallerFilter)) + 
  geom_boxplot() + scale_fill_simpsons() + theme_bw() + guides(fill=F)

ggplot(indel_all_pad, aes(x=Aligner, y=METRIC.Precision, fill=CallerFilter)) + 
  geom_boxplot() + scale_fill_simpsons()

prec_indel_pass <- ggplot(indel_pass_pad, aes(x=Aligner, y=METRIC.Precision, fill=CallerFilter)) + 
  geom_boxplot() + scale_fill_simpsons() + theme_bw() + guides(fill=F)

plot_grid(prec_snp_pass, prec_indel_pass, nrow=2)



# Recall
ggplot(snp_all_pad, aes(x=Aligner, y=METRIC.Recall, fill=CallerFilter)) + 
  geom_boxplot() + scale_fill_simpsons() + theme_bw() + guides(fill=F)

recall_snp_pass <- ggplot(snp_pass_pad, aes(x=Aligner, y=METRIC.Recall, fill=CallerFilter)) + 
  geom_boxplot() + scale_fill_simpsons() + theme_bw() + guides(fill=F)

ggplot(indel_all_pad, aes(x=Aligner, y=METRIC.Recall, fill=CallerFilter)) + 
  geom_boxplot() + scale_fill_simpsons() + theme_bw() + guides(fill=F)

recall_indel_pass <- ggplot(indel_pass_pad, aes(x=Aligner, y=METRIC.Recall, fill=CallerFilter)) + 
  geom_boxplot() + scale_fill_simpsons() + theme_bw() + guides(fill=F)

plot_grid(recall_snp_pass, recall_indel_pass, nrow=2)


# F1 Score
sap <- ggplot(snp_all_pad, aes(x=Aligner, y=METRIC.F1_Score, fill=CallerFilter)) + 
  geom_boxplot() + scale_y_continuous(limits=c(0.3, 1)) + theme_bw() +
  theme(legend.position="bottom", panel.grid=element_blank()) + 
  geom_label(data=generate_label_df(0.4), aes(label=Label), 
             position=position_dodge(width=0.8), label.size=0.075) +
  scale_fill_simpsons()

spp <- ggplot(snp_pass_pad, aes(x=Aligner, y=METRIC.F1_Score, fill=CallerFilter)) + 
  geom_boxplot() + scale_y_continuous(limits=c(0.3, 1)) + theme_bw() + 
  theme(panel.grid=element_blank()) + guides(fill=F) + 
  geom_label(data=generate_label_df(0.4), aes(label=Label), 
             position=position_dodge(width=0.8), label.size=0.075) +
  scale_fill_simpsons()

iap <- ggplot(indel_all_pad, aes(x=Aligner, y=METRIC.F1_Score, fill=CallerFilter)) + 
  geom_boxplot() + scale_y_continuous(limits=c(0.3, 1)) + theme_bw() +
  theme(legend.position="bottom", panel.grid=element_blank()) + 
  geom_label(data=generate_label_df(0.4), aes(label=Label), 
             position=position_dodge(width=0.8), label.size=0.075) +
  scale_fill_simpsons()

ipp <- ggplot(indel_pass_pad, aes(x=Aligner, y=METRIC.F1_Score, fill=CallerFilter)) + 
  geom_boxplot() + scale_y_continuous(limits=c(0.3, 1)) + theme_bw() +
  theme(panel.grid=element_blank())  + guides(fill=F) + 
  geom_label(data=generate_label_df(0.4), aes(label=Label), 
             position=position_dodge(width=0.8), label.size=0.075) +
  scale_fill_simpsons()

print(sap)
print(iap)
plot_grid(spp, ipp, nrow=2)



# protein-coding CDS only
pc_data = bdata[bdata$Subset == "GRCh37_proteincoding_only.bed.gz" & 
                  bdata$Subtype == "*", ]
snp_pass = pc_data[pc_data$Type == 'SNP' & pc_data$Filter == 'PASS', ]
snp_all = pc_data[pc_data$Type == 'SNP' & pc_data$Filter == 'ALL', ]
indel_pass = pc_data[pc_data$Type == 'INDEL' & pc_data$Filter == 'PASS', ]
indel_all = pc_data[pc_data$Type == 'INDEL' & pc_data$Filter == 'ALL', ]

# Precision
ggplot(snp_all, aes(x=Aligner, y=METRIC.Precision, fill=CallerFilter)) + 
  geom_boxplot() + scale_fill_simpsons() + theme_bw() + guides(fill=F)

prec_snp_pc_pass <- ggplot(snp_pass, aes(x=Aligner, y=METRIC.Precision, fill=CallerFilter)) + 
  geom_boxplot() + scale_fill_simpsons() + theme_bw() + guides(fill=F)

ggplot(indel_all, aes(x=Aligner, y=METRIC.Precision, fill=CallerFilter)) + 
  geom_boxplot() + scale_fill_simpsons() + theme_bw() + guides(fill=F)

prec_indel_pc_pass <- ggplot(indel_pass, aes(x=Aligner, y=METRIC.Precision, fill=CallerFilter)) + 
  geom_boxplot() + scale_fill_simpsons() + theme_bw() + guides(fill=F)

plot_grid(prec_snp_pc_pass, prec_indel_pc_pass, nrow=2)



# Recall
ggplot(snp_all, aes(x=Aligner, y=METRIC.Recall, fill=CallerFilter)) + 
  geom_boxplot() + scale_fill_simpsons() + theme_bw() + guides(fill=F)

recall_snp_pc_pass <- ggplot(snp_pass, aes(x=Aligner, y=METRIC.Recall, fill=CallerFilter)) + 
  geom_boxplot() + scale_fill_simpsons() + theme_bw() + guides(fill=F)

ggplot(indel_all, aes(x=Aligner, y=METRIC.Recall, fill=CallerFilter)) + 
  geom_boxplot() + scale_fill_simpsons() + theme_bw() + guides(fill=F)

recall_indel_pc_ass <- ggplot(indel_pass, aes(x=Aligner, y=METRIC.Recall, fill=CallerFilter)) + 
  geom_boxplot() + scale_fill_simpsons() + theme_bw() + guides(fill=F)

plot_grid(recall_snp_pc_pass, recall_indel_pc_ass, nrow=2)



# F1 Score
sa <- ggplot(snp_all, aes(x=Aligner, y=METRIC.F1_Score, fill=CallerFilter)) + 
  geom_boxplot() + scale_y_continuous(limits=c(0.55, 1)) + theme_bw() +
  theme(legend.position="bottom", panel.grid=element_blank()) + 
  geom_label(data=generate_label_df(0.6), aes(label=Label), 
             position=position_dodge(width=0.75), size=2.25) +
  scale_fill_simpsons()

sp <- ggplot(snp_pass, aes(x=Aligner, y=METRIC.F1_Score, fill=CallerFilter)) + 
  geom_boxplot() + scale_y_continuous(limits=c(0.55, 1)) + theme_bw() + 
  theme(panel.grid=element_blank()) + guides(fill=F) + 
  geom_label(data=generate_label_df(0.6), aes(label=Label), 
             position=position_dodge(width=0.75), size=2.25) + 
  ylab('F1 Score (SNP)') +
  scale_fill_simpsons()

ia <- ggplot(indel_all, aes(x=Aligner, y=METRIC.F1_Score, fill=CallerFilter)) + 
  geom_boxplot() + scale_y_continuous(limits=c(0.55, 1)) + theme_bw() +
  theme(legend.position="bottom", panel.grid=element_blank()) + 
  geom_label(data=generate_label_df(0.6), aes(label=Label), 
             position=position_dodge(width=0.75), size=2.25) +
  scale_fill_simpsons()

ip <- ggplot(indel_pass, aes(x=Aligner, y=METRIC.F1_Score, fill=CallerFilter)) + 
  geom_boxplot() + scale_y_continuous(limits=c(0.55, 1)) + theme_bw() +
  theme(panel.grid=element_blank())  + guides(fill=F) + 
  geom_label(data=generate_label_df(0.6), aes(label=Label), 
             position=position_dodge(width=0.75), size=2.25) + 
  ylab('F1 Score (INDEL)') +
  scale_fill_simpsons()

print(sa)
print(ia)
f1 <- plot_grid(sp, ip, nrow=2)
print(f1)
plot_grid(sa, ia, nrow=2)

#######Pairwise comparison of tool performance#################
drawCoolHM = function(df, p_df){
  myPanel_a <- function(x, y, z, ...) {
    panel.levelplot(x,y,z,...)
    panel.text(x, y,  p_df[cbind(x,y)]) ## use handy matrix indexing
  }
  return(levelplot(df, col.regions=jet, 
                   at=seq(min(df), max(df), length.out=100), 
                   aspect='fill', colorkey=list(width=3, labels=list(cex=1.0)),
                   scales=list(x=list(rot=45)), xlab=list(label=''), 
                   ylab=list(label=''), panel=myPanel_a))
}

mypal = colorRampPalette(c('white', '#ecad2f'))
mypal = colorRampPalette(c('white', '#d50f0dff'))
jet = diverge_hsv(100)
#jet = mypal(100)

aligner_pairwise = matrix(0, nrow=4, ncol=4)
rownames(aligner_pairwise) = unique(snp_pass$Aligner)
colnames(aligner_pairwise) = unique(snp_pass$Aligner)
pvals = matrix("", nrow=4, ncol=4)
rownames(pvals) = unique(snp_pass$Aligner)
colnames(pvals) = unique(snp_pass$Aligner)


for (i in unique(snp_pass$Aligner)) {
  for (j in unique(snp_pass$Aligner)) {
    left = as.numeric(snp_pass[snp_pass$Aligner == i, 'METRIC.F1_Score'])
    right = as.numeric(snp_pass[snp_pass$Aligner == j, 'METRIC.F1_Score'])
    difference = median(left - right)
    if (i == j){
      pvals[i, j] = ""
    } else {
      pval = wilcox.test(left, right, paired=T)$p.value
      aligner_pairwise[i, j] = difference
      pvals[i, j] = ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", 
                                                       ifelse(pval < 0.05, "*", "")))
    }
  }
}

snp_aln <- drawCoolHM(aligner_pairwise, pvals)
print(snp_aln)

# Variant callers

caller_pairwise = matrix(0, nrow=6, ncol=6)
rownames(caller_pairwise) = unique(snp_pass$CallerFilter)
colnames(caller_pairwise) = unique(snp_pass$CallerFilter)
pvals_callers = matrix("", nrow=6, ncol=6)
rownames(pvals_callers) = unique(snp_pass$CallerFilter)
colnames(pvals_callers) = unique(snp_pass$CallerFilter)


for (i in unique(snp_pass$CallerFilter)) {
  for (j in unique(snp_pass$CallerFilter)) {
    left = as.numeric(snp_pass[snp_pass$CallerFilter == i, 'METRIC.F1_Score'])
    right = as.numeric(snp_pass[snp_pass$CallerFilter == j, 'METRIC.F1_Score'])
    difference = median(left - right)
    if (i == j){
      pvals_callers[i, j] = ""
    } else {
      pval = wilcox.test(left, right, paired=T)$p.value
      caller_pairwise[i, j] = difference
      pvals_callers[i, j] = ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", 
                                                       ifelse(pval < 0.05, "*", "")))
    }
  }
}

snp_call <- drawCoolHM(caller_pairwise, pvals_callers)
print(snp_call)


aligner_pairwise_id = matrix(0, nrow=4, ncol=4)
rownames(aligner_pairwise_id) = unique(snp_pass$Aligner)
colnames(aligner_pairwise_id) = unique(snp_pass$Aligner)
pvals_aln_id = matrix("", nrow=4, ncol=4)
rownames(pvals_aln_id) = unique(snp_pass$Aligner)
colnames(pvals_aln_id) = unique(snp_pass$Aligner)


for (i in unique(snp_pass$Aligner)) {
  for (j in unique(snp_pass$Aligner)) {
    left = as.numeric(indel_pass[indel_pass$Aligner == i, 'METRIC.F1_Score'])
    right = as.numeric(indel_pass[indel_pass$Aligner == j, 'METRIC.F1_Score'])
    difference = median(left - right)
    if (i == j){
      pvals_aln_id[i, j] = ""
    } else {
      pval = wilcox.test(left, right, paired=T)$p.value
      aligner_pairwise_id[i, j] = difference
      pvals_aln_id[i, j] = ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", 
                                                       ifelse(pval < 0.05, "*", "")))
    }
  }
}

indel_aln <- drawCoolHM(aligner_pairwise_id, pvals_aln_id)
print(indel_aln)

# Variant callers

caller_pairwise_id = matrix(0, nrow=6, ncol=6)
rownames(caller_pairwise_id) = unique(snp_pass$CallerFilter)
colnames(caller_pairwise_id) = unique(snp_pass$CallerFilter)
pvals_callers_id = matrix("", nrow=6, ncol=6)
rownames(pvals_callers_id) = unique(snp_pass$CallerFilter)
colnames(pvals_callers_id) = unique(snp_pass$CallerFilter)


for (i in unique(snp_pass$CallerFilter)) {
  for (j in unique(snp_pass$CallerFilter)) {
    left = as.numeric(indel_pass[indel_pass$CallerFilter == i, 'METRIC.F1_Score'])
    right = as.numeric(indel_pass[indel_pass$CallerFilter == j, 'METRIC.F1_Score'])
    difference = median(left - right)
    if (i == j){
      pvals_callers_id[i, j] = ""
    } else {
      pval = wilcox.test(left, right, paired=T)$p.value
      caller_pairwise_id[i, j] = difference
      pvals_callers_id[i, j] = ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", 
                                                               ifelse(pval < 0.05, "*", "")))
    }
  }
}

indel_call <- drawCoolHM(caller_pairwise_id, pvals_callers_id)
print(indel_call)

#print(f1)

all_data = pc_data[pc_data$Filter == 'ALL', ]
pass_data = pc_data[pc_data$Filter == 'PASS', ]
all_data$PrecisionGain = pass_data$METRIC.Precision - all_data$METRIC.Precision
all_data$RecallLoss = pass_data$METRIC.Recall - all_data$METRIC.Recall
gain_loss = melt(all_data, id.vars=c('Sample', 'ExperimentType', 'Aligner',
        'CallerFilter', 'Type'), measure.vars = c('PrecisionGain', 'RecallLoss'))
head(gain_loss)
mean_diff = aggregate(value ~ variable + CallerFilter + Type, gain_loss, median)
#se_diff = aggregate(value ~ variable + CallerFilter + Type, gain_loss, 
#                    function(x) IQR(x)/sqrt(length(x)))
#mean_diff$se = se_diff$value
mean_diff$y_min = aggregate(value ~ variable + CallerFilter + Type, gain_loss,
                            function(x) quantile(x)[2])$value
mean_diff$y_max = aggregate(value ~ variable + CallerFilter + Type, gain_loss,
                            function(x) quantile(x)[4])$value
mean_diff$variable = factor(rep(c('Precision', 'Recall'), nrow(mean_diff)/2))

raw_f1s = aggregate(METRIC.F1_Score ~ CallerFilter + Type, all_data, median)
print(raw_f1s)
print(mean_diff)

filtr_stats <- ggplot(mean_diff, aes(x=CallerFilter, y=value, 
          fill=variable, ymin=y_min, ymax=y_max)) +
  geom_bar(stat='identity', col='black', position='dodge') + 
  geom_errorbar(width=0.4, position=position_dodge(width=0.9)) + theme_bw() +
  theme(panel.grid = element_blank(), legend.position='top') +
  ylab('Change in metric\nvalue after filtering') +
  scale_fill_manual(values=c(mycol1, mycol3)) +
  facet_wrap(~Type, nrow=1)
print(filtr_stats)


# filter_boxes <- ggplot(gain_loss, aes(x=CallerFilter, y=value, 
#                                      fill=variable)) +
#   geom_boxplot(outlier.shape=NA, position='dodge') + 
#    theme_bw() +
#   theme(panel.grid = element_blank(), legend.position='top') +
#   ylab('Change in metric\nvalue after filtering') +
#   scale_fill_manual(values=c(mycol1, mycol3)) +
#   facet_wrap(~Type, nrow=1) +
#   scale_y_continuous(limits=c(-0.5, 0.25))
# print(filter_boxes)

top_right <- plot_grid(snp_aln, indel_aln, nrow=2, labels=c('b', 'c'))
top <- plot_grid(f1, top_right, ncol=2, rel_widths = c(1, 0.48), labels=c('a', ''))
bottom <- plot_grid(snp_call, indel_call, filtr_stats, nrow=1,
                    rel_widths = c(1, 1, 0.95), labels=c('d', 'e', 'f'))
plot_grid(top, bottom, nrow=2, rel_heights = c(1, 0.7))



pass_data$Pipeline = paste(pass_data$Aligner, pass_data$CallerFilter, sep='_')
head(pass_data)
aggregated_stat = aggregate(METRIC.F1_Score ~ Pipeline + Type, pass_data, median)
aggregated_stat$Precision = aggregate(METRIC.Precision ~ Pipeline + Type, 
                                      pass_data, median)$METRIC.Precision
aggregated_stat$Recall = aggregate(METRIC.Recall ~ Pipeline + Type, 
                                   pass_data, median)$METRIC.Recall

print(aggregated_stat[order(aggregated_stat$METRIC.F1_Score, decreasing = T), ])
# Best performing is BW_DV, next non-DV solution is NO_GH, then the best for Strelka is IS_ST

aggregated_stat2 = aggregate(METRIC.F1_Score ~ Pipeline + Type + ExperimentType,
                             pass_data, median)
aggregated_stat2$Precision = aggregate(METRIC.Precision ~ Pipeline + Type + ExperimentType, 
                                      pass_data, median)$METRIC.Precision
aggregated_stat2$Recall = aggregate(METRIC.Recall ~ Pipeline + Type + ExperimentType, 
                                   pass_data, median)$METRIC.Recall

print(aggregated_stat2[order(aggregated_stat2$METRIC.F1_Score, decreasing = T), ])


wgs_perf <- aggregated_stat2[aggregated_stat2$ExperimentType == 'GENOME', ]
wes_perf = aggregated_stat2[aggregated_stat2$ExperimentType == "EXOME",]
wgs_perf$WES_F1 = wes_perf$METRIC.F1_Score
wgs_perf$WES_WGS_diff = wgs_perf$METRIC.F1_Score - wgs_perf$WES_F1
aggregate(WES_WGS_diff~Type, wgs_perf, median)


####################Analysis of covariates#####################################

snp_all_aln <- ggplot(snp_pass, aes(x=ExperimentType, y=METRIC.F1_Score, col=ExperimentType)) +
  geom_boxplot() + theme_bw() +
  theme(panel.grid=element_blank(), axis.text.x=element_blank()) +
  facet_grid(vars(Aligner), vars(CallerFilter)) + scale_color_npg()

indel_all_aln <- ggplot(indel_pass, aes(x=ExperimentType, y=METRIC.F1_Score, col=ExperimentType)) +
  geom_boxplot() + theme_bw() + 
  theme(panel.grid=element_blank(), axis.text.x=element_blank()) +
  facet_grid(vars(Aligner), vars(CallerFilter)) + scale_color_npg()

plot_grid(snp_all_aln, indel_all_aln, nrow=2)

ggplot(indel_pass, aes(x=ExperimentType, y=METRIC.F1_Score, col=ExperimentType)) +
  geom_boxplot() + geom_jitter() + theme_bw() + theme(panel.grid=element_blank()) +
  facet_grid(vars(Aligner), vars(CallerFilter))


snp_gve <- ggplot(snp_pass[snp_pass$Aligner == 'NO', ], 
       aes(x=ExperimentType, y=METRIC.F1_Score, col=ExperimentType)) +
  geom_boxplot() + #geom_jitter() + 
  theme_bw() + theme(panel.grid=element_blank(), 
                     axis.text.x=element_blank(),
                     axis.title.x=element_blank()) +
  facet_wrap(~CallerFilter, nrow=1) + scale_color_npg() +
  ylab('F1 Score (SNP)') + guides(color=F)


indel_gve <- ggplot(indel_pass[indel_pass$Aligner == 'NO', ], 
       aes(x=ExperimentType, y=METRIC.F1_Score, col=ExperimentType)) +
  geom_boxplot() + #geom_jitter() + 
  theme_bw() + theme(panel.grid=element_blank(), 
                     axis.text.x=element_blank(),
                     axis.title.x=element_blank()) +
  facet_wrap(~CallerFilter, nrow=1) + scale_color_npg() +
  ylab('F1 Score (INDEL)') + guides(color=F)

gve_f1 <- plot_grid(snp_gve, indel_gve, nrow=2)
print(gve_f1)




# Coverage strata
unique(bdata$Subset)

bdata$genome_cov_level = as.numeric(sapply(bdata$Subset, function(x)
     ifelse(grepl('genome_mean_cov', x), strsplit(x, '_')[[1]][4], NA)))

bdata$exome_cov_level = as.numeric(sapply(bdata$Subset, function(x)
  ifelse(grepl('exome_mean_cov', x), strsplit(x, '_')[[1]][4], NA)))


#hist(bdata$genome_cov_level)
# Genome coverage stratification

genome_cov_strata = bdata[(!(is.na(bdata$genome_cov_level))) & 
                            bdata$ExperimentType == 'GENOME' &
                            bdata$Filter == 'PASS' &
                            bdata$Subtype == '*' &
                            bdata$Aligner == 'NO' &
                            bdata$TRUTH.TOTAL > 10, ]

gcs_aggr = aggregate(METRIC.F1_Score ~ genome_cov_level+CallerFilter + Type,
                     genome_cov_strata, mean)
gcs_aggr$se = aggregate(METRIC.F1_Score ~ genome_cov_level+CallerFilter + Type,
                        genome_cov_strata, function(x) sd(x)/sqrt(length(x)))$METRIC.F1_Score


gcs <- ggplot(gcs_aggr, aes(x=genome_cov_level, y=METRIC.F1_Score, 
                     ymin=METRIC.F1_Score-se, ymax=METRIC.F1_Score+se,
                              col=CallerFilter)) +
  geom_line(lwd=1) + #geom_errorbar(width=0.2, lwd=1.1) + 
  geom_point(size=3) + theme_bw() +
  theme(panel.grid=element_blank()) + scale_color_simpsons() +
  facet_wrap(~Type, nrow=1) + scale_y_continuous(limits=c(0.55, 1)) +
  scale_x_continuous(limits=c(-0.1, 2.1)) +
  ylab('F1 Score')

# Exome coverage stratification

exome_cov_strata = bdata[(!(is.na(bdata$exome_cov_level))) & 
                            bdata$ExperimentType == 'EXOME' &
                            bdata$Filter == 'PASS' &
                            bdata$Subtype == '*' &
                            bdata$Aligner == 'NO' &
                            bdata$TRUTH.TOTAL > 10, ]

ecs_aggr = aggregate(METRIC.F1_Score ~ exome_cov_level+CallerFilter + Type,
                     exome_cov_strata, mean)
ecs_aggr$se = aggregate(METRIC.F1_Score ~ exome_cov_level+CallerFilter + Type,
                        exome_cov_strata, function(x) sd(x)/sqrt(length(x)))$METRIC.F1_Score


ecs <- ggplot(ecs_aggr, aes(x=exome_cov_level, y=METRIC.F1_Score, 
                     ymin=METRIC.F1_Score-se, ymax=METRIC.F1_Score+se,
                     col=CallerFilter)) +
  geom_line(lwd=1) + #geom_errorbar(width=0.2, lwd=1.1) + 
  geom_point(size=3) + theme_bw() +
  theme(panel.grid=element_blank()) + scale_color_simpsons() +
  facet_wrap(~Type, nrow=1) + scale_y_continuous(limits=c(0.55, 1)) +
  scale_x_continuous(limits=c(-0.1, 2.1)) +
  ylab('F1 Score')

cov_strata <- plot_grid(gcs, ecs, nrow=2)
print(cov_strata)

gcs_aggr$ExperimentType = 'GENOME'
colnames(gcs_aggr)[1] = 'normcov'
ecs_aggr$ExperimentType = 'EXOME'
colnames(ecs_aggr)[1] = 'normcov'
all_cov_stats = rbind(gcs_aggr, ecs_aggr)

cov_plot <- ggplot(all_cov_stats, aes(x=normcov, y=METRIC.F1_Score, 
                               ymin=METRIC.F1_Score-se, ymax=METRIC.F1_Score+se,
                               col=CallerFilter)) +
  geom_line(lwd=1) + #geom_errorbar(width=0.2, lwd=1.1) + 
  geom_point(size=2) + theme_bw() +
  theme(panel.grid=element_blank(), legend.position = "bottom") + 
  scale_color_simpsons() +
  facet_grid(vars(Type), vars(ExperimentType)) + 
  scale_y_continuous(limits=c(0.55, 1)) +
  scale_x_continuous(limits=c(-0.1, 2.0)) +
  ylab('F1 Score') + guides(color=F)
print(cov_plot)



# By GC - stratifications provided by GIAB (custom GC-stratified data below)

bdata$GC_pct = ifelse(grepl('GRCh37_GC', bdata$Subset),
  as.numeric(sapply(bdata$Subset, 
  function(x) tail(strsplit(strsplit(as.character(x), 
                                         'GC_0')[[1]][2], 
                                '.bed')[[1]][1], n=1))), NA)
hist(bdata$GC_pct)
unique(bdata$GC_pct)

by_gc = bdata[(!(is.na(bdata$GC_pct))) &
                            bdata$Filter == 'PASS' &
                            bdata$Subtype == '*' &
                            bdata$Aligner == 'NO' &
                            bdata$TRUTH.TOTAL > 10, ]

gc_aggr = aggregate(METRIC.F1_Score ~ ExperimentType + Type + CallerFilter + GC_pct,
                    by_gc, mean)
gc_aggr$se = aggregate(METRIC.F1_Score ~ ExperimentType + Type + CallerFilter + GC_pct,
                       by_gc, function(x) sd(x)/sqrt(length(x)))$METRIC.F1_Score

gc_plot <- ggplot(gc_aggr, aes(x=GC_pct, y=METRIC.F1_Score, 
                            ymin=METRIC.F1_Score-se, ymax=METRIC.F1_Score+se,
                            col=CallerFilter)) +
  geom_line(lwd=1) + #geom_errorbar(width=0.2, lwd=1.1) + 
  geom_point(size=2) + theme_bw() +
  theme(panel.grid=element_blank(), legend.position = "bottom") + 
  scale_color_simpsons() +
  facet_grid(vars(Type), vars(ExperimentType)) + 
  scale_y_continuous(limits=c(0.55, 1)) +
  scale_x_continuous(limits=c(0, 100)) +
  ylab('F1 Score') + guides(color=F)
print(gc_plot)

strat_plots <- plot_grid(cov_plot, gc_plot, nrow=1, labels=c('c', 'd'))
print(strat_plots)

exo_data = read.table('./all_exo_q10_mf.tsv', sep='\t', header=T, row.names=1)
hist(exo_data$GC)
hist(exo_data$MCOV_agilent)


# Stratify by MF

bdata$genome_mf = as.numeric(sapply(bdata$Subset, function(x)
  ifelse(grepl('genome_mean_mf', x), strsplit(x, '_')[[1]][4], NA)))

bdata$exome_mf = as.numeric(sapply(bdata$Subset, function(x)
  ifelse(grepl('exome_mean_mf', x), strsplit(x, '_')[[1]][4], NA)))

genome_mf_strata = bdata[(!(is.na(bdata$genome_mf))) & 
                            bdata$ExperimentType == 'GENOME' &
                            bdata$Filter == 'PASS' &
                            bdata$Subtype == '*' &
                            bdata$Aligner == 'NO' &
                            bdata$TRUTH.TOTAL > 10 &
                            bdata$Type == "SNP", ]

gmfs_aggr = aggregate(METRIC.F1_Score ~ genome_mf+CallerFilter + Type,
                     genome_mf_strata, mean)
gmfs_aggr$se = aggregate(METRIC.F1_Score ~ genome_mf+CallerFilter + Type,
                        genome_mf_strata, function(x) sd(x)/sqrt(length(x)))$METRIC.F1_Score


exome_mf_strata = bdata[(!(is.na(bdata$exome_mf))) & 
                           bdata$ExperimentType == 'EXOME' &
                           bdata$Filter == 'PASS' &
                           bdata$Subtype == '*' &
                           bdata$Aligner == 'NO' &
                           bdata$TRUTH.TOTAL > 10 &
                           bdata$Type == "SNP", ]

emfs_aggr = aggregate(METRIC.F1_Score ~ exome_mf+CallerFilter + Type,
                      exome_mf_strata, mean)
emfs_aggr$se = aggregate(METRIC.F1_Score ~ exome_mf+CallerFilter + Type,
                         exome_mf_strata, function(x) sd(x)/sqrt(length(x)))$METRIC.F1_Score

gmfs_aggr$ExperimentType = 'GENOME'
colnames(gmfs_aggr)[1] = 'MF'
emfs_aggr$ExperimentType = 'EXOME'
colnames(emfs_aggr)[1] = 'MF'
all_mf_stats = rbind(gmfs_aggr, emfs_aggr)


mfs <- ggplot(all_mf_stats, aes(x=MF, y=METRIC.F1_Score, 
                            ymin=METRIC.F1_Score-se, ymax=METRIC.F1_Score+se,
                            col=CallerFilter)) +
  geom_line(lwd=1) + #geom_errorbar(width=0.2, lwd=1.1) + 
  geom_point(size=2) + theme_bw() +
  theme(panel.grid=element_blank(), legend.position = "bottom") + 
  scale_color_simpsons() +
  facet_wrap(~ExperimentType, nrow=2) + scale_y_continuous(limits=c(0.4, 1)) +
  scale_x_continuous(limits=c(-0.1, 0.5)) +
  ylab('F1 Score') + guides(color=F)
print(mfs)

f3_bottom <- plot_grid(strat_plots, mfs, rel_widths = c(1, 0.35), labels=c('', 'e'))
print(f3_bottom)



###Padding 

pad_data = bdata[grepl('vicinity', bdata$Subset), ]
pad_data = pad_data[pad_data$Subtype == '*' &
                    pad_data$Filter == "PASS" &
                    pad_data$Aligner == "NO" &
                    pad_data$Caller %in% c('DV', 'ST', 'GH'), ]
unique(pad_data$Subset)
dist_df = data.frame(sets=unique(pad_data$Subset),
                     distancia=c(0, 100, 125, 150, 25, 50, 75))
pad_data$Distance = sapply(as.character(pad_data$Subset),
       function(x) dist_df[dist_df$sets == x, 'distancia']) 

dist_aggr = aggregate(METRIC.F1_Score~Type+ExperimentType+Distance, pad_data, median)
dist_aggr$se = aggregate(METRIC.F1_Score~Type+ExperimentType+Distance, pad_data, 
                         function(x) sd(x)/sqrt(length(x)))$METRIC.F1_Score

distances <- ggplot(dist_aggr, aes(x=Distance, y=METRIC.F1_Score, 
                                ymin=METRIC.F1_Score-se, ymax=METRIC.F1_Score+se,
                                col=ExperimentType)) +
  geom_line(lwd=1) + #geom_errorbar(width=0.2, lwd=1.1) + 
  geom_point(size=2) + theme_bw() +
  theme(panel.grid=element_blank(), legend.position = "bottom") + 
  scale_color_npg() +
  facet_wrap(~Type, nrow=2) + scale_y_continuous(limits=c(0.75, 1)) +
  scale_x_continuous(limits=c(-1, 151)) +
  ylab('F1 Score') + guides(color=F)
print(distances)

f3_top = plot_grid(gve_f1, distances, nrow=1, rel_widths = c(1, 0.35), labels=c('a', 'b'))
plot_grid(f3_top, f3_bottom, nrow=2)




# Distributions of some values

exones_gc_dist <- ggplot(exo_data, aes(x=GC)) + theme_bw() +
  geom_histogram(bins = 20, col='black', fill='gray')

plot_grid(exones_gc_dist, exones_gc_dist, nrow=2)

exones_mcov <- ggplot(exo_data, aes(x=MCOV_agilent)) + theme_bw() +
  geom_histogram(bins = 20, col='black', fill='gray') +
  scale_x_continuous(limits=c(0, 2))
print(exones_mcov)

genomes_mcov <- ggplot(exo_data, aes(x=MCOV_wgs)) + theme_bw() +
  geom_histogram(bins = 20, col='black', fill='gray') +
  scale_x_continuous(limits=c(0, 2))

plot_grid(exones_mcov, genomes_mcov, nrow=2)

genomes_mf <- ggplot(exo_data, aes(x=MF_wgs)) + theme_bw() +
  geom_histogram(bins = 20, col='black', fill='gray') +
  scale_x_continuous(limits=c(-0.1, 0.4))
print(genomes_mf)
plot_grid(genomes_mf, genomes_mf, nrow=2)



# Other stratifications

strata_to_compare = list('GRCh37_AllTandemRepeats_201to10000bp_slop5',
                         c('GRCh37_proteincoding_plus100.bed.gz',
                           'GRCh37_proteincoding_plus50.bed.gz',
                           'GRCh37_proteincoding_plus25.bed.gz', 
                           'GRCh37_HG001_GIABv3.3.2_complexindel10bp_slop50.bed.gz', 
                           '*'))
for (i in strata_to_compare){
  piece_of_sheet = bdata[bdata$Filter == 'PASS' &
                           bdata$Subtype == '*' &
                           bdata$TRUTH.TOTAL > 10 &
                           bdata$Subset %in% c(i, 'GRCh37_proteincoding_only.bed.gz'), ]
  if (length(unique(piece_of_sheet$Subset)) == 1) {
    print(paste('Stratum', i, 'uselsss'))
    next
  }
  p <- ggplot(piece_of_sheet, aes(x=Type, y=METRIC.F1_Score, fill=Subset)) + 
    geom_boxplot() + scale_y_continuous(limits=c(0, 1)) + theme_bw() + 
    theme(panel.grid=element_blank(), legend.position='top') + #guides(fill=F) + 
    #    geom_label(data=generate_label_df(0.1), aes(label=Label), 
    #               position=position_dodge(width=0.75), size=2.25) + 
    ylab('F1 Score (SNP)') +
    scale_fill_simpsons() +
    facet_grid(vars(Aligner), vars(CallerFilter))
  print(p)
}


strata_to_compare = list('GRCh37_HG001_GIABv3.3.2_complexindel10bp_slop50.bed.gz')
for (i in strata_to_compare){
  piece_of_sheet = bdata[bdata$Filter == 'PASS' &
                           bdata$Subtype == '*' &
                           bdata$TRUTH.TOTAL > 10 &
                           bdata$Subset %in% c(i, '*') &
                           bdata$Type == 'INDEL', ]
  if (length(unique(piece_of_sheet$Subset)) == 1) {
    print(paste('Stratum', i, 'uselsss'))
    next
  }
  p <- ggplot(piece_of_sheet, aes(x=ExperimentType, y=METRIC.F1_Score, fill=Subset)) + 
    geom_boxplot() + theme_bw() + 
    theme(panel.grid=element_blank(), legend.position='top') + #guides(fill=F) + 
    #    geom_label(data=generate_label_df(0.1), aes(label=Label), 
    #               position=position_dodge(width=0.75), size=2.25) + 
    ylab('F1 Score (SNP)') +
    scale_fill_simpsons() +
    facet_grid(vars(Aligner), vars(CallerFilter))
  print(p)
}








########## SI - precision and recall pairwise comparison

## SNPs

# Precision

p_aligner_pairwise = matrix(0, nrow=4, ncol=4)
rownames(p_aligner_pairwise) = unique(snp_pass$Aligner)
colnames(p_aligner_pairwise) = unique(snp_pass$Aligner)
p_pvals = matrix("", nrow=4, ncol=4)
rownames(p_pvals) = unique(snp_pass$Aligner)
colnames(p_pvals) = unique(snp_pass$Aligner)


for (i in unique(snp_pass$Aligner)) {
  for (j in unique(snp_pass$Aligner)) {
    left = as.numeric(snp_pass[snp_pass$Aligner == i, 'METRIC.Precision'])
    right = as.numeric(snp_pass[snp_pass$Aligner == j, 'METRIC.Precision'])
    difference = median(left - right)
    if (i == j){
      pvals[i, j] = ""
    } else {
      pval = wilcox.test(left, right, paired=T)$p.value
      p_aligner_pairwise[i, j] = difference
      p_pvals[i, j] = ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", 
                                                       ifelse(pval < 0.05, "*", "")))
    }
  }
}

p_snp_aln <- drawCoolHM(p_aligner_pairwise, p_pvals)
print(p_snp_aln)

p_caller_pairwise = matrix(0, nrow=6, ncol=6)
rownames(p_caller_pairwise) = unique(snp_pass$CallerFilter)
colnames(p_caller_pairwise) = unique(snp_pass$CallerFilter)
p_pvals_callers = matrix("", nrow=6, ncol=6)
rownames(p_pvals_callers) = unique(snp_pass$CallerFilter)
colnames(p_pvals_callers) = unique(snp_pass$CallerFilter)


for (i in unique(snp_pass$CallerFilter)) {
  for (j in unique(snp_pass$CallerFilter)) {
    left = as.numeric(snp_pass[snp_pass$CallerFilter == i, 'METRIC.Precision'])
    right = as.numeric(snp_pass[snp_pass$CallerFilter == j, 'METRIC.Precision'])
    difference = median(left - right)
    if (i == j){
      pvals_callers[i, j] = ""
    } else {
      pval = wilcox.test(left, right, paired=T)$p.value
      p_caller_pairwise[i, j] = difference
      p_pvals_callers[i, j] = ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", 
                                                               ifelse(pval < 0.05, "*", "")))
    }
  }
}

p_snp_call <- drawCoolHM(p_caller_pairwise, p_pvals_callers)
print(p_snp_call)

# Recall

r_aligner_pairwise = matrix(0, nrow=4, ncol=4)
rownames(r_aligner_pairwise) = unique(snp_pass$Aligner)
colnames(r_aligner_pairwise) = unique(snp_pass$Aligner)
r_pvals = matrix("", nrow=4, ncol=4)
rownames(r_pvals) = unique(snp_pass$Aligner)
colnames(r_pvals) = unique(snp_pass$Aligner)


for (i in unique(snp_pass$Aligner)) {
  for (j in unique(snp_pass$Aligner)) {
    left = as.numeric(snp_pass[snp_pass$Aligner == i, 'METRIC.Recall'])
    right = as.numeric(snp_pass[snp_pass$Aligner == j, 'METRIC.Recall'])
    difference = median(left - right)
    if (i == j){
      pvals[i, j] = ""
    } else {
      pval = wilcox.test(left, right, paired=T)$p.value
      r_aligner_pairwise[i, j] = difference
      r_pvals[i, j] = ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", 
                                                         ifelse(pval < 0.05, "*", "")))
    }
  }
}

r_snp_aln <- drawCoolHM(r_aligner_pairwise, r_pvals)
print(r_snp_aln)

r_caller_pairwise = matrix(0, nrow=6, ncol=6)
rownames(r_caller_pairwise) = unique(snp_pass$CallerFilter)
colnames(r_caller_pairwise) = unique(snp_pass$CallerFilter)
r_pvals_callers = matrix("", nrow=6, ncol=6)
rownames(r_pvals_callers) = unique(snp_pass$CallerFilter)
colnames(r_pvals_callers) = unique(snp_pass$CallerFilter)


for (i in unique(snp_pass$CallerFilter)) {
  for (j in unique(snp_pass$CallerFilter)) {
    left = as.numeric(snp_pass[snp_pass$CallerFilter == i, 'METRIC.Recall'])
    right = as.numeric(snp_pass[snp_pass$CallerFilter == j, 'METRIC.Recall'])
    difference = median(left - right)
    if (i == j){
      pvals_callers[i, j] = ""
    } else {
      pval = wilcox.test(left, right, paired=T)$p.value
      r_caller_pairwise[i, j] = difference
      r_pvals_callers[i, j] = ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", 
                                                                 ifelse(pval < 0.05, "*", "")))
    }
  }
}

r_snp_call <- drawCoolHM(r_caller_pairwise, r_pvals_callers)
print(r_snp_call)



## Indels

# Precision

p_aligner_pairwise_id = matrix(0, nrow=4, ncol=4)
rownames(p_aligner_pairwise_id) = unique(snp_pass$Aligner)
colnames(p_aligner_pairwise_id) = unique(snp_pass$Aligner)
p_pvals_aln_id = matrix("", nrow=4, ncol=4)
rownames(p_pvals_aln_id) = unique(snp_pass$Aligner)
colnames(p_pvals_aln_id) = unique(snp_pass$Aligner)


for (i in unique(snp_pass$Aligner)) {
  for (j in unique(snp_pass$Aligner)) {
    left = as.numeric(indel_pass[indel_pass$Aligner == i, 'METRIC.Precision'])
    right = as.numeric(indel_pass[indel_pass$Aligner == j, 'METRIC.Precision'])
    difference = median(left - right)
    if (i == j){
      pvals_aln_id[i, j] = ""
    } else {
      pval = wilcox.test(left, right, paired=T)$p.value
      p_aligner_pairwise_id[i, j] = difference
      p_pvals_aln_id[i, j] = ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", 
                                                              ifelse(pval < 0.05, "*", "")))
    }
  }
}

p_indel_aln <- drawCoolHM(p_aligner_pairwise_id, p_pvals_aln_id)
print(p_indel_aln)

p_caller_pairwise_id = matrix(0, nrow=6, ncol=6)
rownames(p_caller_pairwise_id) = unique(snp_pass$CallerFilter)
colnames(p_caller_pairwise_id) = unique(snp_pass$CallerFilter)
p_pvals_callers_id = matrix("", nrow=6, ncol=6)
rownames(p_pvals_callers_id) = unique(snp_pass$CallerFilter)
colnames(p_pvals_callers_id) = unique(snp_pass$CallerFilter)


for (i in unique(snp_pass$CallerFilter)) {
  for (j in unique(snp_pass$CallerFilter)) {
    left = as.numeric(indel_pass[indel_pass$CallerFilter == i, 'METRIC.Precision'])
    right = as.numeric(indel_pass[indel_pass$CallerFilter == j, 'METRIC.Precision'])
    difference = median(left - right)
    if (i == j){
      pvals_callers_id[i, j] = ""
    } else {
      pval = wilcox.test(left, right, paired=T)$p.value
      p_caller_pairwise_id[i, j] = difference
      p_pvals_callers_id[i, j] = ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", 
                                                                  ifelse(pval < 0.05, "*", "")))
    }
  }
}

p_indel_call <- drawCoolHM(p_caller_pairwise_id, p_pvals_callers_id)
print(p_indel_call)


# Recall


r_aligner_pairwise_id = matrix(0, nrow=4, ncol=4)
rownames(r_aligner_pairwise_id) = unique(snp_pass$Aligner)
colnames(r_aligner_pairwise_id) = unique(snp_pass$Aligner)
r_pvals_aln_id = matrix("", nrow=4, ncol=4)
rownames(r_pvals_aln_id) = unique(snp_pass$Aligner)
colnames(r_pvals_aln_id) = unique(snp_pass$Aligner)


for (i in unique(snp_pass$Aligner)) {
  for (j in unique(snp_pass$Aligner)) {
    left = as.numeric(indel_pass[indel_pass$Aligner == i, 'METRIC.Recall'])
    right = as.numeric(indel_pass[indel_pass$Aligner == j, 'METRIC.Recall'])
    difference = median(left - right)
    if (i == j){
      pvals_aln_id[i, j] = ""
    } else {
      pval = wilcox.test(left, right, paired=T)$p.value
      r_aligner_pairwise_id[i, j] = difference
      r_pvals_aln_id[i, j] = ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", 
                                                                ifelse(pval < 0.05, "*", "")))
    }
  }
}

r_indel_aln <- drawCoolHM(r_aligner_pairwise_id, r_pvals_aln_id)
print(r_indel_aln)

r_caller_pairwise_id = matrix(0, nrow=6, ncol=6)
rownames(r_caller_pairwise_id) = unique(snp_pass$CallerFilter)
colnames(r_caller_pairwise_id) = unique(snp_pass$CallerFilter)
r_pvals_callers_id = matrix("", nrow=6, ncol=6)
rownames(r_pvals_callers_id) = unique(snp_pass$CallerFilter)
colnames(r_pvals_callers_id) = unique(snp_pass$CallerFilter)


for (i in unique(snp_pass$CallerFilter)) {
  for (j in unique(snp_pass$CallerFilter)) {
    left = as.numeric(indel_pass[indel_pass$CallerFilter == i, 'METRIC.Recall'])
    right = as.numeric(indel_pass[indel_pass$CallerFilter == j, 'METRIC.Recall'])
    difference = median(left - right)
    if (i == j){
      pvals_callers_id[i, j] = ""
    } else {
      pval = wilcox.test(left, right, paired=T)$p.value
      r_caller_pairwise_id[i, j] = difference
      r_pvals_callers_id[i, j] = ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", 
                                                                    ifelse(pval < 0.05, "*", "")))
    }
  }
}

r_indel_call <- drawCoolHM(r_caller_pairwise_id, r_pvals_callers_id)
print(r_indel_call)

plot_grid(p_snp_aln, p_indel_aln, p_snp_call, p_indel_call,
          r_snp_aln, r_indel_aln, r_snp_call, r_indel_call, 
          nrow=4, rel_heights = c(0.7, 1, 0.7, 1))


## Figure S6 - exon vicinity for all callers/aligners

pad_data_all = bdata[grepl('vicinity', bdata$Subset), ]
pad_data_all = pad_data_all[pad_data_all$Subtype == '*' &
                              pad_data_all$Filter == "PASS", ]
unique(pad_data_all$Subset)
dist_df = data.frame(sets=unique(pad_data_all$Subset),
                     distancia=c(0, 100, 125, 150, 25, 50, 75))
pad_data_all$Distance = sapply(as.character(pad_data$Subset),
                               function(x) dist_df[dist_df$sets == x, 'distancia']) 

dist_aggr_all = aggregate(METRIC.F1_Score~Type+ExperimentType+Distance+Aligner+CallerFilter, pad_data_all, median)
dist_aggr_all$se = aggregate(METRIC.F1_Score~Type+ExperimentType+Distance+Aligner+CallerFilter, 
                             pad_data_all, function(x) sd(x)/sqrt(length(x)))$METRIC.F1_Score

distances_all_snp <- ggplot(dist_aggr_all[dist_aggr_all$Type == 'SNP', ], 
                            aes(x=Distance, y=METRIC.F1_Score, 
                                ymin=METRIC.F1_Score-se, ymax=METRIC.F1_Score+se,
                                col=ExperimentType)) +
  geom_line(lwd=1) + #geom_errorbar(width=0.2, lwd=1.1) + 
  geom_point(size=2) + theme_bw() +
  theme(panel.grid=element_blank(), legend.position = "bottom") + 
  scale_color_npg() +
  facet_grid(vars(Aligner), vars(CallerFilter)) + 
  scale_y_continuous(limits=c(0.65, 1)) +
  scale_x_continuous(limits=c(-1, 151)) +
  ylab('F1 Score') + guides(color=F)
print(distances_all_snp)


distances_all_indel <- ggplot(dist_aggr_all[dist_aggr_all$Type == 'INDEL', ], 
                              aes(x=Distance, y=METRIC.F1_Score, 
                                  ymin=METRIC.F1_Score-se, ymax=METRIC.F1_Score+se,
                                  col=ExperimentType)) +
  geom_line(lwd=1) + #geom_errorbar(width=0.2, lwd=1.1) + 
  geom_point(size=2) + theme_bw() +
  theme(panel.grid=element_blank(), legend.position = "bottom") + 
  scale_color_npg() +
  facet_grid(vars(Aligner), vars(CallerFilter)) + 
  scale_y_continuous(limits=c(0.4, 1)) +
  scale_x_continuous(limits=c(-1, 151)) +
  ylab('F1 Score') + guides(color=F)
print(distances_all_indel)

plot_grid(distances_all_snp, distances_all_indel, nrow=2)


## FIgure S7/S8/S9 - stratification for all calles/aligners

genome_cov_strata_all = bdata[(!(is.na(bdata$genome_cov_level))) & 
                            bdata$ExperimentType == 'GENOME' &
                            bdata$Filter == 'PASS' &
                            bdata$Subtype == '*' &
                            bdata$TRUTH.TOTAL > 10, ]

gcs_aggr_all = aggregate(METRIC.F1_Score~genome_cov_level+CallerFilter+Type+Aligner,
                     genome_cov_strata_all, mean)
gcs_aggr_all$se = aggregate(METRIC.F1_Score~genome_cov_level+CallerFilter+Type+Aligner,
                        genome_cov_strata_all, function(x) sd(x)/sqrt(length(x)))$METRIC.F1_Score


# Exome coverage stratification

exome_cov_strata_all = bdata[(!(is.na(bdata$exome_cov_level))) & 
                           bdata$ExperimentType == 'EXOME' &
                           bdata$Filter == 'PASS' &
                           bdata$Subtype == '*' &
                           bdata$TRUTH.TOTAL > 10, ]

ecs_aggr_all = aggregate(METRIC.F1_Score~exome_cov_level+CallerFilter+Type+Aligner,
                     exome_cov_strata_all, mean)
ecs_aggr_all$se = aggregate(METRIC.F1_Score~exome_cov_level+CallerFilter+Type+Aligner,
                        exome_cov_strata_all, function(x) sd(x)/sqrt(length(x)))$METRIC.F1_Score

gcs_aggr_all$ExperimentType = 'WGS'
colnames(gcs_aggr_all)[1] = 'normcov'
ecs_aggr_all$ExperimentType = 'WES'
colnames(ecs_aggr_all)[1] = 'normcov'
all_cov_stats_all = rbind(gcs_aggr_all, ecs_aggr_all)

cov_plot_snp_all <- ggplot(all_cov_stats_all[all_cov_stats_all$Type == 'SNP', ],
                           aes(x=normcov, y=METRIC.F1_Score, 
                                      ymin=METRIC.F1_Score-se, ymax=METRIC.F1_Score+se,
                                      col=CallerFilter)) +
  geom_line(lwd=1) + #geom_errorbar(width=0.2, lwd=1.1) + 
  geom_point(size=2) + theme_bw() +
  theme(panel.grid=element_blank(), legend.position = "bottom") + 
  scale_color_simpsons() +
  facet_grid(vars(Aligner), vars(ExperimentType)) + 
  scale_y_continuous(limits=c(0.55, 1)) +
  scale_x_continuous(limits=c(-0.1, 2.0)) +
  ylab('F1 Score') + guides(color=F)

cov_plot_indel_all <- ggplot(all_cov_stats_all[all_cov_stats_all$Type == 'INDEL', ],
                           aes(x=normcov, y=METRIC.F1_Score, 
                               ymin=METRIC.F1_Score-se, ymax=METRIC.F1_Score+se,
                               col=CallerFilter)) +
  geom_line(lwd=1) + #geom_errorbar(width=0.2, lwd=1.1) + 
  geom_point(size=2) + theme_bw() +
  theme(panel.grid=element_blank(), legend.position = "bottom") + 
  scale_color_simpsons() +
  facet_grid(vars(Aligner), vars(ExperimentType)) + 
  scale_y_continuous(limits=c(0.55, 1)) +
  scale_x_continuous(limits=c(-0.1, 2.0)) +
  ylab('F1 Score') + guides(color=F)

s7a <- plot_grid(cov_plot_snp_all, cov_plot_indel_all, ncol=2)




# By GC - stratifications provided by GIAB (custom GC-stratified data below)


by_gc_all = bdata[(!(is.na(bdata$GC_pct))) &
                bdata$Filter == 'PASS' &
                bdata$Subtype == '*' &
                bdata$TRUTH.TOTAL > 10, ]

gc_aggr_all = aggregate(METRIC.F1_Score ~ ExperimentType + Type + 
                          Aligner + CallerFilter + GC_pct,
                    by_gc_all, mean)
gc_aggr_all$se = aggregate(METRIC.F1_Score ~ ExperimentType + Type +
                             Aligner + CallerFilter + GC_pct,
                       by_gc_all, function(x) sd(x)/sqrt(length(x)))$METRIC.F1_Score

gc_plot_snp_all <- ggplot(gc_aggr_all[gc_aggr_all$Type == 'SNP',],
                          aes(x=GC_pct, y=METRIC.F1_Score, 
                               ymin=METRIC.F1_Score-se, ymax=METRIC.F1_Score+se,
                               col=CallerFilter)) +
  geom_line(lwd=1) + #geom_errorbar(width=0.2, lwd=1.1) + 
  geom_point(size=2) + theme_bw() +
  theme(panel.grid=element_blank(), legend.position = "bottom") + 
  scale_color_simpsons() +
  facet_grid(vars(Aligner), vars(ExperimentType)) + 
  scale_y_continuous(limits=c(0.55, 1)) +
  scale_x_continuous(limits=c(0, 100)) +
  ylab('F1 Score') + guides(color=F)

gc_plot_indel_all <- ggplot(gc_aggr_all[gc_aggr_all$Type == 'INDEL',],
                          aes(x=GC_pct, y=METRIC.F1_Score, 
                              ymin=METRIC.F1_Score-se, ymax=METRIC.F1_Score+se,
                              col=CallerFilter)) +
  geom_line(lwd=1) + #geom_errorbar(width=0.2, lwd=1.1) + 
  geom_point(size=2) + theme_bw() +
  theme(panel.grid=element_blank(), legend.position = "bottom") + 
  scale_color_simpsons() +
  facet_grid(vars(Aligner), vars(ExperimentType)) + 
  scale_y_continuous(limits=c(0.55, 1)) +
  scale_x_continuous(limits=c(0, 100)) +
  ylab('F1 Score') + guides(color=F)

s7b <- plot_grid(gc_plot_snp_all, gc_plot_indel_all, ncol=2)

s7ab <- plot_grid(s7a, s7b, nrow=2, labels=c('a', 'b'))


# Stratify by MF

genome_mf_strata_all = bdata[(!(is.na(bdata$genome_mf))) & 
                           bdata$ExperimentType == 'GENOME' &
                           bdata$Filter == 'PASS' &
                           bdata$Subtype == '*' &
                           bdata$TRUTH.TOTAL > 10 &
                           bdata$Type == "SNP", ]

gmfs_aggr_all = aggregate(METRIC.F1_Score~genome_mf+CallerFilter+Type+Aligner,
                      genome_mf_strata_all, mean)
gmfs_aggr_all$se = aggregate(METRIC.F1_Score~genome_mf+CallerFilter+Type+Aligner,
                         genome_mf_strata_all, function(x) sd(x)/sqrt(length(x)))$METRIC.F1_Score


exome_mf_strata_all = bdata[(!(is.na(bdata$exome_mf))) & 
                          bdata$ExperimentType == 'EXOME' &
                          bdata$Filter == 'PASS' &
                          bdata$Subtype == '*' &
                          bdata$TRUTH.TOTAL > 10 &
                          bdata$Type == "SNP", ]

emfs_aggr_all = aggregate(METRIC.F1_Score~exome_mf+CallerFilter+Type+Aligner,
                      exome_mf_strata_all, mean)
emfs_aggr_all$se = aggregate(METRIC.F1_Score~exome_mf+CallerFilter+Type+Aligner,
                         exome_mf_strata_all, function(x) sd(x)/sqrt(length(x)))$METRIC.F1_Score

gmfs_aggr_all$ExperimentType = 'WGS'
colnames(gmfs_aggr_all)[1] = 'MF'
emfs_aggr_all$ExperimentType = 'WES'
colnames(emfs_aggr_all)[1] = 'MF'
all_mf_stats_all = rbind(gmfs_aggr_all, emfs_aggr_all)


s7c <- ggplot(all_mf_stats_all[all_mf_stats_all$Type == 'SNP', ], 
              aes(x=MF, y=METRIC.F1_Score, 
                                ymin=METRIC.F1_Score-se, ymax=METRIC.F1_Score+se,
                                col=CallerFilter)) +
  geom_line(lwd=1) + #geom_errorbar(width=0.2, lwd=1.1) + 
  geom_point(size=2) + theme_bw() +
  theme(panel.grid=element_blank(), legend.position = "bottom") + 
  scale_color_simpsons() +
  facet_grid(vars(Aligner), vars(ExperimentType)) + 
  scale_y_continuous(limits=c(0.25, 1)) +
  scale_x_continuous(limits=c(-0.1, 0.5)) +
  ylab('F1 Score') + guides(color=F)
print(s7c)

pholder <- ggplot(mtcars, aes(x = wt, y = mpg)) + geom_blank() +
  theme_bw() + theme(panel.grid = element_blank(),
                     axis.text = element_blank(), 
                     axis.title = element_blank(),
                     axis.ticks = element_blank(),
                     panel.border = element_blank())
s7cc <- plot_grid(s7c, pholder, nrow=2)

plot_grid(s7ab, s7cc, ncol=2, rel_widths = c(1, 0.5), labels=c('', 'c'))


# Formatting Table S1
pc_pass_data = pc_data[pc_data$Filter == 'PASS', c('Sample', 'ExperimentType',
                                                   'Aligner', 'CallerFilter',
                                                   'Type', 'METRIC.F1_Score',
                                                   'METRIC.Precision', 
                                                   'METRIC.Recall')]
write.table(pc_pass_data, file='main_pipeline_stats.tsv', sep='\t', row.names=F,
            quote=F)
