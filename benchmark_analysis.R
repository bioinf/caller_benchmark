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
library(ggrepel)

setwd('/media/array/callers_proj/')

mycol1 = rgb(110, 224, 255, maxColorValue = 255)
mycol2 = rgb(255, 192, 78, maxColorValue = 255)
mycol3 = rgb(255, 177, 177, maxColorValue = 255)
mycol4 = rgb(221, 221, 221, maxColorValue = 255)
mycol5 = rgb(100, 89, 89, maxColorValue = 255)
mycol6 = rgb(0, 189, 189, maxColorValue = 255)

########################Coverage analysis######################################

cov_data = read.table('./All_HS.metrics', sep='\t', header=T)
head(cov_data)
cov_data$TYPE = ifelse(grepl('EXOME', cov_data$SAMPLE) | grepl('RUSZ', cov_data$SAMPLE), 'WES', 'WGS')
ggplot(cov_data, aes(x=SAMPLE, y=MEAN_BAIT_COVERAGE, fill=TYPE)) +
  geom_bar(stat='identity', col='black', position='dodge') + 
  theme_bw() + theme(panel.grid=element_blank())

ggplot(cov_data, aes(x=SAMPLE, y=PCT_TARGET_BASES_10X, fill=TYPE)) +
  geom_bar(stat='identity', col='black', position='dodge') + 
  theme_bw() + theme(panel.grid=element_blank())

cov_scatter <- ggplot(cov_data[1:14, ], aes(x=MEAN_BAIT_COVERAGE, y=PCT_TARGET_BASES_10X, col=TYPE)) +
  geom_point(size=4) + theme_bw() + scale_color_npg() +
  xlab('Mean CDS coverage (all CDS)') + 
  ylab('% 10x bases')
print(cov_scatter)


cov_scatter_ng <- ggplot(cov_data[15:20, ], aes(x=MEAN_BAIT_COVERAGE, y=PCT_TARGET_BASES_10X, col=TYPE)) +
  geom_point(size=4) + theme_bw() + scale_color_npg() +
  xlab('Mean CDS coverage (all CDS)') + 
  ylab('% 10x bases')
print(cov_scatter_ng)


var_counts = read.table('v4_variant_counts.tsv', sep='\t', header=F)
head(var_counts)
var_counts$sample = sapply(strsplit(var_counts$V1, '_'), function(x) x[1])
var_counts$TYPE = ifelse(str_extract(var_counts$V1, '[A-Z]+OME') == 'EXOME',
                         'WES', 'WGS')
colnames(var_counts) = c('sample', 'TYPE', 'aligner', 'caller', 'V2')

pf_counts <- ggplot(var_counts, aes(x=sample, y=V2/1000, fill=TYPE)) +
  geom_boxplot() + theme_bw() + scale_fill_npg() +
  xlab('Sample') + ylab('PASS variants (x1000)')

plot_grid(cov_scatter, pf_counts, nrow=2)

s1b_rev <- pf_counts + scale_y_continuous(limits=c(0, 35)) +
  theme(axis.text.x = element_text(angle=45, hjust=1), legend.position='bottom')


aggregate(V2~sample+TYPE, var_counts, median)
mean(aggregate(V2~sample, aggregate(V2~sample+TYPE, var_counts, median), 
               IQR)[, 'V2'])


cov_data_nods = read.table('./GIAB_CDS_GIABhiconf_nodownsample.metrics', sep='\t', header=T)
head(cov_data_nods)
cov_data_nods$TYPE = ifelse(grepl('EXOME', cov_data_nods$SAMPLE), 'WES', 'WGS')
cov_data_nods = rbind(cov_data_nods, cov_data[cov_data$SAMPLE == 'HG001_EXOME_BWA', ])
cov_data_nods$SAMPLE[11] = 'HG001_DS_EXOME_BWA'
cov_data_nods$TOTAL_READS = round(cov_data_nods$TOTAL_READS/1000000, digits=0)

s1a <- ggplot(cov_data_nods, aes(x=SAMPLE, y=TOTAL_READS, fill=TYPE)) + 
  geom_bar(stat='identity', col='black', position='dodge') +
  theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1))

s1b <- ggplot(cov_data[grepl('HG00', cov_data$SAMPLE), ], 
              aes(x=SAMPLE, y=PCT_TARGET_BASES_20X, fill=TYPE)) +
  geom_bar(stat='identity', col='black', position='dodge') +
  theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1),
                     legend.position='bottom') +
  scale_y_continuous(limits=c(0, 1)) + xlab('Sample') +
  ylab('% CDS bases with 20x coverage')
s1b

s1c <- ggplot(cov_data_nods, aes(x=SAMPLE, y=MEAN_BAIT_COVERAGE, fill=TYPE)) +
  geom_bar(stat='identity', col='black', position='dodge') +
  theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1))

plot_grid(s1a, s1b, s1c, nrow=2)

# Revision Figure S1
plot_grid(s1b, s1b_rev, labels=c('a', 'b'), nrow=1)



###################Reading and formatting benchmarking data####################

#bdata = read.table('./all_stats_by_strata_rework.tsv', sep='\t', header=T)
bdata = read.table('./all_stats_v4_by_strata_rework.tsv', sep='\t', header=T)

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
unique(bdata$CallerFilter)
unique(bdata$Aligner)


aliases = data.frame(src = c(unique(bdata$CallerFilter), unique(bdata$Aligner)),
            alias = c('CL', 'DV', 'FB', 'G1', 'G2', 'GH', 'OF', 'OS', 'ST', 
                      'BT-LOC', 'BT-E2E', 'BW', 'IS', 'NO'))
bdata$Aligner = sapply(bdata$Aligner, 
        function(x) aliases[aliases$src == as.character(x), 'alias'])
bdata$CallerFilter = sapply(bdata$CallerFilter, 
        function(x) aliases[aliases$src == as.character(x), 'alias'])



  
  
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


##### Main analysis ################################################
# protein-coding CDS only
pc_data = bdata[bdata$Subset == "GRCh37_proteincoding_only.bed.gz" & 
                  bdata$Subtype == "*", ]
pc_data = pc_data[pc_data$Aligner != 'BT-LOC', ]
pc_data$Aligner = as.factor(as.character(pc_data$Aligner))
snp_pass = pc_data[pc_data$Type == 'SNP' & pc_data$Filter == 'PASS', ]
snp_all = pc_data[pc_data$Type == 'SNP' & pc_data$Filter == 'ALL', ]
indel_pass = pc_data[pc_data$Type == 'INDEL' & pc_data$Filter == 'PASS', ]
indel_all = pc_data[pc_data$Type == 'INDEL' & pc_data$Filter == 'ALL', ]

generate_label_df <- function(positioning) {
  df <- data.frame(Aligner = rep(unique(pc_data$Aligner), each = 9),
                   CallerFilter = rep(unique(bdata$CallerFilter), 4),
                   METRIC.Precision = rep(positioning, 36),
                   METRIC.Recall = rep(positioning, 36),
                   METRIC.F1_Score = rep(positioning, 36),
                   Label = rep(c('CL', 'DV', 'FB', 'G1', 'G2', 'GH', 'OF', 'OS', 'ST')))
  return(df)
}


# Precision
ggplot(snp_all, aes(x=Aligner, y=METRIC.Precision, fill=CallerFilter)) + 
  geom_boxplot() + scale_fill_simpsons() + theme_bw() + guides(fill=F)

prec_snp_pc_pass <- ggplot(snp_pass, aes(x=Aligner, y=METRIC.Precision, fill=CallerFilter)) + 
  geom_boxplot() + scale_fill_brewer(palette='Set3') + theme_bw() + guides(fill=F)

ggplot(indel_all, aes(x=Aligner, y=METRIC.Precision, fill=CallerFilter)) + 
  geom_boxplot() + scale_fill_simpsons() + theme_bw() + guides(fill=F)

prec_indel_pc_pass <- ggplot(indel_pass, aes(x=Aligner, y=METRIC.Precision, fill=CallerFilter)) + 
  geom_boxplot() + scale_fill_brewer(palette='Set3') + theme_bw() + guides(fill=F)

plot_grid(prec_snp_pc_pass, prec_indel_pc_pass, nrow=2)



# Recall
ggplot(snp_all, aes(x=Aligner, y=METRIC.Recall, fill=CallerFilter)) + 
  geom_boxplot() + scale_fill_simpsons() + theme_bw() + guides(fill=F)

recall_snp_pc_pass <- ggplot(snp_pass, aes(x=Aligner, y=METRIC.Recall, fill=CallerFilter)) + 
  geom_boxplot() + scale_fill_brewer(palette='Set3') + theme_bw() + guides(fill=F)

ggplot(indel_all, aes(x=Aligner, y=METRIC.Recall, fill=CallerFilter)) + 
  geom_boxplot() + scale_fill_simpsons() + theme_bw() + guides(fill=F)

recall_indel_pc_ass <- ggplot(indel_pass, aes(x=Aligner, y=METRIC.Recall, fill=CallerFilter)) + 
  geom_boxplot() + scale_fill_brewer(palette='Set3')+ theme_bw() + guides(fill=F)

plot_grid(recall_snp_pc_pass, recall_indel_pc_ass, nrow=2)



# F1 Score
sa <- ggplot(snp_all, aes(x=Aligner, y=METRIC.F1_Score, fill=CallerFilter)) + 
  geom_boxplot() + scale_y_continuous(limits=c(0.4, 1)) + theme_bw() +
  theme(legend.position="bottom", panel.grid=element_blank()) + 
  geom_label(data=generate_label_df(0.45), aes(label=Label), 
             position=position_dodge(width=0.75), size=2.25) +
  scale_fill_brewer(palette='Set3')

sp <- ggplot(snp_pass, aes(x=Aligner, y=METRIC.F1_Score, fill=CallerFilter)) + 
  geom_boxplot() + scale_y_continuous(limits=c(0.4, 1)) + theme_bw() + 
  theme(panel.grid=element_blank()) + guides(fill=F) + 
  geom_label(data=generate_label_df(0.45), aes(label=Label), 
             position=position_dodge(width=0.75), size=2.25) + 
  ylab('F1 Score (SNP)') +
  scale_fill_brewer(palette='Set3')

ia <- ggplot(indel_all, aes(x=Aligner, y=METRIC.F1_Score, fill=CallerFilter)) + 
  geom_boxplot() + scale_y_continuous(limits=c(0.4, 1)) + theme_bw() +
  theme(legend.position="bottom", panel.grid=element_blank()) + 
  geom_label(data=generate_label_df(0.45), aes(label=Label), 
             position=position_dodge(width=0.75), size=2.25) +
  scale_fill_brewer(palette='Set3')

ip <- ggplot(indel_pass, aes(x=Aligner, y=METRIC.F1_Score, fill=CallerFilter)) + 
  geom_boxplot() + scale_y_continuous(limits=c(0.4, 1)) + theme_bw() +
  theme(panel.grid=element_blank())  + guides(fill=F) + 
  geom_label(data=generate_label_df(0.45), aes(label=Label), 
             position=position_dodge(width=0.75), size=2.25) + 
  ylab('F1 Score (INDEL)') +
  scale_fill_brewer(palette='Set3')

print(sa)
print(ia)
f1 <- plot_grid(sp, ip, nrow=2)
print(f1)
plot_grid(sa, ia, nrow=2)

#######Pairwise comparison of tool performance#################
lattice.options(
  layout.heights=list(bottom.padding=list(x=0), top.padding=list(x=0)),
  layout.widths=list(left.padding=list(x=0), right.padding=list(x=0))
)

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

caller_pairwise = matrix(0, nrow=9, ncol=9)
rownames(caller_pairwise) = unique(snp_pass$CallerFilter)
colnames(caller_pairwise) = unique(snp_pass$CallerFilter)
pvals_callers = matrix("", nrow=9, ncol=9)
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

caller_pairwise_id = matrix(0, nrow=9, ncol=9)
rownames(caller_pairwise_id) = unique(snp_pass$CallerFilter)
colnames(caller_pairwise_id) = unique(snp_pass$CallerFilter)
pvals_callers_id = matrix("", nrow=9, ncol=9)
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

## Plotting Figure 2

bottom_left <- plot_grid(snp_aln, indel_aln, nrow=2, labels=c('b', 'c'))
bottom <- plot_grid(bottom_left, snp_call, indel_call, 
                    ncol=3, rel_widths = c(1, 1.05, 1), labels=c('', 'd', 'e'))
#bottom <- plot_grid(snp_call, indel_call, filtr_stats, nrow=1,
#                    rel_widths = c(1, 1, 0.95), labels=c('d', 'e', 'f'))
plot_grid(f1, bottom, nrow=2, rel_heights = c(1, 0.8), labels=c('a', ''))


########################################################################################

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

gain_loss$variable = factor(as.character(gain_loss$variable), levels=c('RecallLoss', 'PrecisionGain'))
gain_loss = gain_loss[gain_loss$Aligner != 'BT-E2E', ]
filtr_stats_box <- ggplot(gain_loss, aes(x=CallerFilter, y=abs(value), col=variable)) +
  geom_boxplot(position='dodge') + theme_bw() +
  theme(panel.grid = element_blank(), legend.position='bottom') +
  ylab('Change in metric\nvalue after filtering') +
  scale_color_simpsons() + scale_y_continuous(limits=c(0, 0.4)) +
  facet_grid(vars(Type), vars(ExperimentType))
print(filtr_stats_box)

mean_diff_pub = aggregate(value ~ variable + CallerFilter + Type + ExperimentType, gain_loss, median)
mean_diff_pub


################ Spread across aligners ########################

caller_perf = aggregate(METRIC.F1_Score ~ ExperimentType + Type + Sample + CallerFilter, pc_data, max)
caller_perf$min_f1 <- aggregate(METRIC.F1_Score ~ ExperimentType + Type + Sample + CallerFilter, 
                                pc_data, min)$METRIC.F1_Score
caller_perf$spread = caller_perf$METRIC.F1_Score - caller_perf$min_f1

wbwt <- ggplot(caller_perf, aes(x=CallerFilter, y=spread, fill=CallerFilter)) +
  geom_boxplot() + theme_bw() +
  scale_fill_brewer(palette = "Set3") +
  facet_wrap(~Type) + scale_y_continuous(limits=c(0, 0.3)) +
  ylab('Diff. between aligners\n(with Bowtie2)')
print(wbwt)

caller_perf_b = aggregate(METRIC.F1_Score ~ ExperimentType + Type + Sample + CallerFilter, 
                          pc_data[pc_data$Aligner != 'BT-E2E', ], max)
caller_perf_b$min_f1 <- aggregate(METRIC.F1_Score ~ ExperimentType + Type + Sample + CallerFilter, 
                          pc_data[pc_data$Aligner != 'BT-E2E', ], min)$METRIC.F1_Score
caller_perf_b$spread = caller_perf_b$METRIC.F1_Score - caller_perf_b$min_f1

wobwt <- ggplot(caller_perf_b, aes(x=CallerFilter, y=spread, fill=CallerFilter)) +
  geom_boxplot() + theme_bw() +
  scale_fill_brewer(palette = "Set3") +
  facet_wrap(~Type) + scale_y_continuous(limits=c(0, 0.3)) +
  ylab('Diff. between aligners\n(w/o Bowtie2)')
print(wobwt)

plot_grid(wbwt, wobwt, nrow=2)

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

stats_ordered = aggregated_stat2[order(aggregated_stat2$METRIC.F1_Score, decreasing = T), ]
print(stats_ordered[grepl('NO_', stats_ordered$Pipeline), ])


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
  scale_y_continuous(limits=c(0.7, 1)) +
  ylab('F1 Score (SNP)') + guides(color=F)


indel_gve <- ggplot(indel_pass[indel_pass$Aligner == 'NO', ], 
       aes(x=ExperimentType, y=METRIC.F1_Score, col=ExperimentType)) +
  geom_boxplot() + #geom_jitter() + 
  theme_bw() + theme(panel.grid=element_blank(), 
                     axis.text.x=element_blank(),
                     axis.title.x=element_blank()) +
  facet_wrap(~CallerFilter, nrow=1) + scale_color_npg() +
  scale_y_continuous(limits=c(0.7, 1)) +
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
  theme(panel.grid=element_blank()) + scale_color_brewer(palette = "Set3") +
  facet_wrap(~Type, nrow=1) + scale_y_continuous(limits=c(0, 1)) +
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
  theme(panel.grid=element_blank()) + scale_color_brewer(palette = "Set3") +
  facet_wrap(~Type, nrow=1) + scale_y_continuous(limits=c(0, 1)) +
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
  scale_color_brewer(palette = "Set3") +
  facet_grid(vars(Type), vars(ExperimentType)) + 
  scale_y_continuous(limits=c(0, 1)) +
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
  scale_color_brewer(palette = "Set3") +
  facet_grid(vars(Type), vars(ExperimentType)) + 
  scale_y_continuous(limits=c(0, 1)) +
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
  scale_color_brewer(palette = "Set3") +
  facet_wrap(~ExperimentType, nrow=2) + scale_y_continuous(limits=c(0, 1)) +
  scale_x_continuous(limits=c(-0.1, 0.9)) +
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

distances <- ggplot(dist_aggr[dist_aggr$Distance <= 125, ], aes(x=Distance, y=METRIC.F1_Score, 
                                ymin=METRIC.F1_Score-se, ymax=METRIC.F1_Score+se,
                                col=ExperimentType)) +
  geom_line(lwd=1) + #geom_errorbar(width=0.2, lwd=1.1) + 
  geom_point(size=2) + theme_bw() +
  theme(panel.grid=element_blank(), legend.position = "bottom") + 
  scale_color_npg() +
  facet_wrap(~Type, nrow=2) + scale_y_continuous(limits=c(0.75, 1)) +
  scale_x_continuous(limits=c(-1, 126)) +
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
  scale_x_continuous(limits=c(-0.1, 1))
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

## Bowtie2 local vs e2e

bowties <- bdata[bdata$Subset == "GRCh37_proteincoding_only.bed.gz" & 
        bdata$Subtype == "*", ]
bowties <- bowties[bowties$Aligner %in% c('BT-E2E', 'BT-LOC'), ]


bo_snp_pass <- bowties[bowties$Type == 'SNP' & bowties$Filter == 'PASS', ]
bo_indel_pass <- bowties[bowties$Type == 'INDEL' & bowties$Filter == 'PASS', ]

bo_sp <- ggplot(bo_snp_pass, aes(x=Aligner, y=METRIC.F1_Score, fill=CallerFilter)) + 
  geom_boxplot() + scale_y_continuous(limits=c(0.4, 1)) + theme_bw() + 
  theme(panel.grid=element_blank()) + guides(fill=F)  + 
  ylab('F1 Score (SNP)') +
  scale_fill_brewer(palette='Set3')
bo_sp

bo_ip <- ggplot(bo_indel_pass, aes(x=Aligner, y=METRIC.F1_Score, fill=CallerFilter)) + 
  geom_boxplot() + scale_y_continuous(limits=c(0.4, 1)) + theme_bw() +
  theme(panel.grid=element_blank())  + guides(fill=F) + 
  ylab('F1 Score (INDEL)') +
  scale_fill_brewer(palette='Set3')

bt_box <- plot_grid(bo_sp, bo_ip, nrow=2)
bt_box

aligner_pairwise_bt = matrix(0, nrow=2, ncol=2)
rownames(aligner_pairwise_bt) = unique(bo_snp_pass$Aligner)
colnames(aligner_pairwise_bt) = unique(bo_snp_pass$Aligner)
pvals_bt = matrix("", nrow=2, ncol=2)
rownames(pvals_bt) = unique(bo_snp_pass$Aligner)
colnames(pvals_bt) = unique(bo_snp_pass$Aligner)


for (i in unique(bo_snp_pass$Aligner)) {
  for (j in unique(bo_snp_pass$Aligner)) {
    left = as.numeric(bo_snp_pass[bo_snp_pass$Aligner == i, 'METRIC.F1_Score'])
    right = as.numeric(bo_snp_pass[bo_snp_pass$Aligner == j, 'METRIC.F1_Score'])
    difference = median(left - right)
    if (i == j){
      pvals_bt[i, j] = ""
    } else {
      pval = wilcox.test(left, right, paired=T)$p.value
      aligner_pairwise_bt[i, j] = difference
      pvals_bt[i, j] = ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", 
                                                       ifelse(pval < 0.05, "*", "")))
    }
  }
}

snp_bt <- drawCoolHM(aligner_pairwise_bt, pvals_bt)
print(snp_bt)


aligner_pairwise_bt_id = matrix(0, nrow=2, ncol=2)
rownames(aligner_pairwise_bt_id) = unique(bo_snp_pass$Aligner)
colnames(aligner_pairwise_bt_id) = unique(bo_snp_pass$Aligner)
pvals_bt_id = matrix("", nrow=2, ncol=2)
rownames(pvals_bt_id) = unique(bo_snp_pass$Aligner)
colnames(pvals_bt_id) = unique(bo_snp_pass$Aligner)


for (i in unique(bo_indel_pass$Aligner)) {
  for (j in unique(bo_indel_pass$Aligner)) {
    left = as.numeric(bo_indel_pass[bo_indel_pass$Aligner == i, 'METRIC.F1_Score'])
    right = as.numeric(bo_indel_pass[bo_indel_pass$Aligner == j, 'METRIC.F1_Score'])
    difference = median(left - right)
    if (i == j){
      pvals_bt_id[i, j] = ""
    } else {
      pval = wilcox.test(left, right, paired=T)$p.value
      aligner_pairwise_bt_id[i, j] = difference
      pvals_bt_id[i, j] = ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", 
                                                          ifelse(pval < 0.05, "*", "")))
    }
  }
}

indel_bt <- drawCoolHM(aligner_pairwise_bt_id, pvals_bt_id)
print(indel_bt)

bt_hm <- plot_grid(snp_bt, indel_bt, nrow=1, labels=c('b', 'c'))
bt_hm

plot_grid(bt_box, bt_hm, nrow=2, labels=c('a', ''), rel_heights = c(1, 0.6))



### Some general comparison of v3.3. and v. 4.2.


bdatav3 = read.table('./all_stats_by_strata_rework.tsv', sep='\t', header=T)

bdatav3$Precision.AM = (bdatav3$QUERY.TP + bdatav3$FP.gt)/(bdatav3$QUERY.TP + bdatav3$QUERY.FP)
bdatav3$Recall.AM = (bdatav3$QUERY.TP + bdatav3$FP.gt)/(bdatav3$QUERY.TP + bdatav3$FP.gt + bdatav3$TRUTH.FN)
bdatav3$F1_Score.AM = 2/(1/bdatav3$Precision.AM + 1/bdatav3$Recall.AM)
bdatav3$F1_Score.AM[bdatav3$METRIC.F1_Score == 0] = 0

bdatav3$Aligner = sapply(bdatav3$Aligner, 
                       function(x) aliases[aliases$src == as.character(x), 'alias'])
bdatav3$CallerFilter = sapply(bdatav3$CallerFilter, 
                            function(x) aliases[aliases$src == as.character(x), 'alias'])


pc_data_v4 <- bdata[bdata$Subset == "GRCh37_proteincoding_only.bed.gz" & 
                   bdata$Subtype == "*", ]
pc_data_v4$version = 'v4.2'

pc_data_v3 <- bdatav3[bdatav3$Subset == "GRCh37_proteincoding_only.bed.gz" & 
                      bdatav3$Subtype == "*", ]
pc_data_v3$version = 'v3.3'

pc_data_v4 <- pc_data_v4[, colnames(pc_data_v3)]
version_compare = rbind(pc_data_v3, pc_data_v4)
version_compare = version_compare[version_compare$Filter == 'PASS', ]

version_diff <- aggregate(METRIC.F1_Score~Sample+Type+ExperimentType+Aligner+CallerFilter,
                          version_compare, function(x) x[2] - x[1])
head(version_diff)

version_diff = version_diff[!(version_diff$Aligner %in% c('BT-E2E', 'BT-LOC')), ]

ggplot(version_diff, aes(x=Type, y=METRIC.F1_Score, fill=Type)) + 
  geom_violin(scale='width') + geom_boxplot(outlier.shape=NA, fill='white') +
  geom_hline(yintercept = 0, lty=2, lwd=1) +
  theme_bw() + theme(panel.grid = element_blank()) +
  ylab('Change in F1 value') + scale_fill_brewer(palette='Set1') +
  ylab('Performance difference\n(v4.2 - v3.3)') +
  facet_grid(vars(ExperimentType), vars(CallerFilter))


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

p_caller_pairwise = matrix(0, nrow=9, ncol=9)
rownames(p_caller_pairwise) = unique(snp_pass$CallerFilter)
colnames(p_caller_pairwise) = unique(snp_pass$CallerFilter)
p_pvals_callers = matrix("", nrow=9, ncol=9)
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

r_caller_pairwise = matrix(0, nrow=9, ncol=9)
rownames(r_caller_pairwise) = unique(snp_pass$CallerFilter)
colnames(r_caller_pairwise) = unique(snp_pass$CallerFilter)
r_pvals_callers = matrix("", nrow=9, ncol=9)
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

p_caller_pairwise_id = matrix(0, nrow=9, ncol=9)
rownames(p_caller_pairwise_id) = unique(snp_pass$CallerFilter)
colnames(p_caller_pairwise_id) = unique(snp_pass$CallerFilter)
p_pvals_callers_id = matrix("", nrow=9, ncol=9)
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

r_caller_pairwise_id = matrix(0, nrow=9, ncol=9)
rownames(r_caller_pairwise_id) = unique(snp_pass$CallerFilter)
colnames(r_caller_pairwise_id) = unique(snp_pass$CallerFilter)
r_pvals_callers_id = matrix("", nrow=9, ncol=9)
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
  scale_x_continuous(limits=c(-1, 126)) +
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
  scale_x_continuous(limits=c(-1, 126)) +
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
  scale_color_brewer(palette='Set3') +
  facet_grid(vars(Aligner), vars(ExperimentType)) + 
  scale_y_continuous(limits=c(0, 1)) +
  scale_x_continuous(limits=c(-0.1, 2.0)) +
  ylab('F1 Score') + guides(color=F)

cov_plot_indel_all <- ggplot(all_cov_stats_all[all_cov_stats_all$Type == 'INDEL', ],
                           aes(x=normcov, y=METRIC.F1_Score, 
                               ymin=METRIC.F1_Score-se, ymax=METRIC.F1_Score+se,
                               col=CallerFilter)) +
  geom_line(lwd=1) + #geom_errorbar(width=0.2, lwd=1.1) + 
  geom_point(size=2) + theme_bw() +
  theme(panel.grid=element_blank(), legend.position = "bottom") + 
  scale_color_brewer(palette='Set3') +
  facet_grid(vars(Aligner), vars(ExperimentType)) + 
  scale_y_continuous(limits=c(0, 1)) +
  scale_x_continuous(limits=c(-0.1, 2.0)) +
  ylab('F1 Score') + guides(color=F)

s7a <- plot_grid(cov_plot_snp_all, cov_plot_indel_all, ncol=2)
s7a



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
  scale_color_brewer(palette = "Set3") +
  facet_grid(vars(Aligner), vars(ExperimentType)) + 
  scale_y_continuous(limits=c(0, 1)) +
  scale_x_continuous(limits=c(0, 100)) +
  ylab('F1 Score') + guides(color=F)

gc_plot_indel_all <- ggplot(gc_aggr_all[gc_aggr_all$Type == 'INDEL',],
                          aes(x=GC_pct, y=METRIC.F1_Score, 
                              ymin=METRIC.F1_Score-se, ymax=METRIC.F1_Score+se,
                              col=CallerFilter)) +
  geom_line(lwd=1) + #geom_errorbar(width=0.2, lwd=1.1) + 
  geom_point(size=2) + theme_bw() +
  theme(panel.grid=element_blank(), legend.position = "bottom") + 
  scale_color_brewer(palette = "Set3") +
  facet_grid(vars(Aligner), vars(ExperimentType)) + 
  scale_y_continuous(limits=c(0, 1)) +
  scale_x_continuous(limits=c(0, 100)) +
  ylab('F1 Score') + guides(color=F)

s7b <- plot_grid(gc_plot_snp_all, gc_plot_indel_all, ncol=2)

s7ab <- plot_grid(s7a, s7b, nrow=2, labels=c('a', 'b'))
s7ab

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
  scale_color_brewer(palette = "Set3") +
  facet_grid(vars(Aligner), vars(ExperimentType)) + 
  scale_y_continuous(limits=c(0, 1)) +
  scale_x_continuous(limits=c(-0.1, 0.95)) +
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






### GIAB and non-GIAB concordance

# For Table 4

vc_ng <- read.table('./nongiab_v4_variant_count.tsv', sep='\t', header=F)
aggregate(V4~V1, vc_ng, median)


concord = read.table('./all_concordance.tsv', sep='\t', header=T)
head(concord)

concord$called = rowSums(concord[, c(5:8, 10:13)])
concord = concord[!(concord$HC_2DCNN == 1 & concord$called == 0), ]
#concord = concord[concord$giab == 'True' | concord$sample %in% c('Samp3', 
#        'wes_39.sample_Almaz02',  'wes_39.sample_Almaz05', 'wes_39.sample_Almaz07'), ]
concord$TRUTH = ifelse(is.na(concord$TRUTH), 'non-GIAB', ifelse(concord$TRUTH == 'YES',
                                                  'GIAB (true)', 'GIAB (false)'))
head(concord)

concord$TRUTH = factor(as.character(concord$TRUTH), levels=c('non-GIAB', 'GIAB (true)',
                                                             'GIAB (false)'))

caller_num <- ggplot(concord, aes(x=sample, y=called, fill=TRUTH)) + geom_violin(scale='width') +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  facet_grid(~TRUTH, scales="free_x", space='free_x') + theme_bw() +
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  scale_y_continuous(limits=c(0, 8)) +
  scale_fill_brewer(palette = "Accent") +
  ylab('Number of callers')
caller_num

giab_false <- concord[(concord$TRUTH == 'GIAB (true)' & concord$called < 8) | 
                        concord$TRUTH == 'GIAB (false)', ]
giab_false$ftype = ifelse(giab_false$TRUTH == 'GIAB (true)', 'FN', 'FP')

fn_fp <- ggplot(giab_false, aes(x=sample, y=called, fill=ftype)) + geom_violin(scale='width') +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  facet_wrap(~ftype, nrow=2) + theme_bw() +
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  scale_y_continuous(limits=c(0, 8)) + guides(fill=F) +
  scale_fill_brewer(palette = "Set1") +
  ylab('Number of callers')
fn_fp


aggregate(called~sample, concord, function(x) sum(x == 1))

cm <- melt(concord, id.vars=c('sample', 'giab', 'var_id', 'TRUTH', 'called'))
ps_agg <- aggregate(value~sample+variable, cm, function(x) sum(x > 0))

ps_agg

uq_calls <- melt(concord[concord$called == 1, c(1, 2, 5:8, 10:13)], id.vars=c('sample', 'giab'))
uq_call_df <- aggregate(value~variable+sample+giab, uq_calls, sum)
uq_call_df$variable = sapply(uq_call_df$variable, 
                             function(x) aliases[aliases$src == as.character(x), 'alias'])

uc <- ggplot(uq_call_df, aes(x=variable, y=value, fill=giab)) + 
  geom_boxplot(position=position_dodge(1)) +
  theme_bw() + scale_fill_brewer(palette="Accent") +
#  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  ylab('Unique calls')
uc

uq_noncalls <- melt(concord[concord$called == 7, c(1, 2, 5:8, 10:13)], id.vars=c('sample', 'giab'))
uq_noncall_df <- aggregate(value~variable+sample+giab, uq_noncalls, function(x) length(x) - sum(x))
uq_noncall_df$variable = sapply(uq_noncall_df$variable,
                             function(x) aliases[aliases$src == as.character(x), 'alias'])

unc <- ggplot(uq_noncall_df, aes(x=variable, y=value, fill=giab)) + 
  geom_boxplot(position=position_dodge(1)) +
  theme_bw() + scale_fill_brewer(palette="Accent") +
#  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  ylab('Unique non-calls')
unc

plot_grid(caller_num, plot_grid(uc, unc, nrow=1), nrow=2)


### Running some PCA analysis
giab_c = concord[concord$giab == 'True', c(1:8, 10:14)]
giab_c = giab_c[giab_c$called < 8, ]

head(giab_c)

pcaout_giab <- prcomp(t(as.matrix(giab_c[, 5:12])))
summary(pcaout_giab)
giab_pca_df <- as.data.frame(pcaout_giab$x)
giab_pca_df$caller = sapply(rownames(pcaout_giab$x),
                            function(x) aliases[aliases$src == as.character(x), 'alias'])

giab_pca <- ggplot(giab_pca_df, aes(x=PC1, y=PC2, col=caller, label=caller)) + geom_point(size=4) +
  theme_bw() + xlab('PC1 - 36.4%') + ylab('PC2 = 22.3%') + guides(color=F) +
  geom_label_repel()
giab_pca

# Sane on non-GIAB
nongiab_c = concord[concord$giab == 'False', c(1:8, 10:14)]
nongiab_c = nongiab_c[nongiab_c$called < 8, ]
#nongiab_c = nongiab_c[nongiab_c$sample %in% c('Samp3', 'wes_39.sample_Almaz02', 
#                                             'wes_39.sample_Almaz05', 'wes_39.sample_Almaz07'), ]

head(nongiab_c)

pcaout_nongiab <- prcomp(t(as.matrix(nongiab_c[, 5:12])))
summary(pcaout_nongiab)
nongiab_pca_df <- as.data.frame(pcaout_nongiab$x)
nongiab_pca_df$caller = sapply(rownames(pcaout_nongiab$x),
                            function(x) aliases[aliases$src == as.character(x), 'alias'])


nongiab_pca <- ggplot(nongiab_pca_df, aes(x=PC1, y=-PC2, col=caller, label=caller)) + geom_point(size=4) +
  theme_bw() + xlab('PC1 - 23.4%') + ylab('PC2 = 16.7%') + guides(color=F) +
  geom_label_repel()
nongiab_pca

# Layout

f5_top <- plot_grid(giab_pca, fn_fp, nrow=1, rel_widths = c(0.4, 0.6))
f5_top

f5_bottom_right <- plot_grid(uc, unc, nrow=2)
f5_bottom_right

f5_bottom <- plot_grid(nongiab_pca, f5_bottom_right, nrow=1, rel_widths = c(0.4, 0.6))
f5_bottom

plot_grid(f5_top, f5_bottom, nrow=2)
# 
