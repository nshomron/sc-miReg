
# Suppress warnings
options(warn = -1)

mirna <- read.table('GSE114071_NW_scsmRNA_K562_norm_log2.gct', header=T, stringsAsFactors=F)
mirna <- mirna[,c(1,3:21)]    # Column 2 containing descriptions is removed
head(mirna)

mrna_rpkm <- read.table('GSE114071_RPKM.tsv', header=T, stringsAsFactors=F)
# We retain 19 successfully profiled single cells, excluding K562_HalfCell_06
mrna_rpkm <- mrna_rpkm[,c(1:6,8:21)]
colnames(mrna_rpkm) <- colnames(mirna)
head(mrna_rpkm)

# TargetScan prediction
TS <- read.table('Predicted_Targets_Info.default_predictions.txt', header=T, sep='\t', stringsAsFactors=F)
family_info <- read.table('miR_Family_Info.txt', header=T, sep='\t', stringsAsFactors=F)

# Column 2~20 correspond to 19 single cells
index_cells <- 2:20

# We retained mRNAs expressed in no less than 5 single cells.
gene_filt <- which(apply(mrna_rpkm[,index_cells] > 0, 1, sum) >= 5)

meta <- data.frame(gene=mrna_rpkm[gene_filt, 'Name'])                        # Gene name
meta$mean_rpkm <- apply(mrna_rpkm[gene_filt, index_cells], 1, mean)          # Mean RPKM
meta$sd_rpkm <- apply(mrna_rpkm[gene_filt, index_cells], 1, sd)              # SD RPKM

# Linear regression model for log SD RPKM and log mean RPKM
fit <- lm(log10(sd_rpkm) ~ log10(mean_rpkm), meta)
summary(fit)

# Load required package for plotting and set plot size
library(ggplot2)
options(repr.plot.width=5, repr.plot.height=5)

ggplot(meta, aes(x=log10(mean_rpkm), y=log10(sd_rpkm))) + 
    geom_point(color='gray', size=0.1) +
    geom_smooth(fill='cyan', size=1) + 
    geom_smooth(method=lm, color='red', fill='pink', size=1) +
    labs(x=expression(log[10]~'Mean RPKM'), y=expression(log[10]~'SD RPKM')) +
    annotate('text', x=3, y=-1, label=expression(R^2==0.947), size=5) + 
    theme_bw()

# Calculate Residual SD
meta$rsd_rpkm <- fit$residuals

ggplot(meta, aes(x=log10(mean_rpkm), y=rsd_rpkm)) + 
    geom_point(color='gray', size=0.1) +
    geom_smooth(fill='cyan', size=1) + 
    geom_smooth(method=lm, color='red', fill='pink', size=1) +
    labs(x=expression(log[10]~'Mean RPKM'), y='Residual SD RPKM') +
    theme_bw()

# Get the target mRNAs' expression levels and noises of miRNAs
# mir: vector of string, miRNA names
# meta_var: string, a variable saved in dataframe meta. 

mir2meta_var <- function(mir, meta_var) 
{
    # Convert miRNAs into miRNA families
    mirna_family <- family_info[family_info$MiRBase.ID %in% mir, 'miR.family']
    
    # Predict target mRNAs of miRNA families
    mirna_target <- unique(TS[TS$miR.Family %in% mirna_family, 'Gene.Symbol'])
    
    # Get the expression levels or noises of target mRNAs calculated above
    target_meta_var <- meta[meta$gene %in% mirna_target, meta_var]
    
    return(target_meta_var)
}

# Calculate the mean expression levels of miRNAs
mirna_mean <- apply(mirna[,index_cells], 1, mean)

# Divide miRNAs into four groups according to their expression levels, and get expression levels (mean RPKM) in different groups
min_expression <- min(mirna[,-1])
level_background <- mir2meta_var(mirna[mirna_mean > min_expression & mirna_mean <= -12, 1], 'mean_rpkm')
level_LE_mirna <- mir2meta_var(mirna[mirna_mean > -12 & mirna_mean <= -9, 1], 'mean_rpkm')
level_ME_mirna <- mir2meta_var(mirna[mirna_mean > -9 & mirna_mean <= -6, 1], 'mean_rpkm')
level_HE_mirna <- mir2meta_var(mirna[mirna_mean > -6, 1], 'mean_rpkm')

# Divide miRNAs into four groups according to their expression levels, and get expression noises (RCV RPKM) in different groups
noise_background <- mir2meta_var(mirna[mirna_mean > min_expression & mirna_mean <= -12, 1], 'rsd_rpkm')
noise_LE_mirna <- mir2meta_var(mirna[mirna_mean > -12 & mirna_mean <= -9, 1], 'rsd_rpkm')
noise_ME_mirna <- mir2meta_var(mirna[mirna_mean > -9 & mirna_mean <= -6, 1], 'rsd_rpkm')
noise_HE_mirna <- mir2meta_var(mirna[mirna_mean > -6, 1], 'rsd_rpkm')

# In this dataframe, each row contains the name of a group, the number of miRNA and the number of target mRNA in the group. 
group_size <- data.frame(mirna_expression=c('background','LE','ME','HE'), 
                         mirna_number=c(sum(mirna_mean > min_expression & mirna_mean <= -12), sum(mirna_mean > -12 & mirna_mean <= -9), sum(mirna_mean > -9 & mirna_mean <= -6), sum(mirna_mean > -6)),
                         target_number=c(length(level_background), length(level_LE_mirna), length(level_ME_mirna), length(level_HE_mirna)))
group_size

# In this dataframe, each row contains the group of a target mRNA, its expression level and expression noise. 
meta_var_df <- data.frame(group=rep(group_size$mirna_expression, group_size$target_number),
                          level=c(level_background, level_LE_mirna, level_ME_mirna, level_HE_mirna),
                          noise=c(noise_background, noise_LE_mirna, noise_ME_mirna, noise_HE_mirna))

kruskal.test(level ~ group, data = meta_var_df) 
kruskal.test(noise ~ group, data = meta_var_df) 

# In this dataframe, each row contains names of two groups, statistical significance of difference between expression levels and noises from two groups. 
p_df <- data.frame(x=c('background', 'background', 'background', 'LE', 'LE', 'ME'),       # group A
                   y=c('LE', 'ME', 'HE', 'ME', 'HE', 'HE'),                               # group B
                   p_level=c(ks.test(level_background, level_LE_mirna)$p.value,                
                             ks.test(level_background, level_ME_mirna)$p.value,
                             ks.test(level_background, level_HE_mirna)$p.value,
                             ks.test(level_LE_mirna, level_ME_mirna)$p.value,
                             ks.test(level_LE_mirna, level_HE_mirna)$p.value,
                             ks.test(level_ME_mirna, level_HE_mirna)$p.value),
                   p_noise=c(ks.test(noise_background, noise_LE_mirna)$p.value,                
                             ks.test(noise_background, noise_ME_mirna)$p.value,
                             ks.test(noise_background, noise_HE_mirna)$p.value,
                             ks.test(noise_LE_mirna, noise_ME_mirna)$p.value,
                             ks.test(noise_LE_mirna, noise_HE_mirna)$p.value,
                             ks.test(noise_ME_mirna, noise_HE_mirna)$p.value))

library(ggpubr)
options(repr.plot.width=10, repr.plot.height=5)

plot_ecdf <- ggplot(meta_var_df, aes(x=log10(level), color=group)) + 
             stat_ecdf(geom = "step") + 
             lims(x=quantile(log10(meta_var_df$level),c(0.01,0.99))) + 
             labs(x=expression(log[10]~'Mean RPKM') ,y='Quantile', title="miRNA regulation on target mRNAs' expression levels") + 
             scale_color_discrete(breaks=c('HE','ME','LE','background'), 
                                  name='miRNA mean expression',
                                  labels=expression(log[2]~fraction > -6, -9 < log[2]~fraction <= -6, -12 < log[2]~fraction <= -9, 'All expressed')) + 
             theme_bw() + 
             theme(legend.position=c(0.7,0.2), title=element_text(size=10))
    
plot_p <- ggplot(p_df,aes(x=x,y=y)) + 
          geom_raster(aes(fill = -log10(p_level))) + 
          geom_text(aes(label = round(-log10(p_level),2)), size=7) + 
          scale_fill_gradient(low = "white", high = "red") + 
          scale_x_discrete(limits=c('background','LE','ME'), labels=expression('All expressed', -12 < log[2]~fraction <= -9, -9 < log[2]~fraction <= -6)) + 
          scale_y_discrete(limits=c('LE','ME','HE'), labels=expression(-12 < log[2]~fraction <= -9, -9 < log[2]~fraction <= -6, log[2]~' fraction > -6')) +
          labs(x='',y='',title=expression(-log[10]~'P value (KS-test)')) + 
          theme(legend.position = 'none', axis.text = element_text(angle=30, hjust=1))
          
ggarrange(plot_ecdf, plot_p)

plot_ecdf <- ggplot(meta_var_df, aes(x=noise, color=group)) + 
             stat_ecdf(geom = "step") + 
             lims(x=quantile(meta_var_df$noise,c(0.01,0.99))) + 
             labs(x='Residual SD RPKM' ,y='Quantile', title="miRNA regulation on target mRNAs' expression noises") + 
             scale_color_discrete(breaks=c('HE','ME','LE','background'), 
                                  name='miRNA mean expression',
                                  labels=expression(log[2]~fraction > -6, -9 < log[2]~fraction <= -6, -12 < log[2]~fraction <= -9, 'All expressed')) +
             theme_bw() + 
             theme(legend.position=c(0.7,0.2), title=element_text(size=10))

plot_p <- ggplot(p_df,aes(x=x,y=y)) + 
          geom_raster(aes(fill = -log10(p_noise))) + 
          geom_text(aes(label = round(-log10(p_noise),2)), size=7) + 
          scale_fill_gradient(low = "white", high = "red") + 
          scale_x_discrete(limits=c('background','LE','ME'), labels=expression('All expressed', -12 < log[2]~fraction <= -9, -9 < log[2]~fraction <= -6)) + 
          scale_y_discrete(limits=c('LE','ME','HE'), labels=expression(-12 < log[2]~fraction <= -9, -9 < log[2]~fraction <= -6, log[2]~fraction > -6)) +
          labs(x='',y='',title=expression(-log[10]~'P value (KS-test)')) + 
          theme(legend.position = 'none', axis.text = element_text(angle=30, hjust=1))
          
ggarrange(plot_ecdf, plot_p)

# Read DCA output
mrna_dca <- read.table('GSE114071_DCA.tsv', header=T, stringsAsFactors=F, row.names = NULL)
mrna_dca <- mrna_dca[,c(1:6,8:21)]
colnames(mrna_dca) <- colnames(mirna)
head(mrna_dca)

# Read matrix of effective gene length
mrna_length <- read.table('GSE114071_length.tsv', header=T, stringsAsFactors=F)
mrna_length <- mrna_length[,c(1:6,8:21)]
colnames(mrna_length) <- colnames(mirna)

# Normalize DCA output by dividing gene length
mrna_dca[,-1] <- mrna_dca[,-1] / mrna_length[mrna_length$Name %in% mrna_dca$Name, -1]

# Only retain the matched genes
mrna_dca <- mrna_dca[mrna_dca$Name %in% meta$gene,]
matched_rows <- which(meta$gene %in% mrna_dca$Name)

meta[matched_rows, 'mean_dca'] <- apply(mrna_dca[,index_cells], 1, mean)      # Mean DCA normalized counts
meta[matched_rows, 'sd_dca'] <- apply(mrna_dca[,index_cells], 1, sd)          # SD DCA normalized counts

# Calculate Residual SD
fit <- lm(log10(sd_dca) ~ log10(mean_dca), meta)
meta[matched_rows, 'rsd_dca'] <-  fit$residuals
summary(fit)

plot_left <- ggplot(meta, aes(x=log10(mean_dca), y=log10(sd_dca))) + 
             geom_point(color='gray', size=0.1) +
             geom_smooth(fill='cyan', size=1) + 
             geom_smooth(method=lm, color='red', fill='pink', size=1) +
             labs(x=expression(log[10]~'Mean DCA normalized count'), y=expression(log[10]~'SD DCA normalized count')) +
             annotate('text', x=0, y=-3.5, label=expression(R^2==0.968), size=5) + 
             theme_bw()

plot_right <- ggplot(meta, aes(x=log10(mean_rpkm), y=rsd_rpkm)) + 
              geom_point(color='gray', size=0.1) + theme_bw() +
              geom_smooth(fill='cyan', size=1) + 
              geom_smooth(method=lm, color='red', fill='pink', size=1) +
              labs(x=expression(log[10]~'Mean DCA normalized count'), y='Residual SD DCA normalized count') 

ggarrange(plot_left, plot_right)

# Get expression levels (mean DCA normalized counts) in different groups
noise_background <- mir2meta_var(mirna[mirna_mean > min_expression & mirna_mean <= -12, 1], 'mean_dca')
noise_LE_mirna <- mir2meta_var(mirna[mirna_mean > -12 & mirna_mean <= -9, 1], 'mean_dca')
noise_ME_mirna <- mir2meta_var(mirna[mirna_mean > -9 & mirna_mean <= -6, 1], 'mean_dca')
noise_HE_mirna <- mir2meta_var(mirna[mirna_mean > -6, 1], 'mean_dca')

# Get expression noises (RCV DCA normalized counts) in different groups
noise_background <- mir2meta_var(mirna[mirna_mean > min_expression & mirna_mean <= -12, 1], 'rsd_dca')
noise_LE_mirna <- mir2meta_var(mirna[mirna_mean > -12 & mirna_mean <= -9, 1], 'rsd_dca')
noise_ME_mirna <- mir2meta_var(mirna[mirna_mean > -9 & mirna_mean <= -6, 1], 'rsd_dca')
noise_HE_mirna <- mir2meta_var(mirna[mirna_mean > -6, 1], 'rsd_dca')

# The number of miRNAs and target mRNAs in each group.
group_size <- data.frame(mirna_expression=c('background','LE','ME','HE'), 
              mirna_number=c(sum(mirna_mean > min_expression & mirna_mean <= -12), sum(mirna_mean > -12 & mirna_mean <= -9), sum(mirna_mean > -9 & mirna_mean <= -6), sum(mirna_mean > -6)),
              target_number=c(length(noise_background), length(noise_LE_mirna), length(noise_ME_mirna), length(noise_HE_mirna)))
              
# In this dataframe, each row contains the group of a target mRNA, its expression level and expression noise. 
meta_var_df <- data.frame(group=rep(group_size$mirna_expression, group_size$target_number),
                          level=c(level_background, level_LE_mirna, level_ME_mirna, level_HE_mirna),
                          noise=c(noise_background, noise_LE_mirna, noise_ME_mirna, noise_HE_mirna))
                       
p_df <- data.frame(x=c('background','background','background','LE','LE','ME'),
                   y=c('LE','ME','HE','ME','HE','HE'),
                   p_level=c(ks.test(level_background, level_LE_mirna)$p.value,
                             ks.test(level_background, level_ME_mirna)$p.value,
                             ks.test(level_background, level_HE_mirna)$p.value,
                             ks.test(level_LE_mirna, level_ME_mirna)$p.value,
                             ks.test(level_LE_mirna, level_HE_mirna)$p.value,
                             ks.test(level_ME_mirna, level_HE_mirna)$p.value),
                   p_noise=c(ks.test(noise_background, noise_LE_mirna)$p.value,
                             ks.test(noise_background, noise_ME_mirna)$p.value,
                             ks.test(noise_background, noise_HE_mirna)$p.value,
                             ks.test(noise_LE_mirna, noise_ME_mirna)$p.value,
                             ks.test(noise_LE_mirna, noise_HE_mirna)$p.value,
                             ks.test(noise_ME_mirna, noise_HE_mirna)$p.value))

kruskal.test(level ~ group, data = meta_var_df) 
kruskal.test(noise ~ group, data = meta_var_df) 

plot_ecdf <- ggplot(meta_var_df, aes(x=log10(level), color=group)) + 
             stat_ecdf(geom = "step") + 
             lims(x=quantile(log10(meta_var_df$level), c(0.01,0.99))) + 
             labs(x=expression(log[10]~'Mean DCA normalized count') ,y='Quantile', title="miRNA regulation on target mRNAs' expression levels") + 
             scale_color_discrete(breaks=c('HE','ME','LE','background'), 
                                  name='miRNA Mean Expression',
                                  labels=expression(log[2]~fraction > -6, -9 < log[2]~fraction <= -6, -12 < log[2]~fraction <= -9, 'All expressed')) + 
             theme_bw() + 
             theme(legend.position=c(0.7,0.2), title=element_text(size=10))
    
plot_p <- ggplot(p_df,aes(x=x,y=y)) + 
          geom_raster(aes(fill = -log10(p_level))) + 
          geom_text(aes(label = round(-log10(p_level),2)), size=7) + 
          scale_fill_gradient(low = "white", high = "red") + 
          scale_x_discrete(limits=c('background','LE','ME'), labels=expression('All expressed', -12 < log[2]~fraction <= -9, -9 < log[2]~fraction <= -6)) + 
          scale_y_discrete(limits=c('LE','ME','HE'), labels=expression(-12 < log[2]~fraction <= -9, -9 < log[2]~fraction <= -6, log[2]~fraction > -6)) +
          labs(x='',y='',title=expression(-log[10]~'P value (KS-test)')) + 
          theme(legend.position = 'none', axis.text = element_text(angle=30, hjust=1))
          
ggarrange(plot_ecdf, plot_p)

plot_ecdf <- ggplot(meta_var_df, aes(x=noise, color=group)) + 
             stat_ecdf(geom = "step") + 
             lims(x=quantile(meta_var_df$noise, c(0.01,0.99))) + 
             labs(x='Residual SD DCA normalized count' ,y='Quantile', title="miRNA regulation on target mRNAs' expression noises") + 
             scale_color_discrete(breaks=c('HE','ME','LE','background'), 
                                  name='miRNA Mean Expression',
                                  labels=expression(log[2]~fraction > -6, -9 < log[2]~fraction <= -6, -12 < log[2]~fraction <= -9, 'All expressed')) + 
             theme_bw() + 
             theme(legend.position=c(0.7,0.2), title=element_text(size=10))
    
plot_p <- ggplot(p_df,aes(x=x,y=y)) + 
          geom_raster(aes(fill = -log10(p_noise))) + 
          geom_text(aes(label = round(-log10(p_noise),2)), size=7) + 
          scale_fill_gradient(low = "white", high = "red") + 
          scale_x_discrete(limits=c('background','LE','ME'), labels=expression('All expressed', -12 < log[2]~fraction <= -9, -9 < log[2]~fraction <= -6)) + 
          scale_y_discrete(limits=c('LE','ME','HE'), labels=expression(-12 < log[2]~fraction <= -9, -9 < log[2]~fraction <= -6, log[2]~fraction > -6)) +
          labs(x='',y='',title=expression(-log[10]~'P value (KS-test)')) + 
          theme(legend.position = 'none', axis.text = element_text(angle=30, hjust=1))
          
ggarrange(plot_ecdf, plot_p)

sessionInfo()
