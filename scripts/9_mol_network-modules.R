# ANALYSES ON METABARCODING DATA - COOCCURRENCES NETWORK + MODULES STUDY

remove(list = ls(all.names = TRUE))
set.seed(02121999)

##### LIBRAIRIES
library(tidyverse)
library(vegan)
library(lubridate)
library(igraph)
library(Hmisc)
library(funrar)
library(corrplot)
library(cowplot)
library(rstatix)
library(RColorBrewer)
library(ggpubr)

##### FUNCTION
source("0_network_and_save_relative.R")

##### COLOUR VECTOR
colour <- c('#CDAD00','#458B00','#21918c') # Clusters

################################################################################ DATA IMPORTATION

ASV <- read.csv2("ASV_0-05.csv", header = TRUE, sep = ",")
ASVtot <- read.csv2("ASV_rarefy.csv", header = TRUE, sep = ",")
metadata <- read.csv2("metadata_sample.csv", header = TRUE, sep = ",")
metadata <- metadata[-15,] # Remove 329-Day (October-2016)
clusters <- read.csv2("clusters_pcoa.csv", header = TRUE, sep = ",")
morpho <- read.csv2("img_samplingcoord.csv", header = TRUE, sep = ",")
morpho <- morpho[-25,] # Remove October-Day-2016
ncbi <- read.csv2("results_byhand.csv", header = TRUE, sep = ";")

################################################################################ NETWORK ####

## Thresholds values to filter the significant edges
rho <- 0.6 # Spearman's correlation coefficient threshold
pvalue <- 0.01 # P-value threshold on rho

## Network computation
cn <- network_and_save_relative(ASV, rho, pvalue)

g_cn <- cn$graph
topo_network <- cn$network
write.table(topo_network, file = "mol_toponetwork.csv", sep = ",", quote = FALSE, row.names = F, col.names = T) # Save the dataframe
topo_node <- cn$node

## Modularity 
m_cn <- cluster_leading_eigen(g_cn)

length(m_cn) # 69 modules
sizes(m_cn) # Most of the time 1 ASV per module
modularity(m_cn) # Modularity 0.716

## Remove the modules with only 1 ASV
Small <- which(table(m_cn$membership) < 2)
Keep <- V(g_cn)[!(m_cn$membership %in% Small)]

g_cn2 <- induced_subgraph(g_cn, Keep) # Subgraph with only modules containing > 1 ASV
m_cn2 <- cluster_leading_eigen(g_cn2) # Compute the modules

## Stock in a database
modules <- ASV[Keep,1:28]
modules$Membership <- m_cn2$membership

modules_tax <- left_join(modules, ASVtot[,c(1,29:49)], by = 'Representative_Sequence') # Add the taxonomy from 'rarefy' file
write.table(modules_tax, file = "ASV_Modules.csv", sep = ",", quote = FALSE, row.names = F) # Save dataframe

## Merge with BLAST/NCBI taxonomy
modules_tax$Membership <- as.character(modules_tax$Membership)
mod_tax <- modules_tax %>% filter(Membership %in% c('1','7','8','9','10','15'))
colnames(ncbi)[1] <- 'Representative_Sequence'
mod_ncbi <- left_join(mod_tax, ncbi[,-25], by = 'Representative_Sequence')
mod_ncbi <- unique(mod_ncbi)

mod_ncbi %>% group_by(Membership) %>% summarise(n = n())
mod_tax %>% group_by(Membership) %>% summarise(n = n())
which(mod_ncbi$Representative_Sequence == 'M04837_139_000000000-C4LDR_1_1101_19779_4115')

mod_ncbi <- mod_ncbi[-142,]
which(mod_ncbi$Representative_Sequence == 'M04837_139_000000000-C4LDR_1_1101_19779_4115') # Only one, it works
write.table(mod_ncbi, file = "ASV_ModTax.csv", sep = ",", quote = FALSE, row.names = F) # Save dataframe (NCBI taxo)

################################################################################ 15 MODULES RELATIVE ABUNDANCES ####

## Get the RA% from all the modules
mod_df <- as.data.frame(matrix(ncol = 15, nrow = 27))
mod_ac <- as.data.frame(matrix(ncol = 15, nrow = 27)) # to get absolute counts of reads
mod <- as.data.frame(matrix(ncol = 2, nrow = 15))

for (i in 1:length(unique(m_cn2$membership))) {
  ASV_mod <- mod_ncbi$Representative_Sequence[mod_ncbi$Membership == i]
  mod[i,1] <- i # Module number 
  mod[i,2] <- length(ASV_mod) # Number of ASVs per module
  ASV_modn <- which(ASV$Representative_Sequence %in% ASV_mod)
  for (j in 1:27) {
    ASVsum <- sum(ASV[ASV_modn,j+1])
    mod_ac[j,i] <- ASVsum # Absolute counts of reads per module and sampling
    ASVrel <- 100 * ASVsum / sum(ASV[,j+1]) # RA% of each module (= sum of its ASVs) among the total 225 ASVs
    mod_df[j,i] <- ASVrel
  }
}

colnames(mod) <- c('Module_number','Module_length')
write.table(t(mod), file = "ASV_ModLength.csv", sep = ",", quote = FALSE, row.names = T, col.names = F) # Save dataframe (modules with number of ASVs)

## Add the samplings as the first column
mod_df <- cbind(clusters$Sampling, mod_df)
colnames(mod_df) <- c('Samplings', paste0('Mod', as.character(seq(1,15)))) # Rename all the columns
rownames(mod_df) <- NULL

################################################################################ PLOT - Relative abundances throughout time ####

## Get RA per module throughout time
RA_long <- pivot_longer(mod_df, cols = starts_with('Mod'))
colnames(RA_long)[2:3] <- c('Modules','RA')
RA_long$Time <- substr(RA_long$Samplings, start = 5, stop = 5)
RA_long$Samplings <- substr(RA_long$Samplings, start = 1, stop = 4)

# Create fake rows to represent Mar-2017 (no sampling)
DMar17A <- c('B336','Mod1',0,'D')
DMar17B <- c('B336','Mod7',0,'D')
DMar17C <- c('B336','Mod8',0,'D')
DMar17D <- c('B336','Mod9',0,'D')
DMar17E <- c('B336','Mod10',0,'D')
DMar17F <- c('B336','Mod15',0,'D')
RA_long <- rbind(RA_long, DMar17A, DMar17B, DMar17C,
                 DMar17D, DMar17E, DMar17F)

RA_long$Months <- factor(RA_long$Samplings, 
                         labels =  c('Mar-16','Apr-16','May-16','Jun-16','Jul-16',
                                     'Aug-16','Sep-16','Oct-16','Nov-16','Dec-16','Jan-17','Feb-17',
                                     'Apr-17','May-17','Mar-17'))
RA_long$Months <- ordered(RA_long$Months)
RA_long$Months <- ordered(RA_long$Months, levels = c('Mar-16','Apr-16','May-16','Jun-16','Jul-16',
                                                     'Aug-16','Sep-16','Oct-16','Nov-16','Dec-16','Jan-17','Feb-17',
                                                     'Mar-17','Apr-17','May-17'))

RA_long$Modules <- substr(RA_long$Modules, start = 4, stop = 6)

## Get the RA% for ASVs not in any module
RA_not <- as.data.frame(matrix(ncol = 5, nrow = 27))
colnames(RA_not) <- c('Samplings','Modules','RA','Time','Months')
RA_not$Samplings <- mod_df$Samplings
RA_not$Time <- substr(RA_not$Samplings, start = 5, stop = 5)
RA_not$Time <- factor(RA_not$Time, levels = c('D','N'), labels = c('Day','Night'))
RA_not$Samplings <- substr(RA_not$Samplings, start = 1, stop = 4)
RA_not$Modules <- 'O'
RA_not$RA <- 100 - rowSums(mod_df[,-1])
RA_not$Months <- factor(RA_not$Samplings, labels = c('Mar-16','Apr-16','May-16','Jun-16','Jul-16',
                                                     'Aug-16','Sep-16','Oct-16','Nov-16','Dec-16','Jan-17','Feb-17',
                                                     'Apr-17','May-17'))
RA_not$Months <- ordered(RA_not$Months)
RA_not$Months <- ordered(RA_not$Months, levels = c('Mar-16','Apr-16','May-16','Jun-16','Jul-16',
                                                   'Aug-16','Sep-16','Oct-16','Nov-16','Dec-16','Jan-17','Feb-17',
                                                   'Mar-17','Apr-17','May-17'))

# Filter to keep only modules with at least 10 ASVs
RA_6mod <- RA_long %>% filter(Modules %in% c('1','7','8','9','10','15'))
RA_6mod <- as.data.frame(RA_6mod)
RA_6mod$Modules <- factor(RA_6mod$Modules, levels = c('1','7','8','9','10','15'),
                          labels = LETTERS[1:6])
RA_6mod$Time <- factor(RA_6mod$Time, levels = c('D','N'), labels = c('Day','Night'))
RA_6mod$RA <- as.numeric(RA_6mod$RA)

# Filter to sum all other modules
RA_other <- RA_long %>% filter(Modules %in% c('2','3','4','5','6','11','12','13','14'))
RA_other <- as.data.frame(RA_other)
NMar17 <- c('B336','X',0,'N')
RA_other <- rbind(RA_other, NMar17) # Warning message, but no worries, everything is fine
RA_other$Time <- factor(RA_other$Time, levels = c('D','N'), labels = c('Day','Night'))
RA_other <- RA_other %>% group_by(Samplings, Time) %>%
  mutate(RA = sum(as.numeric(RA))) # Get the sum for all modules (<10ASVs) per cruise and time
RA_other$Modules <- 'X'
RA_other <- unique(RA_other) # Keep one line per cruise and time

RA_other$Months <- factor(RA_other$Samplings, 
                          labels =  c('Mar-16','Apr-16','May-16','Jun-16','Jul-16',
                                      'Aug-16','Sep-16','Oct-16','Nov-16','Dec-16','Jan-17','Feb-17',
                                      'Apr-17','May-17','Mar-17'))
RA_other$Months <- ordered(RA_other$Months)
RA_other$Months <- ordered(RA_other$Months, levels = c('Mar-16','Apr-16','May-16','Jun-16','Jul-16',
                                                       'Aug-16','Sep-16','Oct-16','Nov-16','Dec-16','Jan-17','Feb-17',
                                                       'Mar-17','Apr-17','May-17'))


# Merge all tables (mod>10ASVs, mod<10ASVs = X, not in mod = O)
RA_modoth <- bind_rows(RA_6mod, RA_other)
RA_all <- bind_rows(RA_modoth, RA_not) 
RA_all <- arrange(RA_all, Samplings, Time) # Re-order by cruise and time

remove(RA_6mod, RA_long, RA_modoth, RA_not, RA_other)

# test <- RA_all %>% group_by(Samplings, Time) %>%
#   mutate(test = sum(RA)) # it worked, all 100% at each cruise*time

# Rename letters for the legend
RA_all$Legend <- factor(RA_all$Modules, levels = c('A','B','C','D','E','F','O','X'),
                        labels = c('Module A', 'Module B', 'Module C',
                                   'Module D', 'Module E', 'Module F',
                                   'Other modules', 'Not in module'))

## Plot - Barplots at 100%
axiscol <- c('#21918c','#21918c','#458b00','#CDAD00','#CDAD00',
             '#CDAD00','#CDAD00','#CDAD00','#21918c','#21918c',
             '#21918c','#21918c','dimgrey','#458b00','#21918c') # Colours matching with PCoA clusters

## Facet grid plot
RA_all$RA <- ifelse(RA_all$Time == 'Night', -RA_all$RA, RA_all$RA) # Make negative values for Night samples

# Individual modules
modA <- ggplot(RA_all %>% filter(Legend == 'Module A')) +
  geom_bar(aes(x = Months, y = RA, fill = Time),
           size = 1, color = 'black', stat = 'identity') +
  geom_hline(aes(yintercept = 0), colour = 'black', size = 1) +
  theme_classic(base_size = 30) +
  xlab('Samplings') +
  scale_y_continuous(name = 'Relative counts of reads (%)',
                     limits = c(-40,40), 
                     breaks = seq(-40,40,10), 
                     labels = c(40,30,20,10,0,10,20,30,40)) +
  scale_fill_manual(values = c('white','black')) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = c(0.92,0.94),
        legend.title = element_blank(),
        legend.direction = 'vertical',
        legend.text = element_text(size = 25),
        legend.spacing.y = unit(0, 'cm'),
        legend.margin = margin(t = 0, unit = 'cm'))

modB <-  ggplot(RA_all %>% filter(Legend == 'Module B')) +
  geom_bar(aes(x = Months, y = RA, fill = Time),
           size = 1, color = 'black', stat = 'identity') +
  geom_hline(aes(yintercept = 0), colour = 'black', size = 1) +
  theme_classic(base_size = 30) +
  xlab('Samplings') +
  scale_y_continuous(name = 'Relative counts of reads',
                     limits = c(-30,30), 
                     breaks = seq(-30,30,10), 
                     labels = c(30,20,10,0,10,20,30)) +
  scale_fill_manual(values = c('white','black')) +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        legend.position = 'none',
        legend.direction = 'horizontal',
        legend.text = element_text(size = 25),
        legend.spacing.y = unit(0, 'cm'),
        legend.margin = margin(t = 0, unit = 'cm'))

modC <-  ggplot(RA_all %>% filter(Legend == 'Module C')) +
  geom_bar(aes(x = Months, y = RA, fill = Time),
           size = 1, color = 'black', stat = 'identity') +
  geom_hline(aes(yintercept = 0), colour = 'black', size = 1) +
  theme_classic(base_size = 30) +
  xlab('Samplings') +
  scale_y_continuous(name = 'Relative counts of reads',
                     limits = c(-50,50), 
                     breaks = seq(-50,50,10), 
                     labels = c(50,40,30,20,10,0,10,20,30,40,50)) +
  scale_fill_manual(values = c('white','black')) +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        legend.position = 'none',
        legend.direction = 'horizontal',
        legend.text = element_text(size = 25),
        legend.spacing.y = unit(0, 'cm'),
        legend.margin = margin(t = 0, unit = 'cm'))

modD <-  ggplot(RA_all %>% filter(Legend == 'Module D')) +
  geom_bar(aes(x = Months, y = RA, fill = Time),
           size = 1, color = 'black', stat = 'identity') +
  geom_hline(aes(yintercept = 0), colour = 'black', size = 1) +
  theme_classic(base_size = 30) +
  xlab('Samplings') +
  scale_y_continuous(name = 'Relative counts of reads (%)',
                     limits = c(-60,60), 
                     breaks = seq(-60,60,10), 
                     labels = c(60,50,40,30,20,10,0,10,20,30,40,50,60)) +
  scale_fill_manual(values = c('white','black')) +
  theme(axis.text.x = element_text(colour = axiscol, angle = 90, vjust = 0.5),
        axis.title.x = element_blank(),
        legend.position = 'none',
        legend.direction = 'horizontal',
        legend.text = element_text(size = 25),
        legend.spacing.y = unit(0, 'cm'),
        legend.margin = margin(t = 0, unit = 'cm')) # Warning message, but no worries, everything is fine

modE <-  ggplot(RA_all %>% filter(Legend == 'Module E')) +
  geom_bar(aes(x = Months, y = RA, fill = Time),
           size = 1, color = 'black', stat = 'identity') +
  geom_hline(aes(yintercept = 0), colour = 'black', size = 1) +
  theme_classic(base_size = 30) +
  xlab('Samplings') +
  scale_y_continuous(name = 'Relative counts of reads',
                     limits = c(-40,40), 
                     breaks = seq(-40,40,10), 
                     labels = c(40,30,20,10,0,10,20,30,40)) +
  scale_fill_manual(values = c('white','black')) +
  theme(axis.text.x = element_text(colour = axiscol, angle = 90, vjust = 0.5),
        axis.title = element_blank(),
        legend.position = 'none',
        legend.direction = 'horizontal',
        legend.text = element_text(size = 25),
        legend.spacing.y = unit(0, 'cm'),
        legend.margin = margin(t = 0, unit = 'cm')) # Warning message, but no worries, everything is fine

modF <-  ggplot(RA_all %>% filter(Legend == 'Module F')) +
  geom_bar(aes(x = Months, y = RA, fill = Time),
           size = 1, color = 'black', stat = 'identity') +
  geom_hline(aes(yintercept = 0), colour = 'black', size = 1) +
  theme_classic(base_size = 30) +
  xlab('Samplings') +
  scale_y_continuous(name = 'Relative counts of reads',
                     limits = c(-10,10), 
                     breaks = seq(-10,10,5), 
                     labels = c(10,5,0,5,10)) +
  scale_fill_manual(values = c('white','black')) +
  theme(axis.text.x = element_text(colour = axiscol, angle = 90, vjust = 0.5),
        axis.title = element_blank(),
        legend.position = 'none',
        legend.direction = 'horizontal',
        legend.text = element_text(size = 25),
        legend.spacing.y = unit(0, 'cm'),
        legend.margin = margin(t = 0, unit = 'cm')) # Warning message, but no worries, everything is fine

# All together
plot_grid(modA, modB, modC,
          NULL, NULL, NULL,
          modD, modE, modF,
          ncol = 3, rel_heights = c(1,-0.15,1), 
          align = 'hv', labels = c('A','B','C','','','','D','E','F'),
          label_size = 25, label_x = 0.2) # Figure 5

remove(modA,modB,modC,modD,modE,modF)

################################################################################ 6 MODULES STUDY - Correlations (morpho + export) ####

mod_df <- mod_df[,c('Samplings','Mod1','Mod7','Mod8','Mod9','Mod10','Mod15')]
n_mod_red <- 6

mod_df <- left_join(mod_df, morpho[,-1], by = c('Samplings' = 'cruisename')) # Add the morphological values (coordinates of the centroids of each sample date along each axis)
colnames(mod_df)[8:11] <- paste0('Morpho', as.character(1:4))

mod_df$C200 <- metadata$C200 # Add export values

## Correlation matrix
mod_num <- as.data.frame(lapply(mod_df[,-1], function(x) if(is.character(x)) as.numeric(x) else x))
rownames(mod_num) <- mod_df[,1]

mod_cor <- cor(mod_num, method = 'pearson') # Pearson's correlation coefficients
p_mat <- cor.mtest(mat = mod_num)$p # Associated p-values

corrplot(mod_cor, type = "upper",
         tl.col = "black", tl.srt = 45,
         addCoef.col = "black",
         method = 'color',
         p.mat = p_mat, sig.level = 0.01, insig = "blank",
         diag = FALSE)

remove(p_mat, mod_cor, mod_num)

################################################################################ 6 MODULES STUDY - Day/Night effect ####

mod_df[,8:12] <- as.data.frame(lapply(mod_df[,8:12], function(x) if(is.character(x)) as.numeric(x) else x))
mod_df$Time <- substr(mod_df$Samplings, start = 5, stop = 5)

test_sign <- as.data.frame(matrix(nrow = n_mod_red, ncol = 4))
colnames(test_sign) <- c('Variables', 'Shapiro_P', 'Levene_P', 'Test')

for (i in 1:n_mod_red) {
  test_sign$Variables[i] <- colnames(mod_df)[i+1]
  
  test_sign$Shapiro_P[i] <- shapiro.test(mod_df[,i+1])$p
  test_sign$Levene_P[i] <- levene_test(data = mod_df, formula = mod_df[,i+1] ~ Time)$p
  
  test_sign$Test[i] <- ifelse(test_sign$Shapiro_P[i] > 0.05 & test_sign$Levene_P[i] > 0.05, 'Ttest', 'Wilc')
} # 2 ok for t-test / 4 ok for wilcoxon

wilcoxon <- as.data.frame(matrix(nrow = 11, ncol = 2)) # Easier to do non parametric tests for all of them
colnames(wilcoxon) <- c('Variables','Pvalue')

for (i in 1:11) {
  test <- wilcox.test(mod_df[,i+1] ~ mod_df$Time)
  wilcoxon[i,1] <- colnames(mod_df)[i+1]
  wilcoxon[i,2] <- test$p.value
}  # Warning message, but no worries, everything is fine

daynight <- wilcoxon %>% filter(Pvalue < 0.05) # Only module A/1

################################################################################ 6 MODULES STUDY - PCoA pattern ####

mod_df$PCoA_Clusters <- clusters$AsvClusters

test_sign <- as.data.frame(matrix(nrow = 11, ncol = 4))
colnames(test_sign) <- c('Variables', 'Shapiro_P', 'Levene_P', 'Test')

for (i in 1:11) {
  test_sign$Variables[i] <- colnames(mod_df)[i+1]
  
  test_sign$Shapiro_P[i] <- shapiro.test(mod_df[,i+1])$p
  test_sign$Levene_P[i] <- levene_test(data = mod_df, formula = mod_df[,i+1] ~ as.factor(PCoA_Clusters))$p
  
  test_sign$Test[i] <- ifelse(test_sign$Shapiro_P[i] > 0.05 & test_sign$Levene_P[i] > 0.05, 'Smthg', 'KW')
} 

KW_res <- as.data.frame(matrix(nrow = 11, ncol = 3)) # Easier to do non parametric tests for all of them
colnames(KW_res) <- c('Variables','Pvalue','EffectSize')

for (i in 1:6) {
  KW_res[i,1] <- colnames(mod_df)[i+1]
  KW_res[i,2] <- kruskal_test(data = mod_df, formula = mod_df[,i+1] ~ PCoA_Clusters)$p
  KW_res[i,3] <- kruskal_effsize(data = mod_df, formula = mod_df[,i+1] ~ PCoA_Clusters)$magnitude
}

pcoapattern <- KW_res %>% filter(Pvalue < 0.05) # All modules except Module 1

mod_to_test <- as.data.frame(mod_df[,c('Mod7','Mod8','Mod9','Mod10','Mod15','PCoA_Clusters')])

pwc_7 <- dunn_test(data = mod_to_test, formula = Mod7 ~ PCoA_Clusters, p.adjust.method = "holm") 
pwc_7 <- pwc_7 %>% add_xy_position(x = "PCoA_Clusters")
pwc_8 <- dunn_test(data = mod_to_test, formula = Mod8 ~ PCoA_Clusters, p.adjust.method = "holm") 
pwc_8 <- pwc_8 %>% add_xy_position(x = "PCoA_Clusters")
pwc_9 <- dunn_test(data = mod_to_test, formula = Mod9 ~ PCoA_Clusters, p.adjust.method = "holm") 
pwc_9 <- pwc_9 %>% add_xy_position(x = "PCoA_Clusters")
pwc_10 <- dunn_test(data = mod_to_test, formula = Mod10 ~ PCoA_Clusters, p.adjust.method = "holm") 
pwc_10 <- pwc_10 %>% add_xy_position(x = "PCoA_Clusters")
pwc_15 <- dunn_test(data = mod_to_test, formula = Mod15 ~ PCoA_Clusters, p.adjust.method = "holm") 
pwc_15 <- pwc_15 %>% add_xy_position(x = "PCoA_Clusters")

bp_Mod_7 <- ggboxplot(mod_df, x = "PCoA_Clusters", y = 'Mod7') +
  stat_pvalue_manual(pwc_7, hide.ns = TRUE)
bp_Mod_8 <- ggboxplot(mod_df, x = "PCoA_Clusters", y = 'Mod8') +
  stat_pvalue_manual(pwc_8, hide.ns = TRUE)
bp_Mod_9 <- ggboxplot(mod_df, x = "PCoA_Clusters", y = 'Mod9') +
  stat_pvalue_manual(pwc_9, hide.ns = TRUE)
bp_Mod_10 <- ggboxplot(mod_df, x = "PCoA_Clusters", y = 'Mod10') +
  stat_pvalue_manual(pwc_10, hide.ns = TRUE)
bp_Mod_15 <- ggboxplot(mod_df, x = "PCoA_Clusters", y = 'Mod15') +
  stat_pvalue_manual(pwc_15, hide.ns = TRUE)

plot_grid(bp_Mod_7, bp_Mod_8, bp_Mod_9, bp_Mod_10, bp_Mod_15,
          nrow = 2, ncol = 3) # just to see the significances

remove(bp_Mod_10, bp_Mod_15, bp_Mod_7, bp_Mod_8, bp_Mod_9,
       pwc_10, pwc_15, pwc_7, pwc_8, pwc_9, 
       wilcoxon, KW_res, test, test_sign)
