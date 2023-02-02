# ANALYSES ON METABARCODING DATA - PCOA + CLUSTERING

remove(list = ls())
set.seed(02121999)

##### LIBRAIRIES
library(tidyverse)
library(vegan)
library(lubridate)
library(funrar)
library(NbClust)
library(pairwiseAdonis)
library(ggrepel)
library(cowplot)
library(dendextend)
library(ggdendro)
library(ggVennDiagram)

##### COLOUR VECTOR
colour <- c('#CDAD00','#458B00','#21918c') # Clusters

################################################################################ DATA IMPORTATION

ASV <- read.csv2("ASV_0-05.csv", header = TRUE, sep = ",")
metadata <- read.csv2("metadata_sample.csv", header = TRUE, sep = ",")

################################################################################ PCOA

## Transpose the dataframe in a matrix
abd_t <- as.data.frame(t(ASV))
colnames(abd_t) <- abd_t[1,]
abd_t <- abd_t[-c(1,dim(abd_t)[1]-1,dim(abd_t)[1]),]
abd_t <- as.data.frame(lapply(abd_t, function(x) if(is.character(x)) as.numeric(x) else x))

## Prepare the dataset for the PCoA
rownames(abd_t) <- colnames(ASV)[2:28]
abd_t_r <- make_relative(as.matrix(abd_t)) # Compute relative abundances
abd_h <- decostand(abd_t_r, method = "hellinger") # Hellinger transformation (square-root)
abd_rh <- as.data.frame(abd_h) # Transform into a dataframe

## Calculate ordination
abd_BC <- vegdist(abd_rh, method = 'bray') # Calculating Bray-Curtis distance matrix
PCoA_BC <- cmdscale(abd_BC, eig = TRUE, k = 7) # Ask for 7 dimension

## Extract the eigenvalues 
eigen <- as.data.frame(PCoA_BC$eig)
eigen$relative <- PCoA_BC$eig/sum(abs(PCoA_BC$eig))
eigen <- round(eigen, digits = 4)
eigen %>% head() # Axis 1 = 22.89% and Axis 2 = 14.82% of the variance explained
sum(eigen$relative[1:7]) # 72.59% of the variance explained in the first 7 axes

################################################################################ CLUSTERING ON PCOA COORDINATES

## How many clusters are recommended ?
NbClust(as.matrix(PCoA_BC$points), diss = NULL, distance = 'euclidean',
        min.nc = 2, max.nc = 10, 
        method = "ward.D2", index ="all") # 3 clusters recommended by 7 indexes
par(mfrow = c(1,1))

## Clustering
clust_pcoa <- hclust(vegdist(as.data.frame(PCoA_BC$points), method = 'euclidean'), method = 'ward.D2')
plot(clust_pcoa) # Indeed, 3 main branches

## Stock the coordinates in a dataset
spg <- as.data.frame(PCoA_BC$points)
spg <- rownames_to_column(spg, var = "Sampling")

## Add the clusters values
spg$AsvClusters <- dendextend::cutree(clust_pcoa, k = 3)
spg$AsvClusters <- factor(spg$AsvClusters, levels = c('3','2','1'), labels = c('1','2','3'))

################################################################################ SIGNIFICANCE OF THE CLUSTERING PATTERNS

## Cluster significance - PermANOVA
adonis(vegdist(as.data.frame(PCoA_BC$points), method = 'euclidean') ~ spg$AsvClusters,
       permutations = how(nperm = 1000)) # Significant p-value = 0.001

## Cluster significance - PairWise tests
pairwise.adonis(vegdist(as.data.frame(PCoA_BC$points), method = 'euclidean'),
                spg$AsvClusters, perm = 1000)  # All significantly different from the others, p-value = 0.003

################################################################################ PLOTS

## Add Time variable from the sampling names
spg$Time <- substring(spg$Sampling, first = 5, last = 5)

## Add all the metadata from the 'cruise' column
spg$Cruise <- substring(spg$Sampling, first = 2, last = 4)
metadata$Cruise <- as.character(metadata$Cruise)
spg <- left_join(spg, metadata[,c('Cruise','Year','Month')], by = 'Cruise')
spg$Label <- paste(month(spg$Month, label = TRUE, locale = 'English'), spg$Year, spg$Time, sep = "-")
spg <- distinct(spg)

## Plot - Scatterplot
PCoA12 <- ggplot(spg) + 
  geom_point(aes(x = V1, y = V2, shape = Time, colour = factor(AsvClusters)), size = 5) + 
  geom_text_repel(aes(x = V1, y = V2, label = substr(Label, start = 1, stop = 8),
                      colour = factor(AsvClusters)), size = 8, show.legend = FALSE) +
  # scale_shape_manual(values = c(15, 16), breaks = c('D', 'N'), labels = c('Day', 'Night')) +
  scale_colour_manual(values = colour) +
  xlab(paste("PCoA 1 (", 100*round(eigen$relative[1], 3),"%)", sep ="")) +
  ylab(paste("PCoA 2 (", 100*round(eigen$relative[2], 3),"%)", sep ="")) +
  labs(colour = 'Clusters') +
  theme_classic(base_size = 25) +
  theme(legend.position = c(0.18, 0.05),
        legend.direction = 'horizontal',
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 25),
        legend.spacing.x = unit(0, 'cm'),
        legend.spacing.y = unit(0, 'cm'),
        legend.margin = margin(t = 0, unit = 'cm'))

## Plot - Dendrogram
plot(clust_pcoa, labels = spg$Label)
rect.hclust(tree = clust_pcoa, k = 3, which = 1:3, border = c('#458B00','#CDAD00','#21918c'), cluster = spg$ClustersPCoA)

dend_pcoa <- as.dendrogram(clust_pcoa)
dend_data <- dendro_data(dend_pcoa, type = 'rectangle')
dend_segment <- as.data.frame(dend_data$segments)
dend_label <- as.data.frame(dend_data$labels)
dend_label <- left_join(dend_label, spg[,c('Sampling','Label','AsvClusters')], by = c('label' = 'Sampling'))
dend_label$Varshort <- substr(dend_label$Label, start = 1, stop = 8)
dend_label$Time <- substr(dend_label$label, start = 5, stop = 5)

Dend <- ggplot() + 
  geom_segment(data = dend_segment, aes(y = -x, x = -y, yend = -xend, xend = -yend)) +
  geom_point(data = dend_label, aes(-y, -x, 
                                    shape = Time, colour = factor(AsvClusters)),
             size = 5) + 
  geom_text(data = dend_label, aes(y + 0.05, -x, label = Varshort, colour = factor(AsvClusters)),
            hjust = -0, angle = 0, size = 8) +
  # scale_shape_manual(values = c(15, 16), breaks = c('D', 'N'), labels = c('Day', 'Night')) +
  scale_colour_manual(values = colour) +
  geom_vline(xintercept = -0.55, linetype = 'dashed', size = 1) +
  xlim(-0.8, 0.25) +
  xlab('Height (Euclidean distance)') +
  guides(colour = 'none', shape = 'none') +
  theme_classic(base_size = 25) +
  theme(axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

## Both of them
plot_grid(PCoA12, Dend, nrow = 1, align = 'h') # Figure 3

## Save clusters as csv
write_csv(spg[,-c(2:8,10:13)], "clusters_pcoa.csv", na = "NA", col_names = TRUE) 
