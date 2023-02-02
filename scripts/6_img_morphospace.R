# SCRIPT IMAGING DATA ANALYSES - ENTIRE COMMUNITY

remove(list = ls())
set.seed(02121999)

##### LIBRAIRIES
library(tidyverse)
library(vegan)
library(lubridate)
library(morphr)
library(FactoMineR)
library(ggrepel)
library(RColorBrewer)
library(igraph)
library(Hmisc)
library(cowplot)
library(ggpubr)
library(rstatix)

##### COLOUR VECTORS
descolour <- colorRampPalette(brewer.pal(name = "Paired", n = 8))(4) # Descriptors
descolour[1] <- "#67B3DB"
descolour[3] <- "#eb4471"

################################################################################ IMPORTATION

community <- read.csv2("community_500um.csv", header = TRUE, sep = ",")
metadata <- read.csv2("metadata_sample.csv", header = TRUE, sep = ",")
pcoa_cluster <- read.csv2("clusters_pcoa.csv", header = TRUE, sep = ",")

################################################################################ FORMATTING

community[,c(5:22,33:44)] <- lapply(community[,c(5:22,33:44)], 
                                    function(x) if(is.character(x) | is.integer(x)) as.numeric(x) else x)
community$date <- as.Date(community$date, format = "%Y-%m-%d") # Date settings
community$sampling <- paste(community$cruise, community$sample_time, sep = '-')
community <- select(community, -sample_time)

################################################################################ ADDITION OF VARIABLES

community_meta <- left_join(community, metadata, by = c('sampling' = 'Name')) # Metadata added

images_directory <- "images_cropped/" # Should be replaced by the path to go to the file containing images
community_meta <- community_meta %>% mutate(img_path = str_c(images_directory, object_id,".png")) # Individual images path (! of the .png)

remove(images_directory, community)

################################################################################ PCA ####

## Select morphological descriptors
pca_var <- community_meta[,5:22] %>%
  mutate(across(.fns=mask_extreme, percent=c(0,0.5))) %>%   # Remove the most extreme high values
  mutate(across(.fns=yeo_johnson)) # Normalise using the Yeo-Johnson transformation

## Add supplementary variables
pca_var$Sample_time <- community_meta$Time # Day/Night 
pca_var$Taxonomy <- community_meta$category # Taxonomic classification
pca_var$Cruise <- paste0(month(community_meta$date, lab = TRUE, locale = 'English'),
                         sep = '-',
                         year(community_meta$date)) # Cruises dates
pca_var$Sampling <- paste0(month(community_meta$date, lab = TRUE, locale = 'English'),
                           sep = '-',
                           year(community_meta$date), sep = '-',
                           pca_var$Sample_time) # Cruises dates with Day-Night distinction

## Create a PCoA Cluster variable
pca_var$Clusters <- ifelse(pca_var$Cruise %in% c('May-2016','Apr-2017'), '2',
                           ifelse(pca_var$Cruise %in% c('Jun-2016','Jul-2016','Aug-2016','Sep-2016','Oct-2016'), '1', '3'))

## Variables linked with the ecosystem functioning
pca_var$Export <- as.numeric(community_meta$C200) # Carbon export at 200-m deep

## Individual weight
vec_weighted <- as.numeric(community_meta$obj_conc) # Weight vector = individual density factor  

## Do the PCA
morpho_pca <- FactoMineR::PCA(pca_var, quali.sup = 19:23, quanti.sup = 24,
                              row.w = vec_weighted, graph = FALSE)

morpho_pca$eig
# 4 axes significant, 87.68% explained (42.868 - 23.89 - 14.53 - 6.59)
remove(vec_weighted)

################################################################################ STOCK THE PCA RESULTS ####

## Morphological descriptors
dscptrs_morpho <- as.data.frame(morpho_pca$var$coord[,1:4])
dscptrs_morpho <- rownames_to_column(dscptrs_morpho, 'descriptors')
rownames(dscptrs_morpho) <- NULL

dscptrs_morpho$colour <- factor(c('Size','Transparency','Transparency','Size',
                                  'Size','Size','Shape','Size','Transparency',
                                  'Transparency','Complexity','Transparency',
                                  'Shape','Shape','Shape','Shape','Complexity',
                                  'Complexity'), 
                                levels = c('Size','Transparency','Shape','Complexity'))

dscptrs_morpho[,1] <- c('Area','MeanGrey','StdDevGrey','Perimeter','MajorAxis',
                        'MinorAxis','Circularity','FeretDiameter','MedianGrey',
                        'SkewnessGrey','Fractal','Grey75%','Symmetry','Symmetry75%',
                        'ThicknessRatio','Elongation','Perim/Feret','Perim/Major') # Give more explicit names

## Supplementary variables
quali_morpho <- as.data.frame(morpho_pca$quali.sup$coord[,1:4])
quali_morpho <- rownames_to_column(quali_morpho, 'Variables')
rowC <- c('Carbon',
          unlist(morpho_pca$quanti.sup$coord[,1]),
          unlist(morpho_pca$quanti.sup$coord[,2]),
          unlist(morpho_pca$quanti.sup$coord[,3]),
          unlist(morpho_pca$quanti.sup$coord[,4]))

supp_morpho <- rbind(quali_morpho, rowC)

remove(quali_morpho, rowC)

## Extract the samplings coordinates for modules study
morphcoord <- supp_morpho[c(46:73),]
metadata$Sampling <- paste0(month(metadata$Month, lab = TRUE, locale = 'English'),
                            sep = '-',
                            metadata$Year, sep = '-',
                            metadata$Time) # Cruises dates with Day-Night distinction
metadata$cruisename <- paste0('B', sep = '',
                              metadata$Cruise, sep = '',
                              substr(metadata$Time, start = 1, stop = 1))
morphcoord_c <- left_join(morphcoord, metadata[,c(17:18)], by = c('Variables' = 'Sampling'))

write.table(morphcoord_c, file = "img_samplingcoord.csv", sep = ",", quote = FALSE, row.names = F) # Save dataframe
remove(morphcoord, morphcoord_c)

################################################################################ CORRELATION AXES/VAR ####

# var_sign <- dimdesc(morpho_pca, axes = 2, proba = 0.01)
# 
# signdf <- as.data.frame(var_sign$Dim.1$quanti)
# signdf <- rownames_to_column(signdf, var = 'Variables')
# colnames(signdf)[2:3] <- c('Cor_Dim1', 'Pval_Dim1')
# 
# signdf2 <- as.data.frame(var_sign$Dim.2$quanti)
# signdf2 <- rownames_to_column(signdf2, var = 'Variables')
# signdf <- left_join(signdf, signdf2)
# colnames(signdf)[4:5] <- c('Cor_Dim2', 'Pval_Dim2')
# 
# signdf2 <- as.data.frame(var_sign$Dim.3$quanti)
# signdf2 <- rownames_to_column(signdf2, var = 'Variables')
# signdf <- left_join(signdf, signdf2)
# colnames(signdf)[6:7] <- c('Cor_Dim3', 'Pval_Dim3')
# 
# signdf2 <- as.data.frame(var_sign$Dim.4$quanti)
# signdf2 <- rownames_to_column(signdf2, var = 'Variables')
# signdf <- left_join(signdf, signdf2)
# colnames(signdf)[8:9] <- c('Cor_Dim4', 'Pval_Dim4')
# 
# write.table(signdf, file = "sign_morphospace.csv", sep = ",", quote = FALSE, row.names = F)
# remove(signdf2,signdf)

################################################################################ IMAGES PLOTS - Descriptors ####

supp2 <- supp_morpho[c(1:2,74:77),c(1:5)]
supp2[,c(2:5)] <- lapply(supp2[,c(2:5)], function(x) if(is.character(x) | is.integer(x)) as.numeric(x) else x)
supp2$Legend <- c('Day','Night',
                  'Summer community','Bloom community','Winter mixing community',
                  'Carbon export')
colsupp2 <- c('grey75','black',
              '#CDAD00','#458B00','#21918c',
              'black') # Colour vector for Day, Night, 3 clusters, C export

## PC 12
imagesa <- ggmorph_tile(morpho_pca, community_meta$img_path, dimensions = c(1,2), steps = 16, n_imgs = 5, scale = 0.005) # Images background PC12
imagesb <- ggmorph_tile(morpho_pca, community_meta$img_path, dimensions = c(3,4), steps = 16, n_imgs = 5, scale = 0.005) # Images background PC34

pc12d <- imagesa +
  geom_segment(data = dscptrs_morpho, 
               aes(x = 0, y = 0, xend = Dim.1*5, yend = Dim.2*5, color = colour), 
               size = 1.5, alpha = 0.9,
               arrow = arrow(length = unit(0.2,"cm"))) + # Descriptors arrows
  geom_text_repel(data = dscptrs_morpho,
                  aes(x = Dim.1*5, y = Dim.2*5, label = descriptors, color = colour),
                  fontface = 1, alpha = 1, segment.alpha = 0.2, size = 10, show.legend = F) + # Descriptors labels
  geom_segment(data = supp2[6,], 
               aes(x = 0, y = 0, xend = Dim.1*5, yend = Dim.2*5),
               colour = colsupp2[6],
               size = 1.5, alpha = 0.9,
               arrow = arrow(length = unit(0.2,"cm"))) + # Export arrow
  geom_text_repel(data = supp2[6,], 
                  aes(x = Dim.1*5, y = Dim.2*5 + 0.1, label = Legend),
                  colour = colsupp2[6], 
                  size = 12, show.legend = FALSE, 
                  fontface = 1, alpha = 1, segment.alpha = 0.2) + # Label
  theme_classic(base_size = 30) +
  scale_colour_manual(values = descolour) +
  xlab(paste("PC1 (", round(morpho_pca$eig[1,2], 3),"%)", sep ="")) +
  ylab(paste("PC2 (", round(morpho_pca$eig[2,2], 3),"%)", sep ="")) +
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 25),
        legend.direction = 'vertical',
        legend.position = c(0.15,0.10),
        legend.text = element_text(size = 25),
        legend.title = element_blank()) +
  scale_size(range = c(1, 10))

## PC 34
pc34d <- imagesb +
  geom_segment(data = dscptrs_morpho, 
               aes(x = 0, y = 0, xend = Dim.3*5, yend = Dim.4*5, color = colour), 
               size = 1.5, alpha = 0.9,
               arrow = arrow(length = unit(0.2,"cm"))) + # Descriptors arrows
  geom_text_repel(data = dscptrs_morpho,
                  aes(x = Dim.3*5, y = Dim.4*5, label = descriptors, color = colour),
                  fontface = 1, alpha = 1, segment.alpha = 0.2, size = 10) + # Descriptors labels
  geom_segment(data = supp2[6,], 
               aes(x = 0, y = 0, xend = Dim.3*5, yend = Dim.4*5),
               colour = colsupp2[6],
               size = 1.5, alpha = 0.9, show.legend = F,
               arrow = arrow(length = unit(0.2,"cm"))) + # Export arrow
  geom_text_repel(data = supp2[6,], 
                  aes(x = Dim.3*5 + 2, y = Dim.4*5, label = Legend),
                  colour = colsupp2[6], 
                  size = 12, show.legend = FALSE, 
                  fontface = 1, alpha = 1, segment.alpha = 0.2) + # Label
  theme_classic(base_size = 30) +
  scale_colour_manual(values = descolour) +
  xlab(paste("PC3 (", round(morpho_pca$eig[3,2], 2),"%)", sep ="")) +
  ylab(paste("PC4 (", round(morpho_pca$eig[4,2], 2),"%)", sep ="")) +
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 25)) +
  scale_size(range = c(1, 10)) +
  guides(size = "none", color = "none", shape = "none")

plot_grid(pc12d, pc34d, align = 'hv', nrow = 1) # Figure 1
remove(pc12d, pc34d)

################################################################################ IMAGES PLOTS - Taxonomy ####
nimg <- community_meta %>% 
  group_by(category) %>%
  summarise(n = n())
nimg$logn <- log10(nimg$n)
nimg$Type <- factor(c("Other","Crustacea","Other","Mollusca","Other","Copepoda","Gastropoda",
               "Other","Mollusca","Gelatinous","Crustacea",
               "Crustacea","Gelatinous","Gastropoda","Copepoda","Crustacea","Gelatinous",
               "Other","Other","Crustacea","Gastropoda","Gastropoda","Copepoda",
               "Other","Gelatinous","Gastropoda","Crustacea","Crustacea","Gelatinous"))
supp_morpho <- left_join(supp_morpho, nimg, by = c('Variables' = 'category'))
colnames(supp_morpho)[6] <- 'Number of images'
colnames(supp_morpho)[7] <- 'Number of images, log10'
supp_morpho[,c(2:5)] <- lapply(supp_morpho[,c(2:5)], function(x) if(is.character(x) | is.integer(x)) as.numeric(x) else x)

tax12 <- imagesa +
  geom_point(data = supp_morpho[3:31,],
             aes(x = Dim.1, y = Dim.2, size = `Number of images, log10`, colour = Type), 
             shape = 18, alpha = 0.75) +
  geom_text_repel(data = supp_morpho[3:31,],
                  aes(x = Dim.1, y = Dim.2, label = Variables, colour = Type),
                  max.overlaps = Inf, size = 8, show.legend = FALSE) +
  scale_colour_brewer(palette = "Dark2") +
  scale_radius(limits = c(0, max(supp_morpho[3:31,7])), range = c(3, 8)) +
  theme_classic(base_size = 30) +
  theme(legend.position = c(0.45,0.08),
        legend.direction = 'horizontal',
        legend.spacing.x = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.5, 'cm'),
        legend.margin = margin(t = 0, unit = 'cm'),
        axis.title = element_text(size = 25),
        axis.text = element_text(size = 25)) +
  guides(colour = guide_legend(override.aes = list(size = (6))))

tax34 <- imagesb +
  geom_point(data = supp_morpho[3:31,],
             aes(x = Dim.3, y = Dim.4, size = `Number of images, log10`, colour = Type), 
             shape = 18, alpha = 0.75) +
  geom_text_repel(data = supp_morpho[3:31,],
                  aes(x = Dim.3, y = Dim.4, label = Variables, colour = Type), 
                  max.overlaps = Inf, size = 8) +
  scale_colour_brewer(palette = "Dark2") +
  scale_radius(limits = c(0, NA), range = c(3, 8)) +
  theme_classic(base_size = 30) +
  theme(legend.position = 'none',
        axis.title = element_text(size = 25),
        axis.text = element_text(size = 25))

plot_grid(tax12, tax34, align = 'hv', nrow = 1) # Supp. Figure III

remove(tax12, tax34)

################################################################################ Individual properties ####

ind_morpho <- as.data.frame(morpho_pca$ind$coord[,1:4])

## Create a Cruise variable
ind_morpho$Cruises <- paste0(month(as.Date(community_meta$date, format = "%Y-%m-%d"), lab = TRUE, locale = 'English'),
                             sep = '-',
                             year(as.Date(community_meta$date, format = "%Y-%m-%d")))
ind_morpho$Cruises <- factor(ind_morpho$Cruises, levels = c('Mar-2016','Apr-2016',"May-2016",
                                                            "Jun-2016","Jul-2016","Aug-2016",
                                                            "Sep-2016","Oct-2016","Nov-2016",
                                                            "Dec-2016","Jan-2017","Feb-2017",
                                                            "Apr-2017","May-2017"))
ind_morpho$obj_conc <- community_meta$obj_conc

## Create a PCoA Cluster variable
ind_morpho$Clusters <- ifelse(ind_morpho$Cruises %in% c('May-2016','Apr-2017'), '2',
                              ifelse(ind_morpho$Cruises %in% c("Jun-2016","Jul-2016","Aug-2016","Sep-2016","Oct-2016"), '1', '3'))

supp2$factor <- c(1,1,1,1,1,5)

################################################################################ Clusters distribution along PC1 and PC2 ####

# PC12 morphospace
mainpc12 <- imagesa + 
  geom_point(data = supp2[3:5,], 
             aes(x = Dim.1, y = Dim.2), 
             shape = 15, colour = colsupp2[3:5],
             size = 5, show.legend = FALSE) + # Clusters points
  geom_segment(data = supp2[6,], 
               aes(x = 0, y = 0, xend = Dim.1*5, yend = Dim.2*5),
               colour = colsupp2[6],
               size = 1.5, alpha = 0.9, show.legend = F,
               arrow = arrow(length = unit(0.2,"cm"))) + # Export arrows    
  geom_text_repel(data = supp2[3:5,], 
                  aes(x = Dim.1, y = Dim.2, label = Legend),
                  colour = colsupp2[3:5], 
                  size = 12, show.legend = FALSE, max.overlaps = Inf,
                  fontface = 1, alpha = 1, segment.alpha = 0.2) + # Labels
  geom_text(data = supp2[6,], 
            aes(x = Dim.1*5, y = Dim.2*5 + 0.3, label = Legend),
            colour = colsupp2[2], 
            size = 12, show.legend = FALSE,
            fontface = 1, alpha = 1) + # Labels
  theme_classic(base_size = 30) +
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 25),
        plot.margin = unit(c(0,0,0,0), "pt")) +
  guides(color = 'none', shape = 'none')

# Clusters along PC1
pc1_c <- ggplot(ind_morpho) +
  geom_freqpoly(aes(x = Dim.1,
                    weight = obj_conc,
                    colour = Clusters),
                binwidth = 0.2, position = 'identity',
                alpha = 0.8, size = 1.5) +
  scale_colour_manual(values = colsupp2[3:5]) +
  theme_classic(base_size = 30) +
  ylab(bquote(Density (ind.m^-3))) +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title = element_text(size = 25),
        axis.text = element_text(size = 22),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))

# Clusters along PC2
pc2_c <- ggplot(ind_morpho) + 
  geom_freqpoly(aes(x = Dim.2,
                    weight = obj_conc,
                    colour = Clusters), 
                binwidth = 0.2, position = 'identity',
                alpha = 0.8, size = 1.5) +
  scale_colour_manual(values = colsupp2[3:5]) +
  theme_classic(base_size = 30) +
  xlab('PC2 coordinates') +
  ylab(bquote(Density (ind.m^-3))) +
  theme(legend.position = 'none',
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        axis.title = element_text(size = 25),
        axis.text = element_text(size = 22),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) + 
  rotate()

# All together
plot_grid(pc1_c, NULL, NULL,
          NULL, NULL, NULL,
          mainpc12, NULL, pc2_c,
          ncol = 3, nrow = 3, align = 'hv',
          rel_widths = c(3,-0.2,1), rel_heights = c(1,-0.15,3)) # Figure 4

remove(mainpc12,pc1_c,pc2_c)

################################################################################ Clusters distribution along PC3 and PC4 ####

# PC34 morphospace
mainpc34 <- imagesb + 
  geom_point(data = supp2[3:5,], 
             aes(x = Dim.3*factor, y = Dim.4*factor),
             shape = 15, colour = colsupp2[3:5],
             size = 5, show.legend = FALSE) + # Clusters points
  geom_segment(data = supp2[6,], 
               aes(x = 0, y = 0, xend = Dim.3*factor, yend = Dim.4*factor),
               colour = colsupp2[2],
               size = 1.5, alpha = 0.9, show.legend = F,
               arrow = arrow(length = unit(0.2,"cm"))) + # DW/Export arrows    
  geom_text_repel(data = supp2[3:6,], 
                  aes(x = Dim.3*factor, y = Dim.4*factor, label = Legend),
                  colour = colsupp2[3:6], 
                  size = 12, show.legend = FALSE,
                  fontface = 1, alpha = 1, segment.alpha = 0.2) + # Labels
  theme_classic(base_size = 30) +
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 25),
        plot.margin = unit(c(0,0,0,0), "pt")) +
  guides(color = 'none', shape = 'none')

# Clusters along PC3
pc3_c <- ggplot(ind_morpho) +
  geom_freqpoly(aes(x = Dim.3,
                    weight = obj_conc,
                    colour = Clusters),
                binwidth = 0.2, position = 'identity',
                alpha = 0.8, size = 1.5) +
  scale_colour_manual(values = colsupp2[3:5]) +
  theme_classic(base_size = 30) +
  ylab(bquote(Density (ind.m^-3))) +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title = element_text(size = 25),
        axis.text = element_text(size = 22),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))

# Clusters along PC4
pc4_c <- ggplot(ind_morpho) + 
  geom_freqpoly(aes(x = Dim.4,
                    weight = obj_conc,
                    colour = Clusters), 
                binwidth = 0.2, position = 'identity',
                alpha = 0.8, size = 1.5) +
  scale_colour_manual(values = colsupp2[3:5]) +
  theme_classic(base_size = 30) +
  scale_x_continuous(name = 'PC4 coordiantes',
                     limits = c(-8.5,5)) +
  ylab(bquote(Density (ind.m^-3))) +
  theme(legend.position = 'none',
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        axis.title = element_text(size = 25),
        axis.text = element_text(size = 22),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) + 
  rotate()

# All together
plot_grid(pc3_c, NULL, NULL,
          NULL, NULL, NULL,
          mainpc34, NULL, pc4_c,
          ncol = 3, nrow = 3, align = 'hv',
          rel_widths = c(3,-0.2,1), rel_heights = c(1,-0.15,3)) # Supp. Figure II

remove(mainpc34,pc3_c,pc4_c)
