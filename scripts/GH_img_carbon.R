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

##### ENTIRE COMMUNITY - DATASETS IMPORTATION AND CLEANING

################################################################################ IMPORTATION

community <- read.csv2("community_500um.csv", header = TRUE, sep = ",")
metadata <- read.csv2("metadata_sample.csv", header = TRUE, sep = ",")

################################################################################ FORMATTING

community[,c(5:22,33:44)] <- lapply(community[,c(5:22,33:44)], 
                                    function(x) if(is.character(x) | is.integer(x)) as.numeric(x) else x)
community$date <- as.Date(community$date, format = "%Y-%m-%d")
community$sampling <- paste(community$cruise, community$sample_time, sep = '-')
community <- select(community, -sample_time)

################################################################################ ADDITION OF METADATA

community_meta <- left_join(community, metadata, by = 'sampling')

################################################################################ TRAITS

traitc <- community_meta[,c('stddev','skew','mean','median','hist75',
                            'minor','major','perim','feret','area',
                            'C200','cruise','obj_conc')] # Keep only opacity- and size-related traits

traitr <- traitc[,1:10] * traitc$obj_conc # Multiply each individual value by its density factor d
traitr$obj_conc <- traitc$obj_conc
traitr$C200 <- as.numeric(traitc$C200)
traitr$cruise <- traitc$cruise

traitr_m <- traitr %>%
  group_by(cruise) %>%
  summarise_at(vars(stddev:obj_conc), sum) # Sum all the individual values to get a total per cruise (since we have one export value per cruise)
traitr_m[2:11] <- traitr_m[2:11] / traitr_m$obj_conc # Compute the average value per cruise
traitr_m$C200 <- unique(traitr$C200)

trait <- traitr_m[,-1] %>% pivot_longer(cols = -'C200', names_to = 'traits', values_to = 'value')
trait$type <- ifelse(trait$traits %in% c('stddev','skew','mean','median','hist75'), 'opacity', 'size') # Create a 'type' of traits variable

## Opacity
opacity <- trait %>% filter(type == 'opacity')
opacity$C200log <- log(opacity$C200)

## Size
size <- trait %>% filter(type == 'size') %>% filter(traits != 'obj_conc')
size$C200log <- log(size$C200)

################################################################################ LINEAR RELATIONSHIPS

#### Opacity
lin_opac <- as.data.frame(matrix(ncol = 5, nrow = length(unique(opacity$traits))))
colnames(lin_opac) <- c('trait','pearson','pvaluep','shapiroCE','shapiroMT')

for (i in 1:length(unique(opacity$traits))) { # carbon export transformed values
  j <- unique(opacity$traits)[i]
  lin_opac[i,1] <- j
  fake <- opacity %>% filter(traits == j)
  lin_opac[i,2] <- cor.test(as.numeric(unlist(fake[,3])), as.numeric(unlist(fake[,5])),  method = "pearson")$estimate
  lin_opac[i,3] <- cor.test(as.numeric(unlist(fake[,3])), as.numeric(unlist(fake[,5])),  method = "pearson")$p.value
  lin_opac[i,6] <- shapiro.test(as.numeric(unlist(fake[,5])))$p.value
  lin_opac[i,7] <- shapiro.test(as.numeric(unlist(fake[,3])))$p.value
}
remove(i, j, fake)

lin_opac_s <- lin_opac %>% filter(pvaluep < 0.05) # Linear relationship ok because normality ok and pearson coefficient significant except for skewness

#### Size
lin_size <- as.data.frame(matrix(ncol = 5, nrow = length(unique(size$traits))))
colnames(lin_size) <- c('trait','pearson','pvaluep','shapiroCE','shapiroMT')

for (i in 1:length(unique(size$traits))) { # carbon export transformed values
  j <- unique(size$traits)[i]
  lin_size[i,1] <- j
  fake <- size %>% filter(traits == j)
  lin_size[i,2] <- cor.test(as.numeric(unlist(fake[,5])), as.numeric(unlist(fake[,3])),  method = "pearson")$estimate
  lin_size[i,3] <- cor.test(as.numeric(unlist(fake[,5])), as.numeric(unlist(fake[,3])),  method = "pearson")$p.value
  lin_size[i,6] <- shapiro.test(as.numeric(unlist(fake[,5])))$p.value
  lin_size[i,7] <- shapiro.test(as.numeric(unlist(fake[,3])))$p.value
}
remove(i, j, fake)

lin_size_s <- lin_size %>% filter(pvaluep < 0.05) # No linear relationship 

################################################################################ PLOT

# Traits' names
opacity$name <- factor(opacity$traits, 
                       levels = c('mean','median','stddev','hist75','skew'),
                       labels = c("MeanGrey","MedianGrey","StdDevGrey","Grey75%","SkewnessGrey")) # Make more strait forward names for the plot
size$name <- factor(size$traits, 
                    levels = c('minor','major','perim','feret','area'),
                    labels = c("MinorAxis","MajorAxis","Perimeter","FeretDiameter","Area")) # Make more strait forward names for the plot

# Merge to have only one plot to do
alltraits <- rbind(size, opacity)

# Plot
ggplot(alltraits, aes(x = value, y = C200log)) +
  geom_smooth(data = alltraits %>% filter(type == 'opacity'),
              aes(x = value, y = C200log), 
              method = 'lm', formula = y~x, colour = 'black') +
  geom_point(data = alltraits,
             aes(x = value, y = C200log),
             size = 4, shape = 19) +
  stat_cor(label.x.npc = 0.10,
           label.y.npc = 0.11, 
           p.accuracy = 0.001, r.accuracy = 0.01, size = 8) +
  facet_wrap(~name, scales = 'free_x', ncol = 5, ) +
  ylab('Carbon export at 200 m (ln)') +
  xlab('Mean trait value per cruise') +
  theme_classic(base_size = 25) +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.title = element_blank(),
        legend.text = element_text(size = 25)) # Figure 2
