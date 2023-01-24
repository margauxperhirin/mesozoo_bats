# SCRIPT IMAGING DATA ANALYSES - ENTIRE COMMUNITY

remove(list = ls())
set.seed(02121999)

##### LIBRAIRIES
library(tidyverse)
library(lubridate)
library(corrplot)
library(Hmisc)
library(ggpubr)

##### IMPORTATION

IMG <- read.csv2("community_500um.csv", header = TRUE, sep = ",")
ASVtot <- read.csv2("ASV_rarefy.csv", header = TRUE, sep = ",")

################################################################################ IMAGES - Biomass per sampling and taxonomic category
IMG[,c(5:22,33:44)] <- lapply(IMG[,c(5:22,33:44)], 
                              function(x) if(is.character(x) | is.integer(x)) as.numeric(x) else x)

IMG$date <- as.Date(IMG$date, format = "%Y-%m-%d")
IMG$sampling <- paste0(month(IMG$date, lab = TRUE, locale = 'English'),
                       sep = '-',
                       year(IMG$date), sep = '-',
                       IMG$sample_time) # Cruises dates with Day-Night distinction

IMG$merge <- paste0('B', sep = '',
                    IMG$cruise, sep = '',
                    substr(IMG$sample_time, start = 1, stop = 1)) # Variable to merge both datasets

IMG$category <- ifelse(IMG$category == 'Larvacea', 'Appendicularia', IMG$category) # Rename Appendicularian in Larvacea

IMG$bmESDnorm <- IMG$dwESD*IMG$obj_conc

r_IMG <- IMG %>%
  group_by(sampling) %>%
  summarise(br = sum(bmESDnorm),
            ar = sum(n()*obj_conc))

b_IMG <- IMG %>% 
  group_by(category, sampling) %>%
  summarise(bmESD = sum(bmESDnorm),
            abd = sum(n()*obj_conc))

b_IMG$date <- substr(b_IMG$sampling, start = 1, stop = 8)
b_IMG$date <- factor(b_IMG$date, 
                     levels = c('Mar-2016','Apr-2016','May-2016','Jun-2016','Jul-2016',
                                'Aug-2016','Sep-2016','Oct-2016','Nov-2016','Dec-2016',
                                'Jan-2017','Feb-2017','Apr-2017','May-2017'))
b_IMG$time <- substr(b_IMG$sampling, start = 10, stop = 15)
b_IMG <- b_IMG %>% filter(category != 'Egg') # Remove Egg cateogry
b_IMG <- b_IMG %>% filter(sampling != 'Oct-2016-Day') # remove Oct-Day because no metaB data for it
b_IMG <- unique(b_IMG)

br_IMG <- left_join(b_IMG, r_IMG, by = 'sampling')
br_IMG$abd_r <- br_IMG$abd / br_IMG$ar # Make relative abundances
br_IMG$bmESD_r <- br_IMG$bmESD / br_IMG$br # Make relative biomasses

################################################################################ METAB - Number of reads per sampling and taxonomic category

ASVr <- ASVtot
for (i in 2:28) {
  sum <- sum(as.numeric(ASVr[,i]))
  ASVr[,i] <- ASVr[,i] / sum
} # Make relative abundances

# Filter on phylum
ASVp <- ASVr %>%
  filter(str_detect(Phylum, 'Annelida|Bryozoa|Chaetognatha|Cnidaria|Echinodermata')) # 3167 ASVs
# Filter on subphylum
ASVsp <- ASVr %>%
  filter(str_detect(Subphylum, 'Cephalochordata')) # 54 ASV
# Filter on class - Fish removed when Vertebrata removed
ASVc <- ASVr %>%
  filter(str_detect(Class, 'Ostracoda|Actinopterygii|Bivalvia|Cephalopoda|Appendicularia')) # 5782 ASVs
# Filter on infraclass
ASVic <- ASVr %>%
  filter(str_detect(Infraclass, 'Cirripedia')) # 5 ASV
# Filter on super order
ASVso <- ASVr %>%
  filter(str_detect(Superorder, 'Cladocera')) # 319 ASVs
# Filter on order
ASVo <- ASVr %>%
  filter(str_detect(Order, 'Amphipoda|Calanoida|Cyclopoida|Harpacticoida|Decapoda|Euphausiacea|Doliolida|Salpida')) # 10241 ASVs
# Filter on suborder
ASVos <- ASVr %>%
  filter(str_detect(Suborder, 'Gymnosomata')) # 92 ASVs
# Filter on super-family
ASVsf <- ASVr %>%
  filter(str_detect(Superfamily, 'Pterotracheoidea')) # 0 ASV
# Filter on family
ASVf <- ASVr %>%
  filter(str_detect(Family, 'Cavoliniidae|Creseidae|Limacinidae')) # 272 ASVs

# Bind everything
ASVta <- rbind(ASVp,ASVsp,ASVc,ASVic,ASVso,ASVo,ASVos,ASVsf,ASVf) # 19932 ASVs concerned out of 22814
remove(ASVp,ASVsp,ASVc,ASVic,ASVso,ASVo,ASVos,ASVsf,ASVf,ASVtot)

# Clean and prepare the dataset
ASVta <- ASVta[,-c(1,49)]
ASVta$category <- with(ASVta, ifelse(str_detect(Phylum, 'Annelida') == TRUE, 'Annelida',
                                     ifelse(str_detect(Phylum, 'Bryozoa') == TRUE, 'Bryozoa',
                                            ifelse(str_detect(Phylum, 'Chaetognatha') == TRUE, 'Chaetognatha',
                                                   ifelse(str_detect(Phylum, 'Cnidaria') == TRUE, 'Cnidaria',
                                                          ifelse(str_detect(Phylum, 'Echinodermata') == TRUE, 'Echinodermata',
                                                                 ifelse(str_detect(Subphylum, 'Cephalochordata') == TRUE, 'Cephalochordata',
                                                                        ifelse(str_detect(Class, 'Ostracoda') == TRUE, 'Ostracoda',
                                                                               ifelse(str_detect(Class, 'Actinopterygii') == TRUE, 'Actinopterygii',
                                                                                      ifelse(str_detect(Class, 'Bivalvia') == TRUE, 'Bivalvia',
                                                                                             ifelse(str_detect(Class, 'Cephalopoda') == TRUE, 'Cephalopoda',
                                                                                                    ifelse(str_detect(Class, 'Appendicularia') == TRUE, 'Appendicularia',
                                                                                                           ifelse(str_detect(Infraclass, 'Cirripedia') == TRUE, 'Cirripedia',
                                                                                                                  ifelse(str_detect(Superorder, 'Cladocera') == TRUE, 'Cladocera',
                                                                                                                         ifelse(str_detect(Order, 'Amphipoda') == TRUE, 'Amphipoda',
                                                                                                                                ifelse(str_detect(Order, 'Calanoida') == TRUE, 'Calanoida',
                                                                                                                                       ifelse(str_detect(Order, 'Cyclopoida') == TRUE, 'Cyclopoida',
                                                                                                                                              ifelse(str_detect(Order, 'Harpacticoida') == TRUE, 'Harpacticoida',
                                                                                                                                                     ifelse(str_detect(Order, 'Decapoda') == TRUE, 'Decapoda',
                                                                                                                                                            ifelse(str_detect(Order, 'Euphausiacea') == TRUE, 'Euphausiacea',
                                                                                                                                                                   ifelse(str_detect(Order, 'Doliolida') == TRUE, 'Doliolida',
                                                                                                                                                                          ifelse(str_detect(Order, 'Salpida') == TRUE, 'Salpida',
                                                                                                                                                                                 ifelse(str_detect(Suborder, 'Gymnosomata') == TRUE, 'Gymnosomata',
                                                                                                                                                                                        ifelse(str_detect(Superfamily, 'Pterotracheoidea') == TRUE, 'Pterotracheoidea',
                                                                                                                                                                                               ifelse(str_detect(Family, 'Cavoliniidae') == TRUE, 'Cavoliniidae',
                                                                                                                                                                                                      ifelse(str_detect(Family, 'Creseidae') == TRUE, 'Creseidae',
                                                                                                                                                                                                             ifelse(str_detect(Family, 'Limacinidae') == TRUE, 'Limacinidae',NA)))))))))))))))))))))))))))
ASVta <- ASVta[,-c(28:47)]
r_ASV <- ASVta %>% pivot_longer(cols = !category, names_to = 'cruise', values_to = 'nreads')
r_ASV <- r_ASV %>%
  group_by(category, cruise) %>%
  summarise(sreads = sum(nreads)) # Get a total of reads per cruise and category

r_ASV <- left_join(r_ASV, IMG[,c('sample_time','sampling','merge')], by = c('cruise' = 'merge'))
r_ASV <- unique(r_ASV)
r_ASV <- r_ASV[,-2]
r_ASV$date <- substr(r_ASV$sampling, start = 1, stop = 8)
r_ASV$date <- factor(r_ASV$date, levels = c('Mar-2016','Apr-2016','May-2016','Jun-2016','Jul-2016',
                                            'Aug-2016','Sep-2016','Oct-2016','Nov-2016','Dec-2016',
                                            'Jan-2017','Feb-2017','Apr-2017','May-2017'))

## Final dataset with nreads and biomass
rb <- left_join(r_ASV, br_IMG[,c(1:2,9:10)], by = c('category','sampling'))
rb_na <- na.omit(rb) # 621 observations including 118 without nreads or biomass (need to remove them for the comparison)

################################################################################ Reads - Biomass from images

cor_rb <- as.data.frame(matrix(ncol = 7, nrow = length(unique(rb_na$category)))) # Create a dataframe to save the results
colnames(cor_rb) <- c('category','pearson','pvaluep','shapiroIMG','shapiroASV')

for (i in 1:length(unique(rb_na$category))) {
  j <- unique(rb_na$category)[i]
  cor_rb[i,1] <- j
  rb_fake <- rb_na %>% filter(category == j)
  cor_rb[i,2] <- cor.test(as.numeric(unlist(rb_fake[,2])), as.numeric(unlist(rb_fake[,7])),  method = "pearson")$estimate
  cor_rb[i,3] <- cor.test(as.numeric(unlist(rb_fake[,2])), as.numeric(unlist(rb_fake[,7])),  method = "pearson")$p.value
  cor_rb[i,6] <- shapiro.test(as.numeric(unlist(rb_fake[,2])))$p.value
  cor_rb[i,7] <- shapiro.test(as.numeric(unlist(rb_fake[,7])))$p.value
}
remove(i, j, rb_fake)

cors_rbnolog <- cor_rb %>% filter(pvaluep < 0.05)

################################################################################ Reads - Abundance from images

cor_abd <- as.data.frame(matrix(ncol = 7, nrow = length(unique(rb_na$category)))) # Create a dataframe to save the results
colnames(cor_abd) <- c('category','pearson','pvaluep','shapiroIMG','shapiroASV')

for (i in 1:length(unique(rb_na$category))) {
  j <- unique(rb_na$category)[i]
  cor_abd[i,1] <- j
  fake <- rb_na %>% filter(category == j)
  cor_abd[i,2] <- cor.test(as.numeric(unlist(fake[,2])), as.numeric(unlist(fake[,6])),  method = "pearson")$estimate
  cor_abd[i,3] <- cor.test(as.numeric(unlist(fake[,2])), as.numeric(unlist(fake[,6])),  method = "pearson")$p.value
  cor_abd[i,6] <- shapiro.test(as.numeric(unlist(fake[,2])))$p.value
  cor_abd[i,7] <- shapiro.test(as.numeric(unlist(fake[,6])))$p.value
}
remove(i, j, fake)

cors_abdnolog <- cor_abd %>% filter(pvaluep < 0.05)
