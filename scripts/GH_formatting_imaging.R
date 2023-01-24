# SCRIPT IMAGING DATA FORMATING 

remove(list = ls(all.names = TRUE))
set.seed(02121999)

##### LIBRAIRIES
library(tidyverse)
library(lubridate)

################################################################################ IMPORTATION

images <- read.table("ecotaxa_export.tsv", header = TRUE, sep = "\t")

################################################################################ FORMATTING


colnames(images)[1] <- "object_id" ## Rename the first column
images <- images %>% filter(object_date != "2016-03-23") ## Remove cruise 322 (late March-16)
images$acq_sub_part[images$acq_sub_part == 999999] <- 16 # Correct the missing data (checked with Hannah Gossner)

## Select only mesozooplankton categories
throw <- c('fiber<detritus', 'not-living', 'badfocus<artefact', 'artefact', 'detritus',
           'multiple<other', 'other<living', 'multiple<Copepoda', 'Copepoda',
           'Crustacea', 'Harosa', 'bubble', 'Mollusca')

community <- images %>% filter(!object_annotation_category %in% throw) # 72965 images left 

## Select only 'validated' images
community <- community %>% filter(object_annotation_status == "validated") # 72887 images left 

## Select only variables we need 
keep <- c('object_id', 'object_date', 'object_time', 'object_annotation_category', 'object_area', 
          'object_mean', 'object_stddev', 'object_perim.', 'object_major', 'object_minor', 
          'object_circ.', 'object_feret', 'object_median', 'object_skew', 'object_fractal',
          'object_histcum1', 'object_symetriev', 'object_symetrievc', 'object_thickr', 
          'object_elongation', 'object_perimferet', 'object_perimmajor', 'sample_id', 
          'sample_tot_vol', 'acq_sub_part', 'object_depth_max')

pixel_mm <- unique(community$process_particle_pixel_size_mm) # factor for pixel size in mm
dpi_mm <- 0.00002809 # scan settings dpi to mm

community_small <- community[keep]

## Rename the descriptors into easier names
names(community_small)[1:22] <- c("object_id", "date", "time", "category",
                                  "area","mean","stddev","perim","major","minor","circ",
                                  "feret","median","skew","fractal","hist75","symv",
                                  "symvc","thickr","elongation","perimferet","perimmajor")

## Rename categories into easier names
community_small$category[community_small$category == 'Cnidaria<Metazoa'] <- 'Cnidaria'
community_small$category[community_small$category == 'egg<other'] <- 'Egg'
community_small$category[community_small$category == 'nauplii<Crustacea'] <- 'Nauplii'
community_small$category[community_small$category == 'Bivalvia<Mollusca'] <- 'Bivalvia'
community_small$category[community_small$category == 'Oikopleura'] <- 'Larvacea'

unique(community_small$category) # it works

################################################################################ CREATE NEW VARIABLES
## Cruise number
community_small$cruise <- substring(community_small$sample_id, first = 2, last = 4)

## Day/Night
community_small$time <- as.numeric(paste(substring(community_small$time, first = 1, last = 2), 
                                         substring(community_small$time, first = 4, last = 5), sep = ""))
community_small$sample_time <- NA
community_small$sample_time <- ifelse(community_small$time > 800 & community_small$time < 1900, "Day", community_small$sample_time) # Day = between 8AM and 7PM
community_small$sample_time <- ifelse(community_small$time < 800 | community_small$time > 1900, "Night", community_small$sample_time) # Night = between 7PM and 8AM
community_small$sample_time <- factor(community_small$sample_time, levels = c("Day", "Night"))
summary(community_small$sample_time) # it works

## Sampling date
community_small$date <- as.Date(as.character(community_small$date), format = "%Y-%m-%d")
community_small$year <- year(community_small$date)
community_small$month <- month(community_small$date, label = TRUE, locale = "English") # new variable with the sampling month
community_small$month_number <- month(community_small$date, label = FALSE)
community_small$day_of_year <- yday(community_small$date)
community_small$day <- day(community_small$date)

## Biovolume
community_small$area_mm2 <- community_small$area*dpi_mm 
community_small$major_mm <- community_small$major*pixel_mm 
community_small$minor_mm <- community_small$minor*pixel_mm 
community_small$ESD <- 2*sqrt(community_small$area_mm2 / pi)
community_small$volESD <- (4/3)*pi*((community_small$ESD / 2)^3)
community_small$volume <- (4/3)*pi*((community_small$minor_mm/2)^2)*(community_small$major_mm/2) # Ellipsoidal shape for copepods
community_small$log_vol <- as.numeric(log10(community_small$volume)) # Logarithm of the volume

## Dry weight from individual ESD volumes
# Factors from Maas et al., 2022
community_small$dwESD <- with(community_small, ifelse(category == "Calanoida", 0.055*volESD,
                                                      ifelse(category == "Cyclopoida", 0.055*volESD,
                                                             ifelse(category == "Harpacticoida", 0.055*volESD,
                                                                    ifelse(category == "Chaetognatha", 0.013*volESD,
                                                                           ifelse(category == "Ostracoda", 0.052*volESD,
                                                                                  ifelse(category == "Gastropoda",0.1913*volESD,
                                                                                         ifelse(category == "Heteropoda",0.1913*volESD,
                                                                                                ifelse(category == "Gymnosomata",0.1913*volESD,
                                                                                                       ifelse(category == "Cavoliniidae",0.1913*volESD,
                                                                                                              ifelse(category == "Limacinidae",0.1913*volESD,
                                                                                                                     ifelse(category == "Creseidae",0.1913*volESD,
                                                                                                                            ifelse(category == "Amphipoda",0.034*volESD,
                                                                                                                                   ifelse(category == "Euphausiacea",0.027*volESD,
                                                                                                                                          ifelse(category == "Decapoda",0.034*volESD,
                                                                                                                                                 0.055*volESD)))))))))))))))

community_small$obj_conc <- community_small$acq_sub_part/community_small$sample_tot_vol # Density factor d for each organism
community_small$vol_norm <- community_small$volume*community_small$obj_conc # Normalized volume of each organism, in mm3/m3

## Remove organisms based on their size
community_500um <- community_small %>% filter(major_mm >= 0.5) # 3197 organisms (including 624 copepods) removed
community_500um <- community_500um %>% filter(ifelse(community_500um$category == 'Actinopterygii', major_mm <= 20, major_mm >= 0.5)) # Remove fish > 2 cm or < 0.5 cm long)

## Group by samples and category and compute the number of images per cruise
community_500um <- community_500um %>% 
  group_by(cruise, sample_time) %>% 
  mutate(nimages = n()) 

################################################################################ SAVE DATAFRAMES
write_csv(community_500um, "community_500um.csv", na = "NA", col_names = TRUE) 

