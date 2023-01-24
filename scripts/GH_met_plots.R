# SCRIPT IMAGING DATA ANALYSES - ENTIRE COMMUNITY

remove(list = ls())
set.seed(02121999)

##### LIBRAIRIES
library(tidyverse)
library(lubridate)
library(Hmisc)
library(akima)

##### DATA
env <- read.csv2("env_data_1000m.csv", header = TRUE, sep = ";")
colnames(env)[1] <- 'Cruise'

################################################################################ Temperature
temp <- env[,c(1:8,12:13)]
temp$Date <- as.Date(temp$Date, format = "%d/%m/%Y")
temp$Sampling <- paste0(month(temp$Date, lab = TRUE, locale = 'English'),
                        sep = '-',
                        year(temp$Date)) # Cruises dates 
temp[,-c(2,11)] <- lapply(temp[,-c(2,11)], 
                          function(x) if(is.character(x) | is.integer(x)) as.numeric(x) else x)
temp <- temp %>% filter(Depth <= 500)

## Akima Interpolation
interptemp <- interp(x = temp$DeciYear, y = temp$Depth, temp$Temperature,
                     nx = 250, ny = 250)
inttempdf <- interptemp %>% interp2xyz() %>% as.data.frame() # Create a new dataframe with interpolated values
colnames(inttempdf) <- c('date','depth','temp')

tplot <- ggplot(inttempdf, aes(x = date, y = depth)) +
  geom_tile(aes(fill = temp)) +
  scale_y_reverse(expand = c(0,0)) +
  scale_fill_viridis_c(option = "H",
                       breaks = c(seq(15,30,2.5)),
                       labels = c(seq(15,30,2.5)),
                       limits = c(15,30)) +
  ylab('Depth (m)') +
  scale_x_continuous(limits = c(min(temp$DeciYear), max(temp$DeciYear)), 
                     breaks = unique(temp$DeciYear)[-2], labels = unique(temp$Sampling)) +
  theme_classic(base_size = 28) +
  theme(legend.title = element_blank(),
        legend.key.height = unit(1.25, 'cm'),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.spacing.y = unit(0, 'cm'),
        legend.margin = margin(t = 0, unit = 'cm'))

################################################################################ Fluorescence
fluo <- env[,c(1:8,12,19)]
fluo$Date <- as.Date(fluo$Date, format = "%d/%m/%Y")
fluo$Sampling <- paste0(month(fluo$Date, lab = TRUE, locale = 'English'),
                        sep = '-',
                        year(fluo$Date)) # Cruises dates 
fluo[,-c(2,11)] <- lapply(fluo[,-c(2,11)], 
                          function(x) if(is.character(x) | is.integer(x)) as.numeric(x) else x)
fluo <- fluo %>% filter(Depth <= 500)

## Akima Interpolation
interpfluo <- interp(x = fluo$DeciYear, y = fluo$Depth, fluo$Fluor,
                     nx = 250, ny = 250)
intfluodf <- interpfluo %>% interp2xyz() %>% as.data.frame() # Create a new dataframe with interpolated values
colnames(intfluodf) <- c('date','depth','fluo')

## Plot - Fluorescence
fplot <- ggplot(intfluodf, aes(x = date, y = depth)) +
  geom_tile(aes(fill = fluo)) +
  scale_y_reverse(expand = c(0,0)) +
  scale_fill_viridis_c(option = "H",
                       breaks = c(seq(0,0.25,0.05)),
                       labels = c(seq(0,0.25,0.05)),
                       limits = c(0,0.25)) +
  ylab('Depth (m)')  +
  scale_x_continuous(limits = c(min(fluo$DeciYear), max(fluo$DeciYear)), 
                     breaks = unique(fluo$DeciYear)[-2], labels = unique(fluo$Sampling)) +
  theme_classic(base_size = 28) +
  theme(legend.title = element_blank(),
        legend.key.height = unit(1.25, 'cm'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.spacing.y = unit(0, 'cm'),
        legend.margin = margin(t = 0, unit = 'cm'))

################################################################################ Both
plot_grid(tplot, NULL, fplot,
          align = 'hv', ncol = 1, 
          rel_heights = c(1,-0.33,1))
