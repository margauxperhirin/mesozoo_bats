# SCRIPT META DATA ANALYSES 

remove(list = ls(all.names = TRUE))
set.seed(02121999)

##### LIBRAIRIES
library(tidyverse)
library(lubridate)

################################################################################ IMPORTATION

metadata <- read.table("env_data_1000m.csv", header = TRUE, sep = ";")
export <- read.csv2("export_raw.csv", header = TRUE, sep = ";")

################################################################################ FORMATTING
 
colnames(metadata)[1] <- 'Cruise' # Rename the columns

## Keep only the mean per cruise, using only the first 200-m of the water column
## Except for fluorescence, for which the maximum value is taken
env <- metadata %>%
  filter(Depth < 200) %>%
  filter(Cruise != 322) %>%
  group_by(Cruise, Date) %>% 
  summarise(across(-Fluor, mean), Fluorescence = max(Fluor))
remove(metadata)

## Remove columns we don't need
remove_var <- names(env) %in% c('DeciYear', 'Day', 'Hours', 'Minutes', 'Longitude', 'Latitude', 'Depth', 
                                'Conductivity', 'MLD_dens125', 'MLD_densT2', 'PAR', 'Beam', 'Pressure')
env <- env[!remove_var]

## Rename 'Sig_theta', 'MLD_bvfrq' and 'z_par_1pcnt' by simpler names
names(env)[names(env) == 'Sig_theta'] <- 'Density'
names(env)[names(env) == 'MLD_bvfrq'] <- 'MLD'
names(env)[names(env) == 'z_par_1pcnt'] <- 'PAR'

## Create a sampling column
cruise <- c('321','321','323','323','324','324','325','325','326','326',
            '327','327','328','328','329','329','330','330','331','331',
            '332','332','333','333','334','334','335','335')
time <- c('Day','Night','Day','Night','Day','Night','Day','Night','Day','Night',
          'Day','Night','Day','Night','Day','Night','Day','Night','Day','Night',
          'Day','Night','Day','Night','Day','Night','Day','Night')
sampling <- data.frame(cruise,time)
sampling$name <- paste(sampling$cruise, sampling$time, sep = '-')
env$Cruise <- as.character(env$Cruise)
env <- left_join(sampling, env, by = c('cruise' = 'Cruise'))
colnames(env)[1:3] <- c('Cruise','Time','Name')

################################################################################ ADDITION OF EXPORT VALUES

export <- export[-2,] # Remove Cruise 322

colnames(export)[1] <- 'Date' # Rename the first column
export$Date <- unique(env$Date) # Possible because same order
export[2:10] <- lapply(export[,2:10], function(x) if(is.character(x)) as.numeric(x) else x) # Convert the variables in numeric

################################################################################ MERGING

## Merge metadata and export using 'Date'
metaenv <- left_join(env, export, by = "Date", keep = FALSE, copy = FALSE, na_matches = c("na", "never"))

# Reorder the columns
order <- c('Cruise', 'Time', 'Name', 'Date', 'Year', 'Month', 'Temperature', 'Salinity',
           'O2', 'PAR', 'Density', 'VertZone', 'MLD', 'DCM', 'Fluorescence', 'C200')
metadata <- metaenv[,order]

remove(metaenv, export, env, sampling)

################################################################################ SAVE DATAFRAME
write_csv(metadata, "metadata_sample.csv", na = "NA", col_names = TRUE)
