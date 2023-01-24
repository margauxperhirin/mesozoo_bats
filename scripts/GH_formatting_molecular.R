# SCRIPT MOLECULAR DATA ANALYSES 

remove(list = ls(all.names = TRUE))
set.seed(02121999)

##### LIBRAIRIES
library(tidyverse)
library(vegan)

##### ENTIRE COMMUNITY - DATASETS IMPORTATION AND CLEANING

################################################################################ IMPORTATION

newtax <- base::readLines("taxonomy_04-05.txt")
newtax <- read.table(text = gsub("[()\",\t]", ";", newtax))
abundances <- read.table("ASV_abundances.txt", header = TRUE, sep = "\t")
sequences <- read.table("ASV_sequences.txt", header = FALSE, sep = "\t")

################################################################################ FORMATTING

## Rename the columns
colnames(sequences) <- c("Representative_Sequence", "Fasta", "V3")
sequences <- sequences %>% select(-V3)

newtax2 <- separate(newtax, col = "V1", 
                    into = c("Representative_Sequence",
                             "Kingdom","V1","V1b","Subkingdom","V2","V2b","Infrakingdom","V3","V3b",
                             "Phylum","V4","V4b","Subphylum","V5","V5b","Infraphylum","V6","V6b",
                             "Superclass","V7","V7b","Class","V8","V8b","Subclass","V9","V9b","Infraclass","V10","V10b",
                             "Superorder","V11","V11b","Order","V12","V12b","Suborder","V13","V13b","Infraorder","V14","V14b",
                             "Superfamily","V15","V15b","Family","V16","V16b","Subfamily","V17","V17b",
                             "Genus","V18","V18b","Subgenus","V19","V19b","Species","V20","V20b"), 
                    sep = ";") # Create a column per taxonomic level

newtax3 <- newtax2[, - grep("V", colnames(newtax2))]
remove(newtax, newtax2)

## Keep only the cruises between March 2016 (321) and May 2017 (335) (but not 322) and remove the total column & merge the 3 files
ABDTAX <- left_join(abundances[,c(1,3:4,7:31)] , newtax3, by = "Representative_Sequence") # 477468 ASVs
ASV <- as.data.frame(left_join(ABDTAX, sequences, by = "Representative_Sequence")) # Rows = sequences, columns = abundance/taxonomy/fasta
remove(ABDTAX)

################################################################################ FILTERING

ASV_Animalia <- ASV %>% filter(str_detect(Kingdom, "Animalia")) # Filter zooplankton species, 433827 ASVs left
ASV_NoVertebrata <- ASV_Animalia %>% filter(!str_detect(Subphylum, "Vertebrata")) # Remove vertebrates (including fish), 416348 ASVs left

ASV_NoVertebrata$total <- rowSums(ASV_NoVertebrata[,2:28]) # Compute a new total per sequence/row

ASV_Animalia_Multiple <- ASV_NoVertebrata %>% filter(total > 1) # Remove 0 and singletons, 22814 sequences left
ASV <- ASV_Animalia_Multiple
remove(ASV_Animalia, ASV_Animalia_Multiple, ASV_NoVertebrata)

################################################################################ RAREFICATION

abundances <- ASV[,c(1:28)] # Subset the abundances
abundances_t <- as.data.frame(t(abundances)) # Transposition for vegan
colnames(abundances_t) <- abundances_t[1,] # Rename columns
abundances_t <- abundances_t[-1,] # Remove the 1st row (names in column names)
abundances_t[] <- lapply(abundances_t, function(x) if(is.character(x)) as.numeric(x) else x) # Conversion in numeric values

raremax <- min(rowSums(abundances_t)) # 92150 sequences per sample
abundances_rarefy_t <- as.data.frame(rrarefy(abundances_t, raremax))

abundances_rarefy <- as.data.frame(t(abundances_rarefy_t))
abundances_rarefy <- rownames_to_column(abundances_rarefy, var = "Representative_Sequence")

################################################################################ FINAL MERGING

ABDTAX <- left_join(abundances_rarefy, newtax3, by = "Representative_Sequence")
ASV_rarefy <- as.data.frame(left_join(ABDTAX, sequences, by = "Representative_Sequence")) # Rows = sequences, columns = abundance/taxonomy/fasta

write_csv(ASV_rarefy, "ASV_rarefy.csv", na = "NA", col_names = TRUE) # Backup (used for the taxonomy mostly, except for the images/metaB comparisons)
remove(ABDTAX, raremax)

################################################################################ FILTER 0.05% RA

## Compute RA and filter
abd <- ASV_rarefy[,1:28] # Get the reference sequences name and the number of reads per sample
abd$total <- rowSums(abd[,-1]) # Compute the sum per ASV
abd$relative <- 100 * abd$total / sum(abd$total) # Compute relative abundances per ASV in total
seq_005 <- abd %>% filter(relative > 0.05) # Filter ASVs with RA > 0.05%, 225 ASVs left

################################################################################ SAVE DATASET
write_csv(seq_005, "ASV_0-05.csv", na = "NA", col_names = TRUE)  

