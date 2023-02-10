# mesozoo_bats
[![OCE-17561058](https://img.shields.io/badge/NSF-1756105-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1756105) [![OCE-1829318](https://img.shields.io/badge/NSF-1829318-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1829318) [![OCE-1522206](https://img.shields.io/badge/NSF-1522206-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1522206) 

Scripts and raw data for the manuscript "Morphological and taxonomic diversity of mesozooplankton is an important driver of carbon export fluxes in the ocean"

Authors : Margaux Perhirin $^{1}$ \*, Hannah Gossner $^{2}$, Jessica Godfrey $^{2}$, Rodney Johnson $^{2}$, Leocadio Blanco-Bercial $^{2}$, Sakina-Dorothée Ayata $^{1}$

$^{1}$ Sorbonne Université, UMR 7159 CNRS-IRD-MNHN, LOCEAN-IPSL, Paris, France

$^{2}$ Bermuda Institute of Ocean Sciences, St. Georges, Bermuda

\* corresponding author: margaux.perhirin@locean.ipsl.fr


## How does it work ?

In this GitHub, you can find :
* Data : all the raw data files needed for the scripts to be run, other files will be produced during the analysis
    1. env_data_1000m.csv = all environmental data between surface and 1000 m deep from cruises at BATS
    2. export_raw.csv = export data measured by sediment traps for all the cruises between March 2016 and May 2017
    3. ecotaxa_export.zip.001 = compressed .tsv file containing the dataset exported from EcoTaxa project 149, with morphological features, taxonomic category and metadata associated to each mesozooplankton scanned
    4. ASV_abundances.zip = compressed .txt file containing raw metabarcoding data (number of reads per sequence and sample)
    5. ASV_sequences.zip.001 = compressed .txt file containing fasta sequences for all ASVs sampled between March 2016 and May 2017
    6. taxonomy_04-05.zip = compressed .txt file containing the taxonomy for all ASVs sampled between between March 2016 and May 2017
    7. ASV_0-05.csv = metabarcoding data, number of reads per sample and their total relative abundance, for the 225 most abundant ASVs (total relative abundance on the entire studied period at least equal to 0.05%)
    8. results_byhand.csv = taxonomy retrieved from NCBI BLASTN and associated features (length of the sequences, gaps, max score…) and pieces of information about the species distribution area (possible to find it at BATS or not). 

* Scripts : 
  * the mothur pipeline to treat metabarcoding sequences (ts.txt)
  * all the R scripts necessary to realise the analyses and which should be run in the following order 
    1. formatting_imaging.R
    2. formatting_metadata.R
    3. met_plots.R
    4. formatting_molecular.R
    4. mol_pcoa-clusters.R
    5. img_morphospace.R
    6. img_carbon.R
    7. mol_network-modules.R
    8. img_mol_comparison.R
* Function : a function used in mol_network-modules.R 

You can see the links between the given files (in black), the new ones (in grey), the figures/tables (in green) and the .R scripts (in blue) in the following scheme.

![GH_schema](https://user-images.githubusercontent.com/97614755/214528953-480ba12d-ee95-4835-9909-17d5bcf538cf.jpg)


If you have any questions, please contact the corresponding author. Funded by NSF [OCE-1829318](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1829318), [OCE-17561058](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1756105) and [DBI-1522206](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1522206),  [Simons Foundation International BIOS-SCOPE](https://scope.bios.edu/), ISCD project [FORMAL](https://iscd.sorbonne-universite.fr/research/sponsored-junior-teams/formal-2/), and ANR project [TraitZoo](https://anr.fr/Projet-ANR-22-CE02-0023).

Please cite the manuscript if using this code, partially or in its totality.
