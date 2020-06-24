#this script was used for merging all of the CNA dataframes whcich were single outputs from qdnaseq 
#the script below was used for processing FLO data, this was repeated for oacp4c accordingly
library(tidyverse)
table_2x_v01 <- read_tsv("processing_data/raw_data/copyNumber_2x.txt")
chromosomal_arms <- read_tsv("processing_data/qdnaseq_bins_v03.txt",  col_names = FALSE) #a chromosomal arm column is required for the copy number 
#package to output the heatmap

arm_subselected <- chromosomal_arms %>% select(X4)

chrom_start_end <- table_2x_v01 %>% select(chromosome, start, end)

updated_start_end_arm <-  cbind(chrom_start_end, arm_subselected )

flo_2x_preselection <- table_2x_v01 %>% select(
                                            Merged_blood_FLO_TAK_50cell_2x, Merged_blood_FLO_TAK_10cell_2x,
                                            Merged_blood_FLO_TAK_5cell_2x, Merged_blood_FLO_TAK_1cell_2x,
                                            Merged_blood_FLO_TAK_0.5cell_2x, Merged_blood_FLO_TAK_0.1cell_2x)

flo_2x_selected <- cbind(updated_start_end_arm, flo_2x_preselection)
names(flo_2x_selected)[4] <- "arm"
colnames(flo_2x_selected)

####################1x

table_1x_v01 <- read_tsv("processing_data/raw_data/copyNumber_1x.txt")

flo_1x_preselection <- table_1x_v01 %>% select(
                                               Merged_blood_FLO_TAK_50cell_1x, Merged_blood_FLO_TAK_10cell_1x,
                                               Merged_blood_FLO_TAK_5cell_1x, Merged_blood_FLO_TAK_1cell_1x,
                                               Merged_blood_FLO_TAK_0.5cell_1x, Merged_blood_FLO_TAK_0.1cell_1x)

nrow(flo_2x_selected)
nrow(flo_1x_preselection)
flo_2x_selected_c <- cbind(flo_2x_selected, flo_1x_preselection)
colnames(flo_2x_selected_c)

####################0.5x

table_0.5x_v01 <- read_tsv("processing_data/raw_data/copyNumber_0.5x.txt")

flo_0.5x_preselection <- table_0.5x_v01 %>% select(
  Merged_blood_FLO_TAK_50cell_0.5x, Merged_blood_FLO_TAK_10cell_0.5x,
  Merged_blood_FLO_TAK_5cell_0.5x, Merged_blood_FLO_TAK_1cell_0.5x,
  Merged_blood_FLO_TAK_0.5cell_0.5x, Merged_blood_FLO_TAK_0.1cell_0.5x)

nrow(flo_2x_selected_c)
nrow(flo_1x_preselection)
flo_2x_selected_d <- cbind(flo_2x_selected_c, flo_0.5x_preselection)
colnames(flo_2x_selected_d)

####################0.1x

table_0.1x_v01 <- read_tsv("processing_data/raw_data/copyNumber_0.1x.txt")

flo_0.1x_preselection <- table_0.1x_v01 %>% select(
  Merged_blood_FLO_TAK_50cell_0.1x, Merged_blood_FLO_TAK_10cell_0.1x,
  Merged_blood_FLO_TAK_5cell_0.1x, Merged_blood_FLO_TAK_1cell_0.1x,
  Merged_blood_FLO_TAK_0.5cell_0.1x, Merged_blood_FLO_TAK_0.1cell_0.1x)

nrow(flo_2x_selected_d)
nrow(flo_0.1x_preselection)
flo_2x_selected_e <- cbind(flo_2x_selected_d, flo_0.1x_preselection)
colnames(flo_2x_selected_e)

####################0.01x

table_0.01x_v01 <- read_tsv("processing_data/raw_data/copyNumber_0.01x.txt")

flo_0.01x_preselection <- table_0.01x_v01 %>% select(
  Merged_blood_FLO_TAK_50cell_0.01x, Merged_blood_FLO_TAK_10cell_0.01x,
  Merged_blood_FLO_TAK_5cell_0.01x, Merged_blood_FLO_TAK_1cell_0.01x,
  Merged_blood_FLO_TAK_0.5cell_0.01x, Merged_blood_FLO_TAK_0.1cell_0.01x)

nrow(flo_2x_selected_e)
nrow(flo_0.01x_preselection)
flo_2x_selected_f <- cbind(flo_2x_selected_e, flo_0.01x_preselection)
colnames(flo_2x_selected_f)

####################100cell_allx
allx_100cell <- read_tsv("processing_data/merged_100cell_flo_allx_v02.txt")
nrow(allx_100cell)
nrow(flo_2x_selected_f)

flo_2x_selected_g <- cbind(flo_2x_selected_f, allx_100cell)
colnames(flo_2x_selected_g)

#pivot wider to longer conversion since the copy number package requires the input data to be in the longer format 
longer_flo_2x_selected_g <- flo_2x_selected_g %>% 
                         pivot_longer(cols = Merged_blood_FLO_TAK_50cell_2x:SLX.16220_FLO_TAK_100cell_0.01x,
                             names_to = "sample_name", 
                             values_to = "copynumber") 


write_tsv(longer_flo_2x_selected_g, "longer_flo_2x_selected_g.txt")
#i saved the dataframe with all of the data and then used that as input for the copy number package 

samples_flo_selection <- table_2x_v01 %>% select("SLX-16220_FLO_TAK_100cell_2x",
                                                 Merged_blood_FLO_TAK_50cell_2x, Merged_blood_FLO_TAK_10cell_2x,
                                                 Merged_blood_FLO_TAK_5cell_2x, Merged_blood_FLO_TAK_1cell_2x,
                                                 Merged_blood_FLO_TAK_0.5cell_2x, Merged_blood_FLO_TAK_0.1cell_2x,
                                                 Merged_blood_FLO_TAK_50cell_1x, Merged_blood_FLO_TAK_10cell_1x,
                                                 Merged_blood_FLO_TAK_5cell_1x, Merged_blood_FLO_TAK_1cell_1x,
                                                 Merged_blood_FLO_TAK_0.5cell_1x, Merged_blood_FLO_TAK_0.1cell_1x,
                                                 Merged_blood_FLO_TAK_50cell_0.5x, Merged_blood_FLO_TAK_10cell_0.5x,
                                                 Merged_blood_FLO_TAK_5cell_0.5x, Merged_blood_FLO_TAK_1cell_0.5x,
                                                 Merged_blood_FLO_TAK_0.5cell_0.5x, Merged_blood_FLO_TAK_0.1cell_0.5x,
                                                 Merged_blood_FLO_TAK_50cell_0.1x, Merged_blood_FLO_TAK_10cell_0.1x,
                                                 Merged_blood_FLO_TAK_5cell_0.1x, Merged_blood_FLO_TAK_1cell_0.1x,
                                                 Merged_blood_FLO_TAK_0.5cell_0.1x, Merged_blood_FLO_TAK_0.1cell_0.1x,
                                                 Merged_blood_FLO_TAK_50cell_0.01x, Merged_blood_FLO_TAK_10cell_0.01x,
                                                 Merged_blood_FLO_TAK_5cell_0.01x, Merged_blood_FLO_TAK_1cell_0.01x,
                                                 Merged_blood_FLO_TAK_0.5cell_0.01x, Merged_blood_FLO_TAK_0.1cell_0.01x,
                                                 "SLX-16220_FLO_TAK_100cell_downsample.5_1x","SLX-16220_FLO_TAK_100cell_downsample.25_0.5x", 
                                                 "SLX-16220_FLO_TAK_100cell_downsample.05_0.1x", "SLX-16220_FLO_TAK_100cell_downsample.005_0.01x")
