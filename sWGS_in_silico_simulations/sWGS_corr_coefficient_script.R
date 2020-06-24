#This script was used for running  the Spearman correlation test between the reference file (100%cell at 2x) and different cellularity simulations comparison 
#100cell, 10cell etc mean 100% cellularity, 10% cellularity etc

library(tidyverse)

####read in2 xflo copy number data from qdnaseq (i.e different cellularity simulations at 2x) ####

raw_copynumber_2x <- read_tsv("..//raw_data/copyNumber_2x.txt") #read in the 2x_CN data file
names(raw_copynumber_2x)

#select relavant columns_2xflo
copynumber_2x_v01_flo <- raw_copynumber_2x %>%  select("SLX-16220_FLO_TAK_100cell_2x",
                                                       Merged_blood_FLO_TAK_50cell_2x, Merged_blood_FLO_TAK_10cell_2x,
                                                       Merged_blood_FLO_TAK_5cell_2x, Merged_blood_FLO_TAK_1cell_2x,
                                                       Merged_blood_FLO_TAK_0.5cell_2x, Merged_blood_FLO_TAK_0.1cell_2x)
names(copynumber_2x_v01_flo)


#spearman corr test at 2x_flo
v01_cor_numbers_2x_flo <- cor(copynumber_2x_v01_flo, method = c("spearman")) 

final_corr_numbers_2x_flo <-  v01_cor_numbers_2x_flo %>% 
  as.table() %>% as.data.frame() %>% filter(!duplicated(paste0(pmax(as.character(Var1), 
                                                                    as.character(Var2)), 
                                                               pmin(as.character(Var1), as.character(Var2))))) %>% # keep only unique occurrences, as.character because Var1 and Var2 are factors
  filter(Var2 == "SLX-16220_FLO_TAK_100cell_2x")%>%  #only interested in comparison between 100cell and the other cellularity simmulations
  arrange(desc(Freq))  # sort by Freq

write_tsv(final_corr_numbers_2x_flo, "r_output/final_corr_table/final_corr_numbers_2x_flo.txt") #saved the correlation coefficients for 2x



########################################################################################################################################333
####read in read in 1xflo copy number data from qdnaseq (i.e different cellularity simulations at 1x) a####
raw_copynumber_1x <- read_tsv("raw_data/copyNumber_1x.txt") 
names(raw_copynumber_1x)

#subselect the 100
#create 2x_100cell vector to be inserted in data frame for comparison with other columns

`SLX-16220_FLO_TAK_100cell_2x`<- copynumber_2x_v01_flo$`SLX-16220_FLO_TAK_100cell_2x`

copynumber_1x_v01_flo <- raw_copynumber_1x %>%  mutate( `SLX-16220_FLO_TAK_100cell_2x` = `SLX-16220_FLO_TAK_100cell_2x` ) %>% 
  select ("SLX-16220_FLO_TAK_100cell_2x",
          Merged_blood_FLO_TAK_50cell_1x, Merged_blood_FLO_TAK_10cell_1x,
          Merged_blood_FLO_TAK_5cell_1x, Merged_blood_FLO_TAK_1cell_1x,
          Merged_blood_FLO_TAK_0.5cell_1x, Merged_blood_FLO_TAK_0.1cell_1x)

names(copynumber_1x_v01_flo) #100cell_2x column added confirmed


#exploring the distrubution of the copy number calls for each column
plot(copynumber_1x_v01_flo$`SLX-16220_FLO_TAK_100cell_1x`, copynumber_1x_v01_flo$Merged_blood_FLO_TAK_50cell_1x, xlim= c(-2,2), ylim = c(-2,2))


#spearman corr test at 1x_flo
v01_cor_numbers_1x_flo <- cor(copynumber_1x_v01_flo, method = c("spearman"))

final_corr_numbers_1x_flo <-  v01_cor_numbers_1x_flo %>% 
  as.table() %>% as.data.frame() %>% filter(!duplicated(paste0(pmax(as.character(Var1), 
                                                                    as.character(Var2)), 
                                                               pmin(as.character(Var1), as.character(Var2))))) %>% # keep only unique occurrences, as.character because Var1 and Var2 are factors
  filter(Var2 == "SLX-16220_FLO_TAK_100cell_2x")%>%  #only interested in comparison between 100cell and the other merged cellularities 
  arrange(desc(Freq))  # sort by Freq

write_tsv(final_corr_numbers_1x_flo, "r_output/final_corr_table/final_corr_numbers_1x_flo.txt") #saved the correlation coefficients for 1x



########################################################################################################################################333

####read in read in 0.5xflo copy number data from qdnaseq (i.e different cellularity simulations at 0.5x) 
raw_copynumber_0.5x <- read_tsv("raw_data/copyNumber_0.5x.txt")
names(raw_copynumber_0.5x)

#create 2x_100cell vector to be inserted in data frame for comparison with other columns

`SLX-16220_FLO_TAK_100cell_2x`<- copynumber_2x_v01_flo$`SLX-16220_FLO_TAK_100cell_2x`

copynumber_0.5x_v01_flo <- raw_copynumber_0.5x %>%  mutate( `SLX-16220_FLO_TAK_100cell_2x` = `SLX-16220_FLO_TAK_100cell_2x` ) %>% 
  select ("SLX-16220_FLO_TAK_100cell_2x",
          Merged_blood_FLO_TAK_50cell_0.5x, Merged_blood_FLO_TAK_10cell_0.5x,
          Merged_blood_FLO_TAK_5cell_0.5x, Merged_blood_FLO_TAK_1cell_0.5x,
          Merged_blood_FLO_TAK_0.5cell_0.5x, Merged_blood_FLO_TAK_0.1cell_0.5x)

names(copynumber_1x_v01_flo) #100cell_2x column added confirmed


#exploring the distrubution of the copy number calls for each column
plot(copynumber_1x_v01_flo$`SLX-16220_FLO_TAK_100cell_1x`, copynumber_1x_v01_flo$Merged_blood_FLO_TAK_50cell_1x, xlim= c(-2,2), ylim = c(-2,2))


#spearman corr test at 2x_flo
v01_cor_numbers_0.5x_flo <- cor(copynumber_0.5x_v01_flo, method = c("spearman"))

final_corr_numbers_0.5x_flo <-  v01_cor_numbers_0.5x_flo %>% 
  as.table() %>% as.data.frame() %>% filter(!duplicated(paste0(pmax(as.character(Var1), 
                                                                    as.character(Var2)), 
                                                               pmin(as.character(Var1), as.character(Var2))))) %>% # keep only unique occurrences, as.character because Var1 and Var2 are factors
  filter(Var2 == "SLX-16220_FLO_TAK_100cell_2x")%>%  #only interested in comparison between 100cell and the other merged cellularities 
  arrange(desc(Freq))  # sort by Freq

write_tsv(final_corr_numbers_0.5x_flo, "r_output/final_corr_table/final_corr_numbers_0.5x_flo.txt")

####read in 0.1xflo copy number data####
raw_copynumber_0.5x <- read_tsv("raw_data/copyNumber_0.5x.txt")
names(raw_copynumber_0.5x)

#create 2x_100cell vector to be inserted in data frame for comparison with other columns

`SLX-16220_FLO_TAK_100cell_2x`<- copynumber_2x_v01_flo$`SLX-16220_FLO_TAK_100cell_2x`

copynumber_0.5x_v01_flo <- raw_copynumber_0.5x %>%  mutate( `SLX-16220_FLO_TAK_100cell_2x` = `SLX-16220_FLO_TAK_100cell_2x` ) %>% 
  select ("SLX-16220_FLO_TAK_100cell_2x",
          Merged_blood_FLO_TAK_50cell_0.5x, Merged_blood_FLO_TAK_10cell_0.5x,
          Merged_blood_FLO_TAK_5cell_0.5x, Merged_blood_FLO_TAK_1cell_0.5x,
          Merged_blood_FLO_TAK_0.5cell_0.5x, Merged_blood_FLO_TAK_0.1cell_0.5x)

names(copynumber_0.5x_v01_flo) #100cell_2x column added confirmed


#exploring the distrubution of the copy number calls for each column
plot(copynumber_1x_v01_flo$`SLX-16220_FLO_TAK_100cell_1x`, copynumber_1x_v01_flo$Merged_blood_FLO_TAK_50cell_1x, xlim= c(-2,2), ylim = c(-2,2))


#spearman corr test at 0.5x_flo
v01_cor_numbers_0.5x_flo <- cor(copynumber_0.5x_v01_flo, method = c("spearman"))

final_corr_numbers_0.5x_flo <-  v01_cor_numbers_0.5x_flo %>% 
  as.table() %>% as.data.frame() %>% filter(!duplicated(paste0(pmax(as.character(Var1), 
                                                                    as.character(Var2)), 
                                                               pmin(as.character(Var1), as.character(Var2))))) %>% # keep only unique occurrences, as.character because Var1 and Var2 are factors
  filter(Var2 == "SLX-16220_FLO_TAK_100cell_2x")%>%  #only interested in comparison between 100cell and the other merged cellularities 
  arrange(desc(Freq))  # sort by Freq

write_tsv(final_corr_numbers_0.5x_flo, "r_output/final_corr_table/final_corr_numbers_0.5x_flo.txt") #saved the correlation coefficients for 1x



########################################################################################################################################333

####read in read in 0.1xflo copy number data from qdnaseq (i.e different cellularity simulations at 0.1x) 

raw_copynumber_0.1x <- read_tsv("raw_data/copyNumber_0.1x.txt")
names(raw_copynumber_0.1x)

#create 2x_100cell vector to be inserted in data frame for comparison with other columns

`SLX-16220_FLO_TAK_100cell_2x`<- copynumber_2x_v01_flo$`SLX-16220_FLO_TAK_100cell_2x`

copynumber_0.1x_v01_flo <- raw_copynumber_0.1x %>%  mutate( `SLX-16220_FLO_TAK_100cell_2x` = `SLX-16220_FLO_TAK_100cell_2x` ) %>% 
  select ("SLX-16220_FLO_TAK_100cell_2x",
          Merged_blood_FLO_TAK_50cell_0.1x, Merged_blood_FLO_TAK_10cell_0.1x,
          Merged_blood_FLO_TAK_5cell_0.1x, Merged_blood_FLO_TAK_1cell_0.1x,
          Merged_blood_FLO_TAK_0.5cell_0.1x, Merged_blood_FLO_TAK_0.1cell_0.1x)

names(copynumber_0.1x_v01_flo) #100cell_2x column added confirmed


#exploring the distrubution of the copy number calls for each column
plot(copynumber_0.1x_v01_flo$`SLX-16220_FLO_TAK_100cell_1x`, copynumber_0.1x_v01_flo$Merged_blood_FLO_TAK_50cell_1x, xlim= c(-2,2), ylim = c(-2,2))


#spearman corr test at 0.1x_flo
v01_cor_numbers_0.1x_flo <- cor(copynumber_0.1x_v01_flo, method = c("spearman"))

final_corr_numbers_0.1x_flo <-  v01_cor_numbers_0.1x_flo %>% 
  as.table() %>% as.data.frame() %>% filter(!duplicated(paste0(pmax(as.character(Var1), 
                                                                    as.character(Var2)), 
                                                               pmin(as.character(Var1), as.character(Var2))))) %>% # keep only unique occurrences, as.character because Var1 and Var2 are factors
  filter(Var2 == "SLX-16220_FLO_TAK_100cell_2x")%>%  #only interested in comparison between 100cell and the other merged cellularities 
  arrange(desc(Freq))  # sort by Freq

write_tsv(final_corr_numbers_0.1x_flo, "r_output/final_corr_table/final_corr_numbers_0.1x_flo.txt") #saved the correlation coefficients for 0.1x



########################################################################################################################################333

####read in read in 0.01xflo copy number data from qdnaseq (i.e different cellularity simulations at 0.01x) 

raw_copynumber_0.01x <- read_tsv("raw_data/copyNumber_0.01x.txt")
names(raw_copynumber_0.01x)

#create 2x_100cell vector to be inserted in data frame for comparison with other columns

`SLX-16220_FLO_TAK_100cell_2x`<- copynumber_2x_v01_flo$`SLX-16220_FLO_TAK_100cell_2x`

copynumber_0.01x_v01_flo <- raw_copynumber_0.01x %>%  mutate( `SLX-16220_FLO_TAK_100cell_2x` = `SLX-16220_FLO_TAK_100cell_2x` ) %>% 
  select ("SLX-16220_FLO_TAK_100cell_2x",
          Merged_blood_FLO_TAK_50cell_0.01x, Merged_blood_FLO_TAK_10cell_0.01x,
          Merged_blood_FLO_TAK_5cell_0.01x, Merged_blood_FLO_TAK_1cell_0.01x,
          Merged_blood_FLO_TAK_0.5cell_0.01x, Merged_blood_FLO_TAK_0.1cell_0.01x)

names(copynumber_0.01x_v01_flo) #100cell_2x column added confirmed


#exploring the distrubution of the copy number calls for each column
plot(copynumber_0.01x_v01_flo$`SLX-16220_FLO_TAK_100cell_1x`, copynumber_0.01x_v01_flo$Merged_blood_FLO_TAK_50cell_1x, xlim= c(-2,2), ylim = c(-2,2))


#spearman corr test at 0.01x_flo
v01_cor_numbers_0.01x_flo <- cor(copynumber_0.01x_v01_flo, method = c("spearman"))

final_corr_numbers_0.01x_flo <-  v01_cor_numbers_0.01x_flo %>% 
  as.table() %>% as.data.frame() %>% filter(!duplicated(paste0(pmax(as.character(Var1), 
                                                                    as.character(Var2)), 
                                                               pmin(as.character(Var1), as.character(Var2))))) %>% # keep only unique occurrences, as.character because Var1 and Var2 are factors
  filter(Var2 == "SLX-16220_FLO_TAK_100cell_2x")%>%  #only interested in comparison between 100cell and the other merged cellularities 
  arrange(desc(Freq))  # sort by Freq

write_tsv(final_corr_numbers_0.01x_flo, "r_output/final_corr_table/final_corr_numbers_0.01x_flo.txt") #saved the correlation coefficients for 0.01x



########################################################################################################################################333
#comparing coefficeints for 100% cellularity across all coverages

####read in allcoverageflo copy number data, comparson betweeen 100cell at all coverages####
raw_copynumber_2x <- read_tsv("raw_data/copyNumber_2x.txt")
raw_copynumber_1x <- read_tsv("raw_data/copyNumber_1x.txt")
raw_copynumber_0.5x <- read_tsv("raw_data/copyNumber_0.5x.txt")
raw_copynumber_0.1x <- read_tsv("raw_data/copyNumber_0.1x.txt")
raw_copynumber_0.01x <- read_tsv("raw_data/copyNumber_0.01x.txt")

#important that you turn the vector to dataframe because that is the only way rbind can merge the vectors into columns in one dataframe
`SLX-16220_FLO_TAK_100cell_2x`<- data.frame(raw_copynumber_2x$`SLX-16220_FLO_TAK_100cell_2x`)
`SLX-16220_FLO_TAK_100cell_1x`<- data.frame(raw_copynumber_1x$"SLX-16220_FLO_TAK_100cell_downsample.5_1x")
`SLX-16220_FLO_TAK_100cell_0.5x`<- data.frame(raw_copynumber_0.5x$"SLX-16220_FLO_TAK_100cell_downsample.25_0.5x")
`SLX-16220_FLO_TAK_100cell_0.1x`<- data.frame(raw_copynumber_0.1x$"SLX-16220_FLO_TAK_100cell_downsample.05_0.1x")
`SLX-16220_FLO_TAK_100cell_0.01x`<- data.frame(raw_copynumber_0.01x$"SLX-16220_FLO_TAK_100cell_downsample.005_0.01x")

`SLX-16220_FLO_TAK_100cell_1x`<- data.frame(raw_copynumber_1x$"SLX-16220_FLO_TAK_100cell_downsample.5_1x")

merged_100cell <- cbind(`SLX-16220_FLO_TAK_100cell_2x`, `SLX-16220_FLO_TAK_100cell_1x`, `SLX-16220_FLO_TAK_100cell_0.5x`, `SLX-16220_FLO_TAK_100cell_0.1x`, `SLX-16220_FLO_TAK_100cell_0.01x`)


v01_cor_numbers_100cell_allx_flo <- cor(merged_100cell, method = c("spearman"))

final_corr_numbers_100cell_allx_flo <- v01_cor_numbers_100cell_allx_flo %>% as.table() %>% 
  as.data.frame() %>%      
  filter(!duplicated(paste0(pmax(as.character(Var1), 
                                 as.character(Var2)), 
                            pmin(as.character(Var1), as.character(Var2))))) %>%
  filter(Var2 == "raw_copynumber_2x..SLX.16220_FLO_TAK_100cell_2x.")

write_tsv(final_corr_numbers_100cell_allx_flo, "r_output/final_corr_table/final_corr_numbers_100cell_allx_flo.txt") #saved the correlation coefficients for 100cell across all covergaes
