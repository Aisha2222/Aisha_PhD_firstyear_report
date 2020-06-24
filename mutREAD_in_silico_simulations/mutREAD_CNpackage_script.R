library(tidyverse)
library(copynumber)

#this script was used for processing and ploting the heatmap for mutREAD CN data 


#after adding arm using bedtools, read in the output file to select relevant coluumns
#bedtools intersect -a all_merged_cov_table_woflo_v04.txt -b cytoBand_armLevel_v04.txt -loj >all_merged_cov_table_woflo_v05.txt #syntax for bedtools

#needed to rename columns as column had to be removed when using bedtools to intersect
all_merged_cov_table_woflo_v05 <- read_tsv("r_output/all_merged_cov_table_woflo_v05.txt", col_names = c("chromosome",	"start",	"end", 	"mut_oac_100cell_2x",	"mut_oac_50cell_2x",	
                                                                                                        "mut_oac_10cell_2x",	"mut_oac_5cell_2x",	"mut_oac_1cell_2x",	"mut_oac_0.5cell_2x",	
                                                                                                        "mut_oac_0.1cell_2x",	"mut_oac_100cell_1x",	"mut_oac_50cell_1x",	"mut_oac_10cell_1x",	
                                                                                                        "mut_oac_5cell_1x",	"mut_oac_1cell_1x",	"mut_oac_0.5cell_1x",	"mut_oac_0.1cell_1x",	
                                                                                                        "mut_oac_100cell_0.5x",	"mut_oac_50cell_0.5x",	"mut_oac_10cell_0.5x",	"mut_oac_5cell_0.5x",
                                                                                                        "mut_oac_1cell_0.5x",	"mut_oac_0.5cell_0.5x",	"mut_oac_0.1cell_0.5x",	"mut_oac_100cell_0.1x",	
                                                                                                        "mut_oac_50cell_0.1x",	"mut_oac_10cell_0.1x",	"mut_oac_5cell_0.1x",	"mut_oac_1cell_0.1x",	
                                                                                                        "mut_oac_0.5cell_0.1x",	"mut_oac_0.1cell_0.1x",	"mut_oac_100cell_0.01x",	"mut_oac_50cell_0.01x",
                                                                                                        "mut_oac_10cell_0.01x",	"mut_oac_5cell_0.01x",	"mut_oac_1cell_0.01x",	"mut_oac_0.5cell_0.01x",	
                                                                                                        "mut_oac_0.1cell_0.01x", "chrom_old","start_old", "end_old","arm" )) 
#select the coloumns of interest
all_merged_cov_table_woflo_v06 <- all_merged_cov_table_woflo_v05 %>% select( "chromosome",	"start",	"end", "arm","mut_oac_100cell_2x",	"mut_oac_50cell_2x",	
                                                                             "mut_oac_10cell_2x",	"mut_oac_5cell_2x",	"mut_oac_1cell_2x",	"mut_oac_0.5cell_2x",	
                                                                             "mut_oac_0.1cell_2x",	"mut_oac_100cell_1x",	"mut_oac_50cell_1x",	"mut_oac_10cell_1x",	
                                                                             "mut_oac_5cell_1x",	"mut_oac_1cell_1x",	"mut_oac_0.5cell_1x",	"mut_oac_0.1cell_1x",	
                                                                             "mut_oac_100cell_0.5x",	"mut_oac_50cell_0.5x",	"mut_oac_10cell_0.5x",	"mut_oac_5cell_0.5x",
                                                                             "mut_oac_1cell_0.5x",	"mut_oac_0.5cell_0.5x",	"mut_oac_0.1cell_0.5x",	"mut_oac_100cell_0.1x",	
                                                                             "mut_oac_50cell_0.1x",	"mut_oac_10cell_0.1x",	"mut_oac_5cell_0.1x",	"mut_oac_1cell_0.1x",	
                                                                             "mut_oac_0.5cell_0.1x",	"mut_oac_0.1cell_0.1x",	"mut_oac_100cell_0.01x",	"mut_oac_50cell_0.01x",
                                                                             "mut_oac_10cell_0.01x",	"mut_oac_5cell_0.01x",	"mut_oac_1cell_0.01x",	"mut_oac_0.5cell_0.01x",	
                                                                             "mut_oac_0.1cell_0.01x")

write_tsv(all_merged_cov_table_woflo_v06, "r_output/wider_all_merged_cov_table_woflo_v06_arm.txt")


#pivot transformation 

wider_allcov_mutread_arm <- read_tsv("r_output/wider_all_merged_cov_table_woflo_v06_arm.txt")  #read in wider version of all cov merged table, read in v03 specifically because the second row with just letters was removed 
colnames(wider_allcov_mutread_arm)
#pivot wider to longer conversion

longer_allcov_mutread_arm<- wider_allcov_mutread_arm %>% 
  pivot_longer(cols = mut_oac_100cell_2x:mut_oac_0.1cell_0.01x,
               names_to = "sample_name", 
               values_to = "copynumber")

write_tsv(longer_allcov_mutread_arm, "r_output/longer_all_merged_cov_table_woflo_v06_arm.txt")

#log tranformation of CN values
longer_allcov_mutread <- read_tsv("r_output/longer_all_merged_cov_table_woflo_v06_arm.txt", na =".") 
log_copynumber <- log2(longer_allcov_mutread$copynumber+0.01)
longer_allcov_mutread_v02 <-  cbind(longer_allcov_mutread, log_copynumber)
longer_allcov_mutread_v03 <- longer_allcov_mutread_v02 %>% select("chromosome", "start", "end", "arm", "sample_name","log_copynumber")
write_tsv(longer_allcov_mutread_v03, "r_output/longer_all_merged_cov_table_woflo_v07_arm.txt")


#CNV package
longer_allcov_mutread <- read_tsv("r_output/longer_all_merged_cov_table_woflo_v07_arm.txt")
#rename columns to fit with CNV package column names
distinct(longer_allcov_mutread,chromosome)
longer_allcov_mutread_v02 <-  longer_allcov_mutread %>% 
  rename(
    chrom = chromosome ,
    start.pos = start,
    end.pos = end,
    logR.mean =log_copynumber,
    sampleID = sample_name
  )


longer_allcov_mutread_v02$n.probes=10 #probecolumn is needed for cnv package

longer_allcov_mutread_v03 <- longer_allcov_mutread_v02 [, c(5,1,4,2,3,7,6)] #rearrange the columns to fit the copy number data frame format 

table(longer_allcov_mutread_v03[,1]) #to check if i have NA in the sample id column

unique(longer_allcov_mutread_v03$sampleID)


plotHeatmap(as.data.frame (longer_allcov_mutread_v03), upper.lim = 2, lower.lim = -2, colors= c("dodgerblue","white","red"), sample.cex = 0.6, sample.line = -0.4)
#sample.labels = FALSE)


##############################################100cell comparison of coverage#######################################3
wider_all_merged_cov_table_woflo_v06 <- read_tsv("r_output/wider_all_merged_cov_table_woflo_v06_arm.txt")

oac_100cell <- wider_all_merged_cov_table_woflo_v06 %>% select( "chromosome",	"start",	"end", "arm","mut_oac_100cell_2x",
                                                                "mut_oac_100cell_1x", "mut_oac_100cell_0.5x", "mut_oac_100cell_0.1x",
                                                                "mut_oac_100cell_0.01x")

colnames(oac_100cell)

write_tsv(oac_100cell, "r_output/100cell/wider_100cell.txt")


#pivot transformation 

wider_100cell <- read_tsv("r_output/100cell/wider_100cell.txt")  #read in wider version of all cov merged table, read in v03 specifically because the second row with just letters was removed 
colnames(wider_100cell)
#pivot wider to longer conversion

longer_100cell<- wider_100cell %>% 
  pivot_longer(cols = mut_oac_100cell_2x:mut_oac_100cell_0.01x,
               names_to = "sample_name", 
               values_to = "copynumber")

write_tsv(longer_100cell, "r_output/100cell/longer_100cell_v01.txt")

#log tranformation of CN values
longer_allcov_mutread <- read_tsv("r_output/100cell/longer_100cell_v01.txt", na =".") #dont need conversion of data because x and y are absent from qdna
log_copynumber <- log2(longer_allcov_mutread$copynumber+0.01)
longer_allcov_mutread_v02 <-  cbind(longer_allcov_mutread, log_copynumber)
longer_allcov_mutread_v03 <- longer_allcov_mutread_v02 %>% select("chromosome", "start", "end", "arm", "sample_name","log_copynumber")
write_tsv(longer_allcov_mutread_v03, "r_output/100cell/longer_100cell_v02.txt")

#CNV package for 100cell data 
longer_allcov_mutread <- read_tsv("r_output/100cell/longer_100cell_v02.txt")
#rename columns to fit with CNV package column names
distinct(longer_allcov_mutread,chromosome)
longer_allcov_mutread_v02 <-  longer_allcov_mutread %>% 
  rename(
    chrom = chromosome ,
    start.pos = start,
    end.pos = end,
    logR.mean =log_copynumber,
    sampleID = sample_name
  )


longer_allcov_mutread_v02$n.probes=10 #probecolumn is needed for cnv package

longer_allcov_mutread_v03 <- longer_allcov_mutread_v02 [, c(5,1,4,2,3,7,6)] #rearrange the columns to fit the copy number data frame format 

table(longer_allcov_mutread_v03[,1]) #to check if i have NA in the sample id column

unique(longer_allcov_mutread_v03$sampleID)


plotHeatmap(as.data.frame (longer_allcov_mutread_v03), upper.lim = 2, lower.lim = -2, colors= c("dodgerblue","white","red"), sample.cex = 0.6, sample.line = -0.4)
#sample.labels = FALSE)
