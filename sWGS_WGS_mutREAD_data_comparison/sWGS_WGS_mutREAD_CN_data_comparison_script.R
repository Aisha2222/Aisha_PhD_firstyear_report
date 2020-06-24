#this script was used for comparing CN data from WGS-FREEC, sWGS-qDNAseq and mutREAD -qDNAseq by plotting a heatmap on these data sets using the copy number package
#the script below was used for FLO-1 data but this was reused for the analysis of OACP4C 

library(tidyverse)
library(copynumber)

#remove NA values in mutREAD data using sed '/\NA/d' flo_cnv_50kb_bins3.txt > flo_cnv_50kb_bins3_v02.txt
#load in NA-free file 
cnv_50kb_bins3 <- read_tsv("raw_data/flo_cnv_50kb_bins3_v02.txt")
colnames(cnv_50kb_bins3)

#log2 transform CN values for mutREAD to be able to compare with sWGS and WGS CN data (which are also log2 transformed)
cnv_50kb_bins3_v02 <- cnv_50kb_bins3 %>% select(chromosome, start, end, copynumber)
mutREAD <- log2(cnv_50kb_bins3_v02$copynumber +0.01) #added 0.01 before log2 transformation because there were 0's in the CN column 
                                                    #which gave errors after log2 transformation

#merge the old non-log2 mutREAD data with the new log2 mutread data
cnv_50kb_bins3_v03 <- cbind(cnv_50kb_bins3_v02,mutREAD)

#select the relevant columns needed for downstream analysis 
cnv_50kb_bins3_v04 <- cnv_50kb_bins3_v03 %>% select(chromosome, start, end, mutREAD)


#read in old table with merged sWGS data and FREEC data
old_qdna_freec_table <- read_tsv("raw_data/wider_qdna_freec_flo.txt")

colnames(old_qdna_freec_table) #to confirm that the columns of interest are persent 

#joining of the two tables (sWGS+WGS (denoted as old_qdna_freec_table) and mutREAD( denoted as cnv_50kb_bins3_v04) )

qdna_freec_mutread <- cnv_50kb_bins3_v04 %>% inner_join(old_qdna_freec_table, by=c("chromosome","start", "end"))
colnames(qdna_freec_mutread)

#select relevant columns from merged data (qdna_freec_mutread) 
qdna_freec_mutread_v02<- qdna_freec_mutread  %>% select (chromosome, start, end, arm, "50x",flo_100cell_2x, mutREAD) #note that "50x" represent WGS data

colnames(qdna_freec_mutread_v02) #to confirm that the columns of interest are persent 
write_tsv(qdna_freec_mutread_v02, "r_output/wider_qdna_freec_mutread.txt") #save the newly merged wide table with CNs fro all three methods which 
                                                                          #would be used as input for the CN package heatmap


##################pivot transformation
#neede for compatibility with copy number package format 
wider_qdna_freec_mutread <- read_tsv("r_output/wider_qdna_freec_mutread.txt")
colnames(wider_qdna_freec_mutread)

longer_qdna_freec_mutread <- wider_qdna_freec_mutread %>% 
  pivot_longer(cols = "50x":mutREAD,
               names_to = "sample_name", 
               values_to = "copynumber")

write_tsv(longer_qdna_freec_mutread, "r_output/longer_qdna_freec_mutread.txt") #saved the long version of the merged table which would be used as CN package input

###############actual cnv package script

longer_qdna_freec_mutread <- read_tsv("r_output/longer_qdna_freec_mutread.txt", na =".") #load in data 

#rename columns to fit with CNV package column names
colnames(longer_qdna_freec_mutread)

longer_qdna_freec_mutread_v02 <-  longer_qdna_freec_mutread %>% 
  rename(
    chromosome = chromosome ,
    start.pos = start,
    end.pos = end,
    logR.mean =copynumber,
    sampleID = sample_name
  )

colnames(longer_qdna_freec_mutread_v02)


longer_qdna_freec_mutread_v02$n.probes=10 #probe column is needed for cnv package

longer_qdna_freec_mutread_v03 <- longer_qdna_freec_mutread_v02 [, c(5,1,4,2,3,7,6)] #rearrange the columns to fit the copy number data frame format 

table(longer_qdna_freec_mutread_v03[,1]) #to check if i have NA in the sample id column

unique(longer_qdna_freec_mutread_v03$sampleID)


plotHeatmap(as.data.frame (longer_qdna_freec_mutread_v03), upper.lim = 2, lower.lim = -2, colors= c("dodgerblue","white","red"), sample.cex = 0.6, sample.line = -0.4)
#sample.labels = FALSE)




l