#this script was used for comparing CN data from WGS-FREEC and sWGS-qDNAseq by plotting a heatmap on these data sets using the copy number package
#the script below was used for FLO-1 data but this was reused for the analysis of OACP4C 

library(tidyverse)
library(stringr)
library(copynumber)

#first, read in the raw freec  output (bam_ratio)
#remove the coloumns of interest
#Using excel, convert median ratio to log2 tranformed data and create a new column for start by (old start column -1), old start column becomes end
#you need to create a column for start because bedtools require start and end columns 

flo_freec_v01 <- read_tsv("raw_data/flo_freec.txt")
head(flo_freec_v01)
flo_freec_v02 <- flo_freec_v01  %>% select (Chromosome, Start, MedianRatio)
write_tsv(flo_freec_v02, "r_output/flo_freec_v02.txt")

#after first batch of excel modification, re import data into r and choose columns of interest and then resave the new version
Copy_of_flo_freec_v02 <- read_tsv("r_output/Copy of flo_freec_v02.txt")
flo_freec_v03 <- Copy_of_flo_freec_v02  %>% select (Chromosome, Start, End, log2medianR)
write_tsv(flo_freec_v03, "terminal_processing_bedtools//flo_freec_v03.txt")

#data is next modified using terminal- sorting (sort -V), removal of #NUM! (sed '/\#NUM!/d') , removal NA (sed '/\NA/d' ), 
#intersection with qdna bin (which has arm column neeeded for copy number package) and cuting of relevant columns
#read in modified intersected data 
intersected_freec_flo_v02<- read_tsv("terminal_processing_bedtools/intersected_flo_freec_v02.txt" , col_names =FALSE)

#remove mutliple columns for the same segment introduced after bedtools intersection 
intersected_freec_flo_v03 <- intersected_freec_flo_v02 %>%
  distinct()
write_tsv(intersected_freec_flo_v03, "terminal_processing_bedtools/intersected_freec_flo_v03.txt") #save the modified file as a new version 
#remove NA using  terminal as stated above 

#next, load flo_CN_wideformat_v02.txt which is the qdnaseq data on FLO (denoted as qdna_table) that is 
#to be compared to the FLO sWGS data from qdna seq (denoted as freec_table)

qdna_table <- read_tsv("../2020-04-27_20_30_40cell_CNVpackage/r_output/flo_CN_wideformat_v02.txt")
freec_table <- read_tsv("terminal_processing_bedtools/intersected_freec_flo_v04.txt", col_names = c("chromosome","start", "end", "arm", "50x"))

#merge the qdna and freec table together such that the CN data for each method is in the same data frame 
qdna_freec <- freec_table %>% inner_join(qdna_table, by=c("chromosome","start", "end", "arm"))
colnames(qdna_freec)

#pivot transformation is needed to convert from wide to long since the CN package is only compatible with long format data input
longer_qdna_freec<- qdna_freec %>% 
  pivot_longer(cols = '50x':flo_0.1cell_0.01x,
               names_to = "sample_name", 
               values_to = "copynumber")

write_tsv(qdna_freec, "r_output/wider_qdna_freec_flo.txt") #save the wide format file for the records

write_tsv(longer_qdna_freec, "r_output/longer_qdna_freec_flo.txt") #save the long format to be used as input for the CN package heatmap plot

############################## copy number heatmap syntax####################################################################3
longer_qdna_freec_flo <- read_tsv("r_output/longer_qdna_freec_flo.txt", na =".") #dont need conversion of chromosome coloumn to character 
                                                                                  #because x and y chromosomes are absent from qdnaseq profile
colnames(longer_qdna_freec_flo)

#rename columns to fit with CNV package column names

longer_qdna_freec_flo_v02 <-  longer_qdna_freec_flo %>% 
  rename(
    chrom = chromosome ,
    start.pos = start,
    end.pos = end,
    logR.mean =copynumber,
    sampleID = sample_name
  )
colnames(longer_qdna_freec_flo_v02)


longer_qdna_freec_flo_v02$n.probes=10 #probecolumn is needed for cnv package

longer_qdna_freec_flo_v03 <- longer_qdna_freec_flo_v02 [, c(5,1,4,2,3,7,6)] #rearrange the columns to fit the copy number data frame format
#sampleID chrom arm start.pos end.pos n.probes logR.mean BAF.mean

table(longer_qdna_freec_flo_v03[,1]) #to check if i have NA in the sample id column

unique(longer_qdna_freec_flo_v03$sampleID)


plotHeatmap(as.data.frame (longer_qdna_freec_flo_v03), upper.lim = 2, lower.lim = -2, colors= c("dodgerblue","white","red"), sample.cex = 0.6, sample.line = -0.4, sample.labels= FALSE)
#sample.labels = FALSE)
plotHeatmap(as.data.frame (longer_qdna_freec_flo_v03), upper.lim = 2, lower.lim = -2, colors= c("dodgerblue","white","red"), sample.cex = 0.6, sample.line = -0.4)

















