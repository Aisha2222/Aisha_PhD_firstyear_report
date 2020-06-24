library(tidyverse)
library(copynumber)

#this script was used for obtaining the heatmap from the copy number package for both FLO and OACP4C CN data across all of the cellularity simulations  


###flo_tak 
longer_flo_v01 <- read_tsv("processing_data/longer_flo_2x_selected_g.txt", na =".") 

#renamed columns to fit with CNV package column names

longer_flo_v02 <-  longer_flo_v01 %>% 
  rename(
    chrom = chromosome ,
    start.pos = start,
    end.pos = end,
    logR.mean =copynumber,
    sampleID = sample_name
  )

colnames(longer_flo_2x_selected_g)

longer_flo_v02$n.probes=10 #probecolumn is needed for cnv package

longer_flo_v03 <- longer_flo_v02 [, c(5,1,4,2,3,7,6)] #rearrange the columns to fit the copy number data frame format 

table(longer_flo_v03[,1]) #to check if i have NA in the sample id column

unique(longer_flo_v03$sampleID)


plotHeatmap(as.data.frame (longer_flo_v03), upper.lim = 1, lower.lim = -1, colors= c("dodgerblue","white","red"), sample.cex = 0.5 )


###oacp4c_neb
longer_oacp4c_v01 <- read_tsv("processing_data/longer_oacp4c_2x_selected_g.txt", na =".") #dont need conversion of data because x and y are absent from qdna

#rename columns to fit with CNV package column names

longer_oacp4c_v02 <-  longer_oacp4c_v01 %>% 
  rename(
    chrom = chromosome ,
    start.pos = start,
    end.pos = end,
    logR.mean =copynumber,
    sampleID = sample_name
  )

colnames(longer_oacp4c_v02)

longer_oacp4c_v02$n.probes=10 #probecolumn is needed for cnv package

longer_oacp4c_v03 <- longer_oacp4c_v02 [, c(5,1,4,2,3,7,6)] #rearrange the columns to fit the copy number data frame format 

table(longer_oacp4c_v03[,1]) #to check if i have NA in the sample id column

unique(longer_oacp4c_v03$sampleID)


plotHeatmap(as.data.frame (longer_oacp4c_v03), upper.lim=2, lower.lim = -1, colors= c("dodgerblue","white","red"), sample.cex = 0.5 )



