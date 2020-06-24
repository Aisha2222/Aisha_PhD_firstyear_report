#this script was used for merging mutREAD FLO and OACP4C data at different percentages (OACP4C spiked in FLO-1)
#code has been extensively annotated in the script used for running sWGS merging

library(parallel)
library(tidyverse)

merging_table_oacp4c <- read_tsv("merging_table/oacp4c_downsampled_fractions_merge.txt", col_types = "cccccc")
merging_table_flo <- read_tsv("merging_table/n_downsamp_fractions_merge.txt", col_types = "cccccc"  )


n = 2 # column id starting at 2
nn = 1 # row number staring at 1
merging_table_oacp4c[nn,n] #this is to check that what I am doing makes sense 

all_comm <- vector() 
for (n in 2:ncol(merging_table_flo) ) { 
    for (nn in 1:nrow(merging_table_flo)) { 
    name_flo<-paste0("FLO_MUTREAD_", merging_table_flo[nn,1], "cell_downsample", merging_table_flo[nn,n], "_", colnames(merging_table_flo)[n] ,".bam" )
    name_oacp4c_mutread<-paste0("OACP4C_MUTREAD_", merging_table_oacp4c[nn,1], "cell_downsample", merging_table_oacp4c[nn,n], "_", colnames(merging_table_oacp4c)[n] ,".bam" )
    
    name_OUT <- paste0("./out_oacp4c/Merged_FLO_OACP4C_mutread_", merging_table_flo[nn,1], "cell_", colnames(merging_table_flo)[n] ,".bam" )
    command <- paste0("samtools merge ", 
                      name_OUT,
                      " ",
                      name_flo,
                      " ",
                      name_oacp4c_mutread)

    all_comm <-c(all_comm, command) 
  } 
  
}  

mclapply(all_comm, system, mc.cores = 4) 

#next step-check that your coverages are correct using 
#java -jar picard.jar CollectWgsMetrics I=OACP4C_MUTREAD_100cell_downsample.005_0.01x.bam 
#O=OACP4C_MUTREAD_100cell_downsample.005_0.01x_metric_file R=homo_sapiens.fa









