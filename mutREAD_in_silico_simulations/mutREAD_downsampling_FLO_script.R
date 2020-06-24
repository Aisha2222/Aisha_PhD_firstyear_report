#this is a scipt used for downsamplong mutREAD flo data 
#code has been extensively annotated in the script used for running sWGS downsampling

library(tidyverse)
library(parallel)

bams <- list.files("raw_data/", full.names = T)
bams <- bams[grepl('bam$', bams)]
ds_fraction_table <- read_tsv('ds_fraction_table/n_downsamp_fractions_coverage.txt') #%>% dplyr::rename(Coverage = `Coverage-for file`)

### dummy variables 
bam_files <- bams[1]  
## we have cellurity as a character first 
cellularity_number <- '50'
#50 as first column because 100cell is 0 for all coverages because cellurity numbers is with respect to OACP4C not flo (oapc4c would be spiked with flo)
cell <- "cell"
down_sampled_bams_files <- function(bam_files, cellularity_number){ 
  seed_number <- '21'
   sample_name <-gsub('_100cell_2x.bam$', '', gsub('raw_data//', '', bam_files))                           
  cellularity_column<- ds_fraction_table %>% dplyr::select(Coverage, cellularity_number) 
  for (n in 1:nrow(cellularity_column)){ 
    
    picked_fraction <- cellularity_column %>% 
      select(cellularity_number) %>% 
      slice(n) %>% 
      pull() %>% 
      as.character() %>% 
      gsub('^0', '', .) 
    
    picked_coverage <- cellularity_column %>% 
      select(Coverage) %>% 
      slice(n) %>% 
      pull()
    
    ###run samtools
    
    command <- paste0('samtools view -bs ', 
                      seed_number, 
                      picked_fraction, 
                      ' ', 
                      bam_files,   
                      ' > ',
                      getwd(),
                      '/',
                      sample_name, 
                      '_' , 
                      cellularity_number,
                      'cell', 
                      '_' , 
                      'downsample', 
                      picked_fraction,
                      '_', 
                      picked_coverage, 
                      '.bam')
    
    
    ## this does not run the command but prints it only.
    #cat(command, '\n')
    
    ## uncomment the system function to actually run samtools
    system(command)
    
  }
  
}  

### For loop version, for better understanding, runs 1 job 1 core at a time
###
## next level, run down_sampled_bams_files for every column in fraction table 
## column 1 is Coverage information column, lets skip it 
#for (colnum in 3:ncol(ds_fraction_table)) {
# 
#   ## run command for each column in ds_fraction_table
#   column.name <- colnames(ds_fraction_table[colnum])
#   #print(column.name)
# 
#   down_sampled_bams_files(bam_files, column.name)
# 
# }

### Run it in parallel

parallel.version <- mclapply(3:ncol(ds_fraction_table), function(colnum){
  
  ## get column name using the number of columns in ds_fraction_table 
  column.name <- colnames(ds_fraction_table[colnum])
  
  ## run the function to generate the bam files
  down_sampled_bams_files(bam_files, column.name)
}, mc.cores = 4 ## runs this function on multiple cores DANGER!!!!
)
