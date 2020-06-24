#this is a scipt used for downsamplong mutREAD oacp4c data 
#code has been extensively annotated in the script used for running sWGS downsampling
library(tidyverse)
library(parallel)

bams <- list.files("raw_data/", full.names = T)
bams <- bams[grepl('bam$', bams)] #picked the files that specifically end with .bam 
ds_fraction_table <- read_tsv("ds_fraction_table/oacp4c_downsamp_fractions_coverage.txt", col_types ="ccccccccc")
bam_files <- bams[2]
#selected the second file which is equivalent to the OACP4C_NEB file 
cellularity_number <- '100'
cell <- "cell"
down_sampled_bams_files <- function(bam_files, cellularity_number){ 
  seed_number <- '21'
  sample_name <- gsub('_100cell_2x.bam$', '', gsub('raw_data//', '', bam_files)) 
  
  cellularity_column<- ds_fraction_table %>% select(Coverage, cellularity_number) 
  
  for (n in 1:nrow(cellularity_column)){ 
    
    picked_fraction <- cellularity_column %>% 
      select(cellularity_number) %>% 
      slice(n) %>% 
      pull() %>% 
      as.character() %>% 
      gsub('^0', '', .) %>% 
      gsub(',', '.', .)
    
    picked_coverage <- cellularity_column %>% 
      select(Coverage) %>% 
      slice(n) %>% 
      pull()
    
    command <- paste0('samtools view -bs ', 
                      seed_number, 
                      picked_fraction, 
                      ' ', 
                      bam_files,   
                      ' > ',
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
    
    #cat(command) #when checking if sample name is right for first column
    system(command)
    
  }
  
}  

### Run it in parallel

parallel.version <- mclapply(2:ncol(ds_fraction_table), function(colnum){
  
  ## get column name using the number of columns in ds_fraction_table 
  column.name <- colnames(ds_fraction_table[colnum])
  
  ## run the function to generate the bam files
  down_sampled_bams_files(bam_files, column.name)
}, mc.cores = 4## runs this function on multiple cores DANGER!!!!
)