library(tidyverse)
library(parallel)

#This is a script used for downsampling flo takara bam files in various fractions 
#loaded in the the files in my directoru
list.files("../bamfiles_at2x_v04//") -> bams
bams <- bams[grepl('bam$', bams)] #sub selected the files that specifically end with .bam from working directory and  saved them into the bams object

ds_fraction_table <- read_tsv("../downsamp_fraction_table_v04/flo_downsamp_fractions_coverage.txt", col_types ="ccccccccc")
#loaded in the downsampling fraction table with all of the fractions (represeneted as decimal points) that were to be used for downsampling 
bam_files <- bams[2] 
#selected the second file from the list of bam files saved in the bams object which is equivalent to the specific FLO_Takara file 

cellularity_number <- '100' 
#next, I picked the cellualrity percentage column in the ds_fraction_table that I was interested in, in this particular case, I was 
#interested in the 100% column, so I assingned the number '100' to an object 
#It was particularly important to assign the cellularity column names (which are essentially number) to a character object beacause numbers
#cannot be directly inserted into the code functions that I planned on using

cell <- "cell" #this is needed in order to ensure cell comes after the cellularity percentage in the function, cell means cellularity percentage 

#below is a function which downsamples the bam files using the downsampled fractions from the ds_fraction_table, 
#this automatically generates the downsampled bam files one after the next instead of downsampling the bam files at different fractions
#one by one 

down_sampled_bams_files <- function(bam_files, cellularity_number){ 
#the entire function would be passed and saved as an object named 'down_sampled_bams_files'
#this means that any command after the curly bracket '{' will be aplied to the arguments in the preceeding curved brackets which in this
#this case are 'bam_files' and 'cellularity_column' 
  seed_number <- '42' #needed for the samtools function for downsampling, "42" is the seeding number which samtools will use for downsampling
#gsub('_100cell_2x.bam$', '', gsub('.+bamfiles_at2x_v04//', '', bam_files)) is a code which is specifiying how I want the downsampled file 
#to be named, below is an attempt to explain what the code means
#the '.+' has to be present at the front of the string to be eliminated and the '//' has to be at the end of the string
#gsub removes anything ending with the word in quotation and replaces it with the secpnd quotation
#in the first gsub within the bracket, it replaces the first quotation with an empty space in the second quotation (more like deleting) 
#replaced with an empty space because I did not want those extra details in the name of the final file that is created 
#the sample name isnt perfect yet with the first gsub in the bracket 
#to make this accurate, I removed the .bam that ends the sample name, I did this using gsub function again the only difference 
#is that the  first '' will have bam with "$" the dollar sign after .bam means replace any text ending on .bam  as the file name with the texts 
#in the  second quotation '' which is essentially nothing
#after getting the sample_name right just the way I want the downsampled file to be named, I then passed the syntax into the sample_name object
#the object "sample_name" will be used in the function later on 

  cellularity_column<- ds_fraction_table %>% select(Coverage, cellularity_number) 
  #this picks the first column of the ds-fraction table and saves it in the cellularity column object which would be used in downstream function 
  
  for (n in 1:nrow(cellularity_column)){ 
    
    picked_fraction <- cellularity_column %>% 
      select(cellularity_number) %>% 
      slice(n) %>% 
      pull() %>% 
      as.character() %>% 
      gsub('^0', '', .) %>% 
      gsub(',', '.', .)
#this first selects the row 100 (which was the cellularity number object) of the cellularity column
#picks out the firts row with 'slice'
#converts it tocharacter such that the 0 can be removed afterwards with gsub. then assign this to an object called "picked_fraction""
#picked fraction object reflects the fraction that the downsampling should be done to 
    
    picked_coverage <- cellularity_column %>% 
      select(Coverage) %>% 
      slice(n) %>% 
      pull()

#this chunk of code first selects coverage row from theh cellularity column table, 
#picks out the first row of this column using slice and pulls it out and then assigned it to "picked_coverage"
#picked coverage object reflects the coverage that is to  downsampled to refelcted in the name of the output downsampled file
###the picked fraction and coverage were then be used in the downstream samtools downsampling syntax

#actual samtools downsampling syntax
#
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

### i then ran with mclappy to run the function in multiple cores and be able to generate multiple downsampled files simultaenously 
#rather than downsampling one by one. Essentially, running the parallel.version function increases the processing speed 

parallel.version <- mclapply(2:ncol(ds_fraction_table), function(colnum){
  
  ## gets column name using the number of columns in ds_fraction_table, starting column is 3 because 2 has 0's
  column.name <- colnames(ds_fraction_table[colnum])
  
  ## ran the function to generate the bam files
  down_sampled_bams_files(bam_files, column.name)
}, mc.cores = 4 ## runs this function on multiple cores
)