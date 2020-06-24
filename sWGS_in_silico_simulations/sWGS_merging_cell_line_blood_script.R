library(parallel)
library(tidyverse)
#this script was used for merging downsampled flo/oacp4c bam file with downsampled blood bam files at different cellularity percentages
#first I manually formated the merging table in excel to have the numbers represented as ".number"
#i then read the two tables for the blood and cell line separately 
# added in "col_types = "ccccc", only then can R read it as characters not introducing 0 to the ".number" every single time
merging_table_flo <- read_tsv("../downsamp_fraction_table_v04/flo_downsampled_fractions_merge.txt", col_types = "cccccc")
merging_table_blood <- read_tsv("../downsamp_fraction_table_v04/n_downsamp_fractions_merge.txt", col_types = "cccccc"  )


#the n which is the iteration value used for making the for loop is 
#first assigned to specific numbers to test that you function worked

n = 2 # column id starting at 2
nn = 1 # row number staring at 1
merging_table_flo[nn,n] #this is to check that what I am doing makes sense 

all_comm <- vector() #the all_com object has all the commands in it after it is run so 
#it has to be emptied before running it for another batch 
#you do this by assinging the all_com command to an empty vector "vector()" which basically empties the command
#otherwise one could create double files,

#i have here a nested for loop which means for every column in merging_table_blood, perform a for loop of row 
#the for loop of row is for every row in the column picked, perform this sets of instructions
#the sets of instructions are the disciption of the two input bam files names assigned to an object 
#and the last instruction is the description of the output file name assigned to an object
#the naming was done using paste0() which will leave no space between each string(set of words), if i wanted space between the words, I would have to specify as " "
#if i used paste() , it will automatically create spaces between the words 

for (n in 2:ncol(merging_table_blood) ) { #another way to have written this 2:5 since ncol is function for no of columns in a table, this means for every column starting from column 2,

    for (nn in 1:nrow(merging_table_blood)) { #for every row in the col picked starting from row 1, perform the following sets of instructions
    # Create the variable with the name of blood fraction (downsamples)
    name_blood<-paste0("LP6007415-DNA_A01_", merging_table_blood[nn,1], "cell_downsample", merging_table_blood[nn,n], "_", colnames(merging_table_blood)[n] ,".bam" )
    #print (name_blood) #this is to check that the name format looks good 
    #merging_table_flo[nn,1], a two dimensional table with only one thing changing which is the row characters/numbers of the first column, represents the celllularity numbers and picks the cell number into sample name
    #merging_table_flo[nn,n], a two dimensional table with two things changing both the row and column, represents d/s(downsample) fraction and coverage but  picks out d/s fraction into the sample name 
    #colnames(merging_table_flo)[n], a one dimensional table with only one thing changing which is the col characters of the first row, basically the heading, coverages 2x,1x,0.5x,0.1x
    name_flo_tak<-paste0("SLX-16220_FLO_TAK_", merging_table_flo[nn,1], "cell_downsample", merging_table_flo[nn,n], "_", colnames(merging_table_flo)[n] ,".bam" )
    
    name_OUT <- paste0("./out_flo/Merged_blood_FLO_TAK_", merging_table_blood[nn,1], "cell_", colnames(merging_table_blood)[n] ,".bam" )
    #print(name_OUT)
    #print(paste0(name_blood, " ", name_blood))
    #i added this ./out/ before my output file name because i wanted the merged files to be saved in the output folder in the same directory
    #for my output file name, just two things needed to be specified : the cellularity number, first row names (merging_table_flo[nn,1]) and the coverage,, column names/headingcolnames(merging_table_flo)[n]  
    command <- paste0("samtools merge ", #since i had assigned my input bam files names to an object, i did not have to paste the entire long name into my command,I therefore used those object names 
                                          #directly in my command syntax which made it easier for me to understand what is happening 
                      name_OUT,
                      " ",
                      name_blood,
                      " ",
                      name_flo_tak)
    
    all_comm <-c(all_comm, command) #what this does is that the commands for each merging operation is 
                                    #catenated with the all_com(which is the empty vector that I created at the beginning of the script)
                                    #i assigned this c(all_comm, command) back to all_com to update(fill up) the empty vectors with each of the commands that my for loop generated
    
  } 
  
}  

mclapply(all_comm, system, mc.cores = 4) #mclapply requires the command/object to come first, then the function and finally mc.cores to run on 4 processors at a time

#next step-check that your coverages are correct using collectwgsmetrics on picard tools