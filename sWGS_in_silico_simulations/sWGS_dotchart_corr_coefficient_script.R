###this is a script for merging all of the the corr tables for each coverage from previous script and then plotting a Correlation coefficient dotchart based on the
###merged correlation coefficients data 

library(tidyverse)

#step 1-load in correlation numbers for each coverage
corrtable_2x_flo <- read_tsv("r_output/final_corr_table/final_corr_numbers_2x_flo.txt")
corrtable_1x_flo <- read_tsv("r_output/final_corr_table/final_corr_numbers_1x_flo.txt")
corrtable_0.5x_flo <- read_tsv("r_output/final_corr_table/final_corr_numbers_0.5x_flo.txt")
corrtable_0.1x_flo <- read_tsv("r_output/final_corr_table/final_corr_numbers_0.1x_flo.txt")
corrtable_0.01x_flo <- read_tsv("r_output/final_corr_table/final_corr_numbers_0.01x_flo.txt")
corrtable_100cell_allx_flo <- read_tsv("r_output/final_corr_table/final_corr_numbers_100cell_allx_flo.txt")

#step 2 merge the correlation tables vertically using rbind, do not sort by freq hence why arrange is hastagged 
merged_corr_table <- rbind(corrtable_2x_flo, corrtable_1x_flo, corrtable_0.5x_flo, corrtable_0.1x_flo, corrtable_0.01x_flo, corrtable_100cell_allx_flo ) %>%
  filter(!duplicated(paste0(pmax(as.character(Var1), 
                                 as.character(Var2)), 
                            pmin(as.character(Var1), as.character(Var2))))) 

merged_corr_table%>% #removed duplicated row 
  select(Var1, Freq) 


#the dataa fram row name is automatically "1", "2", "3", rename row name to fit with file names  because row name is needed for the dot chart graph
merged_corr_table<- merged_corr_table %>%
  remove_rownames() %>%
  column_to_rownames(var = "Var1")
names(merged_corr_table) #confirm to see what your row names look like

#since file names have been converted to row name, have to introduce a new column called Var1 where the coverage information can be found

merged_corr_table$Var1 <- c("2x", "2x","2x","2x","2x","2x","2x",
                            "1x", "1x", "1x", "1x", "1x", "1x",
                            "0.5x","0.5x","0.5x","0.5x","0.5x","0.5x",
                            "0.1x","0.1x","0.1x","0.1x","0.1x","0.1x",
                            "0.01x","0.01x","0.01x","0.01x","0.01x","0.01x",
                            "100%cellularity", "100%cellularity", "100%cellularity", "100%cellularity", "100%cellularity")

grp_var1 <- factor(merged_corr_table$Var1, levels = c("2x", "1x", "0.5x", "0.1x", "0.01x", "100%cellularity")) #now have to turn the coverage into factor
my_cols <- c("#ff0000",  "#56B4E9", "#e69f00", "#009e73", "#ff00ee", "#ff0000") #choose the color palletes that you would like to use

#the dot chart syntax
dotchart(merged_corr_table$Freq, labels = row.names(merged_corr_table),
         groups = grp_var1, gcolor = my_cols,
         color = my_cols[grp_var1],
         cex = 0.6,  pch = 19, xlab = "spearman_correlation")

write_tsv(merged_corr_table, "r_output/merged_corr_table.txt") #saved merged corr coefficients from all coverges for flo


######################################################################################################################
####oacp4c dotchart####
#step 1-load in all the cosine tables
corrtable_2x_oacp4c <- read_tsv("r_output/final_corr_table/final_corr_numbers_2x_oacp4c.txt")
corrtable_1x_oacp4c <- read_tsv("r_output/final_corr_table/final_corr_numbers_1x_oacp4c.txt")
corrtable_0.5x_oacp4c <- read_tsv("r_output/final_corr_table/final_corr_numbers_0.5x_oacp4c.txt")
corrtable_0.1x_oacp4c <- read_tsv("r_output/final_corr_table/final_corr_numbers_0.1x_oacp4c.txt")
corrtable_0.01x_oacp4c <- read_tsv("r_output/final_corr_table/final_corr_numbers_0.01x_oacp4c.txt")
corrtable_100cell_allx_oacp4c <- read_tsv("r_output/final_corr_table/final_corr_numbers_100cell_allx_oacp4c.txt")

#step 2 merge the cosine tables vertically using rbind, do not sort by freq hence why arrange is hastagged 
merged_corr_table <- rbind(corrtable_2x_oacp4c, corrtable_1x_oacp4c, corrtable_0.5x_oacp4c, corrtable_0.1x_oacp4c, corrtable_0.01x_oacp4c, corrtable_100cell_allx_oacp4c ) %>%
  filter(!duplicated(paste0(pmax(as.character(Var1), 
                                 as.character(Var2)), 
                            pmin(as.character(Var1), as.character(Var2))))) 

merged_corr_table%>% #removed duplicated row 
  select(Var1, Freq) 


#the data fram row name is automatically "1", "2", "3", rename row name to fit with file names  because rown name is needed for 
#the ggplot graph
merged_corr_table<- merged_corr_table %>%
  remove_rownames() %>%
  column_to_rownames(var = "Var1")
names(merged_corr_table) #confirm to see what your row names look like

#since file names have been converted to row name, have to introduce a new column called Var1 where the coverage information can be found

merged_corr_table$Var1 <- c("2x", "2x","2x","2x","2x","2x","2x",
                            "1x", "1x", "1x", "1x", "1x", "1x",
                            "0.5x","0.5x","0.5x","0.5x","0.5x","0.5x",
                            "0.1x","0.1x","0.1x","0.1x","0.1x","0.1x",
                            "0.01x","0.01x","0.01x","0.01x","0.01x","0.01x",
                            "100%cellularity", "100%cellularity", "100%cellularity", "100%cellularity", "100%cellularity")

grp_var1 <- factor(merged_corr_table$Var1, levels = c("2x", "1x", "0.5x", "0.1x", "0.01x", "100%cellularity")) #now have to turn the coverage into factor
my_cols <- c("#ff0000",  "#56B4E9", "#e69f00", "#009e73", "#ff00ee", "#ff0000") #choose the color palletes that you would like to use

#the dot chart syntax
dotchart(merged_corr_table$Freq, labels = row.names(merged_corr_table),
         groups = grp_var1, gcolor = my_cols,
         color = my_cols[grp_var1],
         cex = 0.6,  pch = 19, xlab = "spearman_correlation")

write_tsv(merged_corr_table, "r_output/merged_corr_table_oacp4c.txt") #saved merged corr coefficients from all coverges for oacp4c

