#this script was used for computing corr coefficients for all mutREAD cellularity simulations 
#code has been extensively annotated in the script used for running sWGS corr coefficients calculation

library(tidyverse)
#pivot tranform log2 mutREAD CN data into wider format becasue correlation test can only work for that format

longer_non_log_allcov <- read_tsv("r_output/longer_all_merged_cov_table_woflo_v07_arm.txt")
colnames(longer_non_log_allcov)

wider_log_allcov <- longer_non_log_allcov %>% pivot_wider(names_from = "sample_name", 
                                              
                                                   values_from = "log_copynumber")

write_tsv(wider_log_allcov, "r_output/qdna_Seq/wider_log_allcov_table.txt")


### spearman correlation coefficient computation 
all_cov_table <- read_tsv("r_output/qdna_Seq/wider_log_allcov_table.txt")
colnames(all_cov_table)

selected_all_cov_table <- all_cov_table %>% select (! c(chromosome, start, end, arm) )
colnames(selected_all_cov_table)

corr_numbers_allcov_table <- cor(selected_all_cov_table, method = c("spearman"))

final_corr_numbers_allcov_table<-  corr_numbers_allcov_table %>% 
  as.table() %>% as.data.frame() %>% filter(!duplicated(paste0(pmax(as.character(Var1), 
                                                                    as.character(Var2)), 
                                                               pmin(as.character(Var1), as.character(Var2))))) %>% # keep only unique occurrences, as.character because Var1 and Var2 are factors
  filter(Var2 == "mut_oac_100cell_2x")%>%  #only interested in comparison between 100cell and the other merged cellularities 
  arrange(desc(Freq))  # sort by Freq

write_tsv(final_corr_numbers_allcov_table, "r_output/qdna_Seq/final_corr_numbers_allcov_table.txt") #saved the corr coefficient table to be used for dotchart




#######for 100cell data comparing coverage

#pivot tranform log data into wider format becasue correlation test can only work for that format-100cell

longer_non_log_100cell <- read_tsv("r_output/100cell/longer_100cell_v02.txt")
colnames(longer_non_log_100cell)

wider_log_100cell <- longer_non_log_100cell %>% pivot_wider(names_from = "sample_name", 
                                                          
                                                          values_from = "log_copynumber")

write_tsv(wider_log_100cell, "r_output/dotch/wider_log_100cell.txt") #saved the wide format of table which would be used for downstream analysis 

### spearman correlation coefficient computation 
all_cov_100cell <- read_tsv("r_output/qdna_Seq/wider_log_100cell.txt")
colnames(all_cov_100cell)

selected_all_cov_100cell <- all_cov_100cell %>% select (! c(chromosome, start, end, arm) )
colnames(selected_all_cov_100cell)


corr_numbers_all_cov_100cell <- cor(selected_all_cov_100cell, method = c("spearman"))

final_corr_numbers_all_cov_100cell<-  corr_numbers_all_cov_100cell %>% 
  as.table() %>% as.data.frame() %>% filter(!duplicated(paste0(pmax(as.character(Var1), 
                                                                    as.character(Var2)), 
                                                               pmin(as.character(Var1), as.character(Var2))))) %>% # keep only unique occurrences, as.character because Var1 and Var2 are factors
  filter(Var2 == "mut_oac_100cell_2x")%>%  #only interested in comparison between 100cell and the other merged cellularities 
  arrange(desc(Freq))  # sort by Freq

write_tsv(final_corr_numbers_all_cov_100cell, "r_output/qdna_Seq/final_corr_numbers_100cell_allcov_table.txt") #saved the corr coefficient table to be used for dotchart
