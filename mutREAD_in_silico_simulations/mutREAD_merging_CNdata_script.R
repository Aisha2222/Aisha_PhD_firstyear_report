#this script was used for merging mutREAD CN data for simulations within and across all coverage groups 

library(tidyverse)

#read in raw merged oacp4c flo qDNAseq  data at 2x 
oac_100cell_2x <- read_tsv ("2x/OACP4C_MUTREAD_100cell_2x.bam_cnv_50kb_bins3.txt", col_names = c("feature","chromosome","start","end","mut_oac_100cell_2x"))
oac_50cell_2x <- read_tsv("2x/Merged_FLO_OACP4C_mutread_50cell_2x.bam_cnv_50kb_bins3.txt", col_names = c("feature","chromosome","start","end","mut_oac_50cell_2x"))
oac_10cell_2x <- read_tsv("2x/Merged_FLO_OACP4C_mutread_10cell_2x.bam_cnv_50kb_bins3.txt", col_names = c("feature","chromosome","start","end","mut_oac_10cell_2x"))
oac_5cell_2x <- read_tsv("2x/Merged_FLO_OACP4C_mutread_5cell_2x.bam_cnv_50kb_bins3.txt", col_names = c("feature","chromosome","start","end","mut_oac_5cell_2x"))
oac_1cell_2x <- read_tsv("2x/Merged_FLO_OACP4C_mutread_1cell_2x.bam_cnv_50kb_bins3.txt", col_names = c("feature","chromosome","start","end","mut_oac_1cell_2x"))
oac_0.5cell_2x <- read_tsv("2x/Merged_FLO_OACP4C_mutread_0.5cell_2x.bam_cnv_50kb_bins3.txt", col_names = c("feature","chromosome","start","end","mut_oac_0.5cell_2x"))
oac_0.1cell_2x <- read_tsv("2x/Merged_FLO_OACP4C_mutread_0.1cell_2x.bam_cnv_50kb_bins3.txt", col_names = c("feature","chromosome","start","end","mut_oac_0.1cell_2x"))
flo_0cell_2x <- read_tsv("2x/FLO_MUTREAD_100cell_2x.bam_cnv_50kb_bins3.txt", col_names = c("feature","chromosome","start","end","mut_flo_0cell_2x"))

#MERGING all 2x in one table
two_cell_table <- oac_100cell_2x %>% inner_join(oac_50cell_2x, by=c("feature", "chromosome","start", "end"))
three_10cell_table <- two_cell_table %>% inner_join(oac_10cell_2x, by=c("feature", "chromosome","start", "end"))
four_5cell_table <- three_10cell_table %>% inner_join(oac_5cell_2x, by=c("feature", "chromosome","start", "end"))
five_1cell_table <- four_5cell_table %>% inner_join(oac_1cell_2x, by=c("feature", "chromosome","start", "end"))
six_0.5cell_table <- five_1cell_table %>% inner_join(oac_0.5cell_2x, by=c("feature", "chromosome","start", "end"))
seven_0.1cell_table <- six_0.5cell_table %>% inner_join(oac_0.1cell_2x, by=c("feature", "chromosome","start", "end"))
eight_0cellflo_table <- seven_0.1cell_table %>% inner_join(flo_0cell_2x, by=c("feature", "chromosome","start", "end"))

#saving the table as 2x stand-alone
merged_2x_data_wo_flo <- write_tsv( seven_0.1cell_table, "r_output/merged_2x_data_wo_flo.txt")

merged_2x_data_w_flo <- write_tsv(eight_0cellflo_table, "r_output/merged_2x_data_w_flo.txt")

                        
                                                
#read in raw merged oacp4c flo qDNAseq  data at 1x 
oac_100cell_1x <- read_tsv ("1x/OACP4C_MUTREAD_100cell_downsample.5_1x.bam_cnv_50kb_bins3.txt", col_names = c("feature","chromosome","start","end","mut_oac_100cell_1x"))
oac_50cell_1x <- read_tsv ("1x/Merged_FLO_OACP4C_mutread_50cell_1x.bam_cnv_50kb_bins3.txt", col_names = c("feature","chromosome","start","end","mut_oac_50cell_1x"))
oac_10cell_1x <- read_tsv ("1x/Merged_FLO_OACP4C_mutread_10cell_1x.bam_cnv_50kb_bins3.txt", col_names = c("feature","chromosome","start","end","mut_oac_10cell_1x"))
oac_5cell_1x <- read_tsv ("1x/Merged_FLO_OACP4C_mutread_5cell_1x.bam_cnv_50kb_bins3.txt", col_names = c("feature","chromosome","start","end","mut_oac_5cell_1x"))
oac_1cell_1x <- read_tsv ("1x/Merged_FLO_OACP4C_mutread_1cell_1x.bam_cnv_50kb_bins3.txt", col_names = c("feature","chromosome","start","end","mut_oac_1cell_1x"))
oac_0.5cell_1x <- read_tsv ("1x/Merged_FLO_OACP4C_mutread_0.5cell_1x.bam_cnv_50kb_bins3.txt", col_names = c("feature","chromosome","start","end","mut_oac_0.5cell_1x"))
oac_0.1cell_1x <- read_tsv ("1x/Merged_FLO_OACP4C_mutread_0.1cell_1x.bam_cnv_50kb_bins3.txt", col_names = c("feature","chromosome","start","end","mut_oac_0.1cell_1x"))
flo_0cell_1x <- read_tsv ("1x/FLO_MUTREAD_0cell_downsample.5_1x.bam_cnv_50kb_bins3.txt", col_names = c("feature","chromosome","start","end","mut_flo_0cell_1x"))


#MERGING all 1x in one table
two_cell_table <- oac_100cell_1x %>% inner_join(oac_50cell_1x, by=c("feature", "chromosome","start", "end"))
three_10cell_table <- two_cell_table %>% inner_join(oac_10cell_1x, by=c("feature", "chromosome","start", "end"))
four_5cell_table <- three_10cell_table %>% inner_join(oac_5cell_1x, by=c("feature", "chromosome","start", "end"))
five_1cell_table <- four_5cell_table %>% inner_join(oac_1cell_1x, by=c("feature", "chromosome","start", "end"))
six_0.5cell_table <- five_1cell_table %>% inner_join(oac_0.5cell_1x, by=c("feature", "chromosome","start", "end"))
seven_0.1cell_table <- six_0.5cell_table %>% inner_join(oac_0.1cell_1x, by=c("feature", "chromosome","start", "end"))
eight_0cellflo_table <- seven_0.1cell_table %>% inner_join(flo_0cell_1x, by=c("feature", "chromosome","start", "end"))

colnames(seven_0.1cell_table)
colnames(eight_0cellflo_table)

#saving the table as 1x stand-alone
merged_1x_data_wo_flo <- write_tsv( seven_0.1cell_table, "r_output/merged_1x_data_wo_flo.txt")

merged_1x_data_w_flo <- write_tsv(eight_0cellflo_table, "r_output/merged_1x_data_w_flo.txt")



#read in raw merged oacp4c flo qDNAseq  data at 0.5x 
oac_100cell_0.5x <- read_tsv ("0.5x/OACP4C_MUTREAD_100cell_downsample.25_0.5x.bam_cnv_50kb_bins3.txt", col_names = c("feature","chromosome","start","end","mut_oac_100cell_0.5x"))
oac_50cell_0.5x <- read_tsv ("0.5x/Merged_FLO_OACP4C_mutread_50cell_0.5x.bam_cnv_50kb_bins3.txt", col_names = c("feature","chromosome","start","end","mut_oac_50cell_0.5x"))
oac_10cell_0.5x <- read_tsv ("0.5x/Merged_FLO_OACP4C_mutread_10cell_0.5x.bam_cnv_50kb_bins3.txt",col_names = c("feature","chromosome","start","end","mut_oac_10cell_0.5x"))
oac_5cell_0.5x <- read_tsv ("0.5x/Merged_FLO_OACP4C_mutread_5cell_0.5x.bam_cnv_50kb_bins3.txt", col_names = c("feature","chromosome","start","end","mut_oac_5cell_0.5x"))
oac_1cell_0.5x <- read_tsv ("0.5x/Merged_FLO_OACP4C_mutread_1cell_0.5x.bam_cnv_50kb_bins3.txt", col_names = c("feature","chromosome","start","end","mut_oac_1cell_0.5x"))
oac_0.5cell_0.5x <- read_tsv ("0.5x/Merged_FLO_OACP4C_mutread_0.5cell_0.5x.bam_cnv_50kb_bins3.txt", col_names = c("feature","chromosome","start","end","mut_oac_0.5cell_0.5x"))
oac_0.1cell_0.5x <- read_tsv ("0.5x/Merged_FLO_OACP4C_mutread_0.1cell_0.5x.bam_cnv_50kb_bins3.txt", col_names = c("feature","chromosome","start","end","mut_oac_0.1cell_0.5x"))
flo_0cell_0.5x <- read_tsv ("0.5x/FLO_MUTREAD_0cell_downsample.25_0.5x.bam_cnv_50kb_bins3.txt", col_names = c("feature","chromosome","start","end","mut_flo_0cell_0.5x"))

#MERGING all 0.5x in one table
two_cell_table <- oac_100cell_0.5x %>% inner_join(oac_50cell_0.5x, by=c("feature", "chromosome","start", "end"))
three_10cell_table <- two_cell_table %>% inner_join(oac_10cell_0.5x, by=c("feature", "chromosome","start", "end"))
four_5cell_table <- three_10cell_table %>% inner_join(oac_5cell_0.5x, by=c("feature", "chromosome","start", "end"))
five_1cell_table <- four_5cell_table %>% inner_join(oac_1cell_0.5x, by=c("feature", "chromosome","start", "end"))
six_0.5cell_table <- five_1cell_table %>% inner_join(oac_0.5cell_0.5x, by=c("feature", "chromosome","start", "end"))
seven_0.1cell_table <- six_0.5cell_table %>% inner_join(oac_0.1cell_0.5x, by=c("feature", "chromosome","start", "end"))
eight_0cellflo_table <- seven_0.1cell_table %>% inner_join(flo_0cell_0.5x, by=c("feature", "chromosome","start", "end"))

colnames(seven_0.1cell_table)
colnames(eight_0cellflo_table)

#saving the table as 0.5x stand-alone
merged_0.5x_data_wo_flo <- write_tsv( seven_0.1cell_table, "r_output/merged_0.5x_data_wo_flo.txt")

merged_0.5x_data_w_flo <- write_tsv(eight_0cellflo_table, "r_output/merged_0.5x_data_w_flo.txt")




#read in raw merged oacp4c flo qDNAseq  data at 0.1x 
oac_100cell_0.1x <- read_tsv ("0.1x/OACP4C_MUTREAD_100cell_downsample.05_0.1x.bam_cnv_50kb_bins3.txt", col_names = c("feature","chromosome","start","end","mut_oac_100cell_0.1x"))
oac_50cell_0.1x <- read_tsv ("0.1x/Merged_FLO_OACP4C_mutread_50cell_0.1x.bam_cnv_50kb_bins3.txt", col_names = c("feature","chromosome","start","end","mut_oac_50cell_0.1x"))
oac_10cell_0.1x <- read_tsv ("0.1x/Merged_FLO_OACP4C_mutread_10cell_0.1x.bam_cnv_50kb_bins3.txt", col_names = c("feature","chromosome","start","end","mut_oac_10cell_0.1x"))
oac_5cell_0.1x <- read_tsv ("0.1x/Merged_FLO_OACP4C_mutread_5cell_0.1x.bam_cnv_50kb_bins3.txt", col_names = c("feature","chromosome","start","end","mut_oac_5cell_0.1x"))
oac_1cell_0.1x <- read_tsv ("0.1x/Merged_FLO_OACP4C_mutread_1cell_0.1x.bam_cnv_50kb_bins3.txt", col_names = c("feature","chromosome","start","end","mut_oac_1cell_0.1x"))
oac_0.5cell_0.1x <- read_tsv ("0.1x/Merged_FLO_OACP4C_mutread_0.5cell_0.1x.bam_cnv_50kb_bins3.txt",col_names = c("feature","chromosome","start","end","mut_oac_0.5cell_0.1x"))
oac_0.1cell_0.1x <- read_tsv ("0.1x/Merged_FLO_OACP4C_mutread_0.1cell_0.1x.bam_cnv_50kb_bins3.txt", col_names = c("feature","chromosome","start","end","mut_oac_0.1cell_0.1x"))
flo_0cell_0.1x <- read_tsv ("0.1x/FLO_MUTREAD_0cell_downsample.05_0.1x.bam_cnv_50kb_bins3.txt", col_names = c("feature","chromosome","start","end","mut_flo_0cell_0.1x"))


#MERGING all 0.1x in one table
two_cell_table <- oac_100cell_0.1x %>% inner_join(oac_50cell_0.1x, by=c("feature", "chromosome","start", "end"))
three_10cell_table <- two_cell_table %>% inner_join(oac_10cell_0.1x, by=c("feature", "chromosome","start", "end"))
four_5cell_table <- three_10cell_table %>% inner_join(oac_5cell_0.1x, by=c("feature", "chromosome","start", "end"))
five_1cell_table <- four_5cell_table %>% inner_join(oac_1cell_0.1x, by=c("feature", "chromosome","start", "end"))
six_0.5cell_table <- five_1cell_table %>% inner_join(oac_0.5cell_0.1x, by=c("feature", "chromosome","start", "end"))
seven_0.1cell_table <- six_0.5cell_table %>% inner_join(oac_0.1cell_0.1x, by=c("feature", "chromosome","start", "end"))
eight_0cellflo_table <- seven_0.1cell_table %>% inner_join(flo_0cell_0.1x, by=c("feature", "chromosome","start", "end"))

colnames(seven_0.1cell_table)
colnames(eight_0cellflo_table)

#saving the table as 0.5x stand-alone
merged_0.1x_data_wo_flo <- write_tsv( seven_0.1cell_table, "r_output/merged_0.1x_data_wo_flo.txt")

merged_0.1x_data_w_flo <- write_tsv(eight_0cellflo_table, "r_output/merged_0.1x_data_w_flo.txt")



#read in raw merged oacp4c flo qDNAseq  data at 0.01x 
oac_100cell_0.01x <- read_tsv ("0.01x/OACP4C_MUTREAD_100cell_downsample.005_0.01x.bam_cnv_50kb_bins3.txt", col_names = c("feature","chromosome","start","end","mut_oac_100cell_0.01x"))
oac_50cell_0.01x <- read_tsv ("0.01x/Merged_FLO_OACP4C_mutread_50cell_0.01x.bam_cnv_50kb_bins3.txt", col_names = c("feature","chromosome","start","end","mut_oac_50cell_0.01x"))
oac_10cell_0.01x <- read_tsv ("0.01x/Merged_FLO_OACP4C_mutread_10cell_0.01x.bam_cnv_50kb_bins3.txt", col_names = c("feature","chromosome","start","end","mut_oac_10cell_0.01x"))
oac_5cell_0.01x <- read_tsv ("0.01x/Merged_FLO_OACP4C_mutread_5cell_0.01x.bam_cnv_50kb_bins3.txt", col_names = c("feature","chromosome","start","end","mut_oac_5cell_0.01x"))
oac_1cell_0.01x <- read_tsv ("0.01x/Merged_FLO_OACP4C_mutread_1cell_0.01x.bam_cnv_50kb_bins3.txt", col_names = c("feature","chromosome","start","end","mut_oac_1cell_0.01x"))
oac_0.5cell_0.01x <- read_tsv ("0.01x/Merged_FLO_OACP4C_mutread_0.5cell_0.01x.bam_cnv_50kb_bins3.txt", col_names = c("feature","chromosome","start","end","mut_oac_0.5cell_0.01x"))
oac_0.1cell_0.01x <- read_tsv ("0.01x/Merged_FLO_OACP4C_mutread_0.1cell_0.01x.bam_cnv_50kb_bins3.txt", col_names = c("feature","chromosome","start","end","mut_oac_0.1cell_0.01x"))
flo_0cell_0.01x <- read_tsv ("0.01x/FLO_MUTREAD_0cell_downsample.005_0.01x.bam_cnv_50kb_bins3.txt",col_names = c("feature","chromosome","start","end","mut_flo_0cell_0.01x"))


#MERGING all 0.01x in one table
two_cell_table <- oac_100cell_0.01x %>% inner_join(oac_50cell_0.01x, by=c("feature", "chromosome","start", "end"))
three_10cell_table <- two_cell_table %>% inner_join(oac_10cell_0.01x, by=c("feature", "chromosome","start", "end"))
four_5cell_table <- three_10cell_table %>% inner_join(oac_5cell_0.01x, by=c("feature", "chromosome","start", "end"))
five_1cell_table <- four_5cell_table %>% inner_join(oac_1cell_0.01x, by=c("feature", "chromosome","start", "end"))
six_0.5cell_table <- five_1cell_table %>% inner_join(oac_0.5cell_0.01x, by=c("feature", "chromosome","start", "end"))
seven_0.1cell_table <- six_0.5cell_table %>% inner_join(oac_0.1cell_0.01x, by=c("feature", "chromosome","start", "end"))
eight_0cellflo_table <- seven_0.1cell_table %>% inner_join(flo_0cell_0.01x, by=c("feature", "chromosome","start", "end"))

colnames(seven_0.1cell_table)
colnames(eight_0cellflo_table)

#saving the table as 0.5x stand-alone
merged_0.01x_data_wo_flo <- write_tsv( seven_0.1cell_table, "r_output/merged_0.01x_data_wo_flo.txt")

merged_0.01x_data_w_flo <- write_tsv(eight_0cellflo_table, "r_output/merged_0.01x_data_w_flo.txt")


#re read data for each coverage 
cov_2x <- read_tsv("r_output/merged_2x_data_wo_flo.txt")
cov_1x <- read_tsv("r_output/merged_1x_data_wo_flo.txt")
cov_0.5x <- read_tsv("r_output/merged_0.5x_data_wo_flo.txt")
cov_0.1x <- read_tsv("r_output/merged_0.1x_data_wo_flo.txt")
cov_0.01x <- read_tsv("r_output/merged_0.01x_data_wo_flo.txt")

#merg each table 
two_cov_table <- cov_2x %>% inner_join(cov_1x, by=c("feature", "chromosome","start", "end"))
three_cov_table <- two_cov_table %>% inner_join(cov_0.5x, by=c("feature", "chromosome","start", "end"))
colnames(three_cov_table)
four_cov_table <- three_cov_table %>% inner_join(cov_0.1x, by=c("feature", "chromosome","start", "end"))
colnames(four_cov_table)
five_cov_table <- four_cov_table %>% inner_join(cov_0.01x, by=c("feature", "chromosome","start", "end"))
colnames(five_cov_table)

all_merged_cov_table <- five_cov_table %>% select(!feature)
colnames(all_merged_cov_table)

write_tsv(all_merged_cov_table, "r_output/all_merged_cov_table_woflo.txt") #saved the dataframe with all of the data 








