



library('DESeq2')
library("Hmisc")
library("gplots")

#********************separate normaliztion with esat ****************************************


source('~/Link to DORON/summary of filters/Gauss/load_and_normalize_esat_Gauss.R')
source("~/Link to DORON/summary of filters/Dorons_project/sources/help_tools.R")


#******run 222*********
csv_path_222_esat = "~/Link to DORON/summary of filters/Dorons_project/222/esat_res/"
file_names_ord_222 = c("C11","C12","B2", "B3","A2","A3","C2","C3","D1", "D2","F1","F2","rna")
condition_names_222 = c("cont","cont", "30_min","30_min","60_min","60_min","120_min","120_min","240_min","240_min","480_min","480_min","spike")

file_names_ord_222_no_spike = c("C11","C12","B2", "B3","A2","A3","C2","C3","D1", "D2","F1","F2")
condition_names_222_no_spike = c("cont","cont", "30_min","30_min","60_min","60_min","120_min","120_min","240_min","240_min","480_min","480_min")


# ESAT - deseq2
#orig filtred is a reference to remove non-coding and zero count reads
data_222_esat_orig_filtered <- load_data_esat(csv_path_222_esat, file_names_ord_222)
data_222_esat_t_normalized = normalize_data_sep_esat(csv_path_222_esat, file_names_ord_222_no_spike, condition_names_222_no_spike, 'pca_res_esat_filtered.png',data_222_esat_orig_filtered[,1:ncol(data_222_esat_orig_filtered)-1],"exp_222")



colnames(data_222_esat_orig_filtered) <- condition_names_222
colnames(data_222_esat_t_normalized) <- condition_names_222_no_spike


#******run 21 *********
#run 10 without the spike


csv_path_21_10_esat = "~/Link to DORON/summary of filters/Dorons_project/21_10/esat_res/"
file_names_ord_021_10 = c("S19","S20","S21", "S22","S23","S24","S25","S26","S27", "S28","S29","S30","S31","S32","S33", "S34","S35","S36","rna")
condition_names_21_10 = c("cont","cont", "cont","30_min","30_min","30_min","60_min","60_min","60_min","120_min","120_min","120_min","240_min","240_min","240_min","480_min","480_min","480_min","spike")

file_names_ord_21_10_no_spike = c("S19","S20","S21", "S22","S23","S24","S25","S26","S27", "S28","S29","S30","S31","S32","S33", "S34","S35","S36")
condition_names_21_10_no_spike = c("cont","cont", "cont","30_min","30_min","30_min","60_min","60_min","60_min","120_min","120_min","120_min","240_min","240_min","240_min","480_min","480_min","480_min")

# ESAT - deseq2
data_21_10_esat_orig_filtered = load_data_esat(csv_path_21_10_esat, file_names_ord_021_10)
data_21_10_esat_t_normalized  = normalize_data_sep_esat(csv_path_21_10_esat, file_names_ord_21_10_no_spike, condition_names_21_10_no_spike, 'pca_res_esat_filtered.png',data_21_10_esat_orig_filtered[,1:ncol(data_21_10_esat_orig_filtered)-1],"exp_21_10")


colnames(data_21_10_esat_orig_filtered) <- condition_names_21_10
colnames(data_21_10_esat_t_normalized) <- condition_names_21_10_no_spike


#******run 21 *********
#run 14 with the spike and S11(C2) outlier
csv_path_21_14_esat = "~/Link to DORON/summary of filters/Dorons_project/21-14/esat_res/"
file_names_ord_21_14 = c("S1","S2","S3", "S4","S5","S6","S7", "S8","S9","S10","S12","S13", "S14","S15","S16","S17","S18")
condition_names_21_14 = c("cont","cont","cont","30_min","30_min","30_min","60_min","60_min","60_min","120_min","240_min","240_min","240_min","480_min","480_min","480_min","spike")


file_names_ord_21_14_no_spike = c("S1","S2","S3", "S4","S5","S6","S7", "S8","S9","S10","S12","S13", "S14","S15","S16","S17")
condition_names_21_14_no_spike = c("cont","cont","cont","30_min","30_min","30_min","60_min","60_min","60_min","120_min","240_min","240_min","240_min","480_min","480_min","480_min")


# ESAT - deseq2
data_21_14_esat_orig_filtered = load_data_esat(csv_path_21_14_esat,file_names_ord_21_14)
data_21_14_esat_t_normalized = normalize_data_sep_esat(csv_path_21_14_esat,file_names_ord_21_14_no_spike,  condition_names_21_14_no_spike, 'pca_res_without_outliers_filtered.png',data_21_14_esat_orig_filtered[,1:ncol(data_21_14_esat_orig_filtered)-1],"exp_21_14")



colnames(data_21_14_esat_orig_filtered) <- condition_names_21_14 
colnames(data_21_14_esat_t_normalized) <- condition_names_21_14_no_spike


#******run 90*********
csv_path_90_esat = "/media/liorz/myVolume/DORON/run009/esat_res/"
file_names_ord_90 = c("S18","S19","S21","S22", "S23","S24","S25","S26","S29","S30", "S31","S32", "S34","S35")
#
file_names_ord_90_full = c("S18","S19","S20","S21","S22", "S23","S24","S25","S26","S27","S28","S29","S30", "S31","S32", "S33","S34","S35")
condition_names_90_full = c("cont","cont","cont","30_min","30_min","30_min","60_min","60_min","60_min","120_min","120_min","120_min","240_min","240_min","240_min","480_min","480_min","480_min")

#file_names_ord_90 = c("C11", "C12", "C13" ,"A1","A2","A3","B1","B2","B3","C1",  "C2",  "C3",  "D1",  "D2",  "D3",  "E1",  "E2"  ,"E3" )
condition_names_90 = c("cont","cont","30_min","30_min","30_min","60_min","60_min","60_min","120_min","240_min","240_min","240_min","480_min","480_min")


# ESAT - deseq2
#orig filtred is a reference to remove non-coding and zero count reads
data_90_esat_orig_filtered <- load_data_esat(csv_path_90_esat, file_names_ord_90)
data_90_esat_t_normalized = normalize_data_sep_esat(csv_path_90_esat, file_names_ord_90, condition_names_90, 'pca_res_esat_filtered_90.png',data_90_esat_orig_filtered,"exp_90")

data_90_full <- load_data_esat(csv_path_90_esat, file_names_ord_90_full)
data_90_norm_full = normalize_data_sep_esat(csv_path_90_esat, file_names_ord_90_full, condition_names_90_full, 'pca_res_esat_filtered_90.png',data_90_full,"exp_90")


# #******run 89*********
csv_path_89_esat = "/media/liorz/myVolume/DORON/run89/esat_res/"
file_names_ord_89_only_one_480 = c("C11", "C12", "C13","A1"  ,"A2",  "A3",  "B1",  "B2",  "B3", "C1" ,"C3",  "D1"  ,  "D3",  "E1"  )
condition_names_89_pnly_one_480 = c("cont","cont","cont","30_min","30_min","30_min","60_min","60_min","60_min","120_min","120_min","240_min","240_min","480_min")

file_names_ord_89 = c("C11", "C12", "C13","A1"  ,"A2",  "A3",  "B1",  "B2",  "B3", "C1" ,"C3",  "D1"  ,  "D3",  "E1"  ,"E2")
condition_names_89 = c("cont","cont","cont","30_min","30_min","30_min","60_min","60_min","60_min","120_min","120_min","240_min","240_min","480_min","480_min")



# ESAT - deseq2
#orig filtred is a reference to remove non-coding and zero count readss
data_89_esat_orig_filtered <- load_data_esat(csv_path_89_esat, file_names_ord_89)
data_89_esat_t_normalized = normalize_data_sep_esat(csv_path_89_esat, file_names_ord_89, condition_names_89, 'pca_res_esat_filtered_89.png',data_89_esat_orig_filtered,"exp_89")
data_89_esat_t_normalized_only_one_480 = normalize_data_sep_esat(csv_path_89_esat, file_names_ord_89, condition_names_89, 'pca_res_esat_filtered_89.png',data_89_esat_orig_filtered,"exp_89")

write.csv(data_89_esat_t_normalized,"data_89_esat_t_normalized.csv")
