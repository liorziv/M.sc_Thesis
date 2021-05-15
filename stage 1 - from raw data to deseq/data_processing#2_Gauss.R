source("~/Link to DORON/summary of filters/Dorons_project/sources/analysis and data processing functions.R")

### fold change + mean ##############################
#222

ret_val_222_t <-get_data_fold_change_mat(data_222_esat_t_normalized,2/3,1.49)
fold_change_222_t <-do.call(rbind,ret_val_222_t[1])
mean_by_rep_222 <- do.call(rbind,ret_val_222_t[3])
all_rep_FC_222 <-  do.call(rbind,ret_val_222_t[4])

#90
ret_val_90_t <-get_data_fold_change_mat(data_90_esat_t_normalized,2/3,1.49)
fold_change_90_t <-do.call(rbind,ret_val_90_t[1])
mean_by_rep_90 <- do.call(rbind,ret_val_90_t[3])
all_rep_FC_90 <-  do.call(rbind,ret_val_90_t[4])


#21_10
ret_val_21_10_t <-get_data_fold_change_mat(data_21_10_esat_t_normalized,2/3,1.49)
fold_change_21_10_t <-do.call(rbind,ret_val_21_10_t[1])
mean_by_rep_21_10 <- do.call(rbind,ret_val_21_10_t[3])
all_rep_FC_21_10 <- do.call(rbind,ret_val_21_10_t[4])


#21_14
ret_val_21_14_t <-get_data_fold_change_mat(data_21_14_esat_t_normalized,2/3,1.49)
fold_change_21_14_t <-do.call(rbind,ret_val_21_14_t[1])
mean_by_rep_21_14 <- do.call(rbind,ret_val_21_14_t[3])
all_rep_FC_21_14 <- do.call(rbind,ret_val_21_14_t[4])


#89
ret_val_89_t <-get_data_fold_change_mat(data_89_esat_t_normalized,2/3,1.49)
fold_change_89_t <-do.call(rbind,ret_val_89_t[1])
mean_by_rep_89 <- do.call(rbind,ret_val_89_t[3])
all_rep_FC_89 <- do.call(rbind,ret_val_89_t[4])

