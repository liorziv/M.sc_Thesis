#load and preper raw data
##1
##KCl 14
data_90_esat_t_normalized <- as.matrix(load_csv_file("C:/Users/Lior/Desktop/Gauss_model_fit/sources/data_90_esat_t_normalized.csv"))
data_90_esat_t_normalized<-convert_mat_to_numeric(data_90_esat_t_normalized)
mean_by_rep_90 <- get_mean_of_repeats(data_90_esat_t_normalized)
rep_fc_90 <- calc_fold_change(data_90_esat_t_normalized,mean_by_rep_90[,1])
mean_fc_90 <- calc_fold_change(mean_by_rep_90,mean_by_rep_90[,1])
peak_place_and_value_rep_90 <- find_peak_place_and_value(mean_fc_90,rownames(mean_fc_90))


##Pic 10
data_21_10_esat_t_normalized <- as.matrix(
  load_csv_file("C:/Users/Lior/Desktop/Gauss_model_fit/sources/data_21_10_esat_t_normalized.csv"))
data_21_10_esat_t_normalized <- convert_mat_to_numeric(data_21_10_esat_t_normalized)
mean_by_rep_21_10 <- get_mean_of_repeats(data_21_10_esat_t_normalized)
var_by_rep_21_10 <- get_var_of_repeats(data_21_10_esat_t_normalized)
rep_fc_21_10 <- calc_fold_change(data_21_10_esat_t_normalized,mean_by_rep_21_10[,1])
mean_fc_21_10 <- calc_fold_change(mean_by_rep_21_10,mean_by_rep_21_10[,1])
peak_place_and_value_rep_21_10 <- find_peak_place_and_value(mean_fc_21_10,rownames(mean_fc_21_10))

##Pic 10
data_222_esat_t_normalized <- as.matrix(
  load_csv_file("C:/Users/Lior/Desktop/Gauss_model_fit/sources/data_222_esat_t_normalized.csv"))
data_222_esat_t_normalized <- convert_mat_to_numeric(data_222_esat_t_normalized)
mean_by_rep_222 <- get_mean_of_repeats(data_222_esat_t_normalized)
var_by_rep_222 <- get_var_of_repeats(data_222_esat_t_normalized)
rep_fc_222 <- calc_fold_change(data_222_esat_t_normalized,mean_by_rep_222[,1])
mean_fc_222 <- calc_fold_change(mean_by_rep_222,mean_by_rep_222[,1])
peak_place_and_value_rep_222 <- find_peak_place_and_value(mean_fc_222,rownames(mean_fc_222))


##Pic 14 1
data_21_14_esat_t_normalized <- convert_mat_to_numeric(load_csv_file(paste0(data_path,"data_21_14_esat_t_normalized.csv")))
mean_by_rep_21_14<- get_mean_of_repeats(data_21_14_esat_t_normalized)
rep_fc_21_14 <- calc_fold_change(data_21_14_esat_t_normalized,mean_by_rep_21_14[,1])
mean_fc_21_14 <- calc_fold_change(mean_by_rep_21_14,mean_by_rep_21_14[,1])
var_by_rep_21_14 <- get_var_of_repeats(data_21_14_esat_t_normalized)
peak_place_and_value_rep_21_14 <- find_peak_place_and_value(mean_fc_21_14,rownames(mean_fc_21_14))

##Pic 14 2
data_89_esat_t_normalized <- convert_mat_to_numeric(load_csv_file(paste0(data_path,"data_89_esat_t_normalized.csv")))
mean_by_rep_89<- get_mean_of_repeats(data_89_esat_t_normalized)
rep_fc_89 <- calc_fold_change(data_89_esat_t_normalized,mean_by_rep_89[,1])
mean_fc_89 <- calc_fold_change(mean_by_rep_89,mean_by_rep_89[,1])
var_by_rep_89 <- get_var_of_repeats(data_89_esat_t_normalized)
peak_place_and_value_rep_89 <- find_peak_place_and_value(mean_fc_89,rownames(mean_fc_89))


## get the GS model params and p-values and pick genes be them
options(scipen = 999)

params_222 <- as.matrix(read.csv("C:/Users/Lior/Desktop/Gauss_model_fit/sources/all_model_res_222.csv", header = FALSE))
rownames(params_222) <- params_222[,1]
colnames(params_222) <- params_222[1,]
params_222 <- params_222[c(-1),c(-1)]


non_zero_genes_222 <- names(which(apply(data_222_esat_t_normalized,1,min) >0.99))
non_zero_above_1.2_222 <- names(which(peak_place_and_value_rep_222[non_zero_genes_222,1]>=1.2))

#get real idxs
real_data_222_idx <- get_real_y_idxs(params_222,"real_y")
model_data_222_idx <- get_real_y_idxs(params_222,"real_y + 10")
#get all beneth 0.05 from the data
res<-get_under_0.05(params_222,real_data_222_idx,data_222_esat_t_normalized,12,0.1,non_zero_above_1.2_222)
real_data_222 <- do.call(rbind,res[1]) 
under_0.1_222 <- as.character(do.call(rbind,res[2]))
tot_i_r_222_up <- names(which(peak_place_and_value_rep_222[under_0.1_222,2]> 1))


params_90 <- as.matrix(read.csv("C:/Users/Lior/Desktop/Gauss_model_fit/sources/all_model_res_90.csv", header = FALSE))
rownames(params_90) <- params_90[,1]
colnames(params_90) <- params_90[1,]
params_90 <- params_90[c(-1),c(-1)]


non_zero_genes_90 <- names(which(apply(data_90_esat_t_normalized,1,min) >0))
non_zero_above_1.2_90 <- names(which(peak_place_and_value_rep_90[non_zero_genes_90,1]>=1.2))

#get real idxs
real_data_90_idx <- get_real_y_idxs(params_90,"real_y")
model_data_90_idx <- get_real_y_idxs(params_90,"real_y + 10")

#get all beneth 0.05 from the data
res<-get_under_0.05(params_90,real_data_90_idx,data_90_esat_t_normalized,12,0.1,non_zero_above_1.2_90)
real_data_90 <- do.call(rbind,res[1]) 
under_0.05_90 <- as.character(do.call(rbind,res[2]))


#load 21_10 res
params_21_10 <- as.matrix(read.csv("C:/Users/Lior/Desktop/Gauss_model_fit/sources/all_model_res_21_10.csv", header = FALSE))
rownames(params_21_10) <- params_21_10[,1]
colnames(params_21_10) <- params_21_10[1,]
params_21_10 <- params_21_10[c(-1),c(-1)]


#get real idxs
non_zero_genes_21_10 <- names(which(apply(data_21_10_esat_t_normalized,1,min) >0))
non_zero_above_1.2_21_10 <- names(which(peak_place_and_value_rep_21_10[non_zero_genes_21_10,1]>=1.2))

real_data_21_10_idx <- get_real_y_idxs(params_21_10,'real_y')
#get all beneth 0.05 and returns the fdr data
res<-get_under_0.05(params_21_10,real_data_21_10_idx,data_21_10_esat_t_normalized,12,0.1,non_zero_above_1.2_21_10)
real_data_21_10 <- do.call(rbind,res[1]) 
under_0.1_21_10 <- as.character(do.call(rbind,res[2]))

tot_i_r_21_10_up <- names(which(peak_place_and_value_rep_21_10[under_0.1_21_10,2]> 1))


#load 21_14 res
params_21_14 <- as.matrix(read.csv("C:/Users/Lior/Desktop/Gauss_model_fit/sources/all_res_21_14.csv", header = FALSE))
rownames(params_21_14) <- params_21_14[,1]
colnames(params_21_14) <- params_21_14[1,]
params_21_14 <- params_21_14[c(-1),c(-1)]
params_21_14 <- as.matrix(params_21_14[,c(2:6,1)])


#get real idxs
non_zero_genes_21_14 <- names(which(apply(data_21_14_esat_t_normalized,1,min) >0))
non_zero_above_1.2_21_14 <- names(which(peak_place_and_value_rep_21_14[non_zero_genes_21_14,1]>=1.2))

real_data_21_14_idx <- 1:nrow(params_21_14) #get_real_y_idxs(params_21_14,'real_y')

#get all beneth 0.05 and returns the fdr data
res<-get_under_0.05(params_21_14,real_data_21_14_idx,data_21_14_esat_t_normalized,12,0.1,non_zero_above_1.2_21_14)
real_data_21_14 <- do.call(rbind,res[1]) 
under_0.1_21_14 <- as.character(do.call(rbind,res[2]))
length(under_0.1_21_14)

tot_i_r_21_14_up <- names(which(peak_place_and_value_rep_21_14[under_0.1_21_14,2]> 1))


#load 21_14 res no 120
params_21_14_no_120 <- as.matrix(read.csv("C:/Users/Lior/Desktop/Gauss_model_fit/sources/all_model_res_21_14_no_120.csv", header = FALSE))
rownames(params_21_14_no_120) <- params_21_14_no_120[,1]
colnames(params_21_14_no_120) <- params_21_14_no_120[1,]
params_21_14_no_120 <- params_21_14_no_120[c(-1),c(-1)]



#get real idxs
real_data_21_14_no_120_idx <- get_real_y_idxs(params_21_14_no_120,'real_y')
#get all beneth 0.05 and returns the fdr data
res<-get_under_0.05(params_21_14_no_120,real_data_21_14_no_120_idx,data_21_14_esat_t_normalized,12,0.05,non_zero_genes_21_14)

real_data_21_14_no_120 <- do.call(rbind,res[1]) 
under_0.05_21_14_no_120 <- as.character(do.call(rbind,res[2]))
length(under_0.05_21_14_no_120)

diff_21_14_no_120 <- setdiff(under_0.05_21_14,under_0.05_21_14_no_120)


params_89 <- as.matrix(read.csv("C:/Users/Lior/Desktop/Gauss_model_fit/sources/all_model_res_89_480_2.csv", header = FALSE))
rownames(params_89) <- params_89[,1]
colnames(params_89) <- params_89[1,]
params_89 <- params_89[c(-1),c(-1)]



#get real idxs
non_zero_genes_89 <- names(which(apply(data_89_esat_t_normalized,1,min) >0))
non_zero_above_1.2_89 <- names(which(peak_place_and_value_rep_89[non_zero_genes_89,1]>=1.2))

# no need for this only real no model data
# real_data_89_idx <- get_real_y_idxs(params_89,'real_y')
#get all beneth 0.05 and returns the fdr data
real_data_89_idx <- get_real_y_idxs(params_89,'real_y')
res<-get_under_0.05(params_89,real_data_89_idx,data_89_esat_t_normalized,12,0.1,non_zero_above_1.2_89)
real_data_89 <- do.call(rbind,res[1]) 
under_0.1_89 <- as.character(do.call(rbind,res[2]))
length(under_0.1_89)


tot_i_r_89_up <- names(which(peak_place_and_value_rep_89[under_0.1_89,2]> 1))

#### union of all

exp <- paste(c('KCl 10','Picrotoxin 10','Picrotoxin 14 1','Picrotoxin 14 2','KCl 14'),'Days')
gene_pool <- Reduce(union,list(tot_i_r_89_up,tot_i_r_21_10_up,tot_i_r_222_up))
gene_pool <- setdiff(gene_pool, "Lce3c")

gene_pool_pic_14 <- Reduce(union,list(tot_i_r_89_up,tot_i_r_21_14_up))

##reduced the 14 days

p_low_90_up <- find_new_common(as.matrix(bg_90[gene_pool,]),0.05,non_zero_above_1.2_90)
p_low_90_up <- names(which(peak_place_and_value_rep_90[p_low_90_up,2]> 1))

p_low_21_10_up <- find_new_common(as.matrix(cbind(peak_place_and_value_rep_21_10[gene_pool,],
                                                  real_data_21_10[gene_pool,7],
                                                  real_data_21_10[gene_pool,7])),0.05,non_zero_above_1.2_21_10)
p_low_21_10_up <- names(which(peak_place_and_value_rep_21_10[p_low_21_10_up,2]> 1))

tmp<-p_low_222_up


p_low_89_up <- find_new_common(as.matrix(cbind(peak_place_and_value_rep_89[gene_pool,],
                                               real_data_89[gene_pool,7],
                                               real_data_89[gene_pool,7])),0.05,non_zero_above_1.2_89)
p_low_89_up <- names(which(peak_place_and_value_rep_89[p_low_89_up,2]> 1))

p_low_222_up <- find_new_common(as.matrix(cbind(peak_place_and_value_rep_222[gene_pool,],
                                                real_data_222[gene_pool,7],
                                                real_data_222[gene_pool,7])),0.15,non_zero_above_1.2_222)
p_low_222_up <- names(which(peak_place_and_value_rep_222[p_low_222_up,2]> 1))



p_low_89_only_14_up <- find_new_common(as.matrix(cbind(peak_place_and_value_rep_89[gene_pool_pic_14,],
                                               real_data_89[gene_pool_pic_14,7],
                                               real_data_89[gene_pool_pic_14,7])),0.05,non_zero_above_1.2_89)
p_low_89_only_14_up <- names(which(peak_place_and_value_rep_89[p_low_89_only_14_up,2]> 1))

p_low_21_14_only_14_up <- find_new_common(as.matrix(cbind(peak_place_and_value_rep_21_14[gene_pool_pic_14,],
                                                       real_data_21_14[gene_pool_pic_14,7],
                                                       real_data_21_14[gene_pool_pic_14,7])),0.05,non_zero_above_1.2_21_14)
p_low_21_14_only_14_up <- names(which(peak_place_and_value_rep_21_14[p_low_21_14_only_14_up,2]> 1))

###########33 only to compare pic 14


all_peak_value <- matrix(0,length(gene_pool),3)
rownames(all_peak_value) <- gene_pool
colnames(all_peak_value)<- exp[1:3]
all_peak_value[p_low_222_up,1] <- peak_place_and_value_rep_222[p_low_222_up,1]
all_peak_value[p_low_21_10_up,2] <- peak_place_and_value_rep_21_10[p_low_21_10_up,1]
all_peak_value[p_low_89_up,3] <- peak_place_and_value_rep_89[p_low_89_up,1]


all_peak_time <- all_peak_value
all_peak_time[p_low_222_up,1] <- peak_place_and_value_rep_222[p_low_222_up,2]
all_peak_time[p_low_21_10_up,2] <- peak_place_and_value_rep_21_10[p_low_21_10_up,2]
all_peak_time[p_low_89_up,3] <- peak_place_and_value_rep_89[p_low_89_up,2]
colnames(all_peak_time)<-c(paste('peak_time',exp[1:3]))



all_var_at_peak <- all_peak_value
all_var_at_peak[p_low_222_up,1] <- get_peak_std(var_by_rep_222,peak_place_and_value_rep_222,p_low_222_up)
all_var_at_peak[p_low_21_10_up,2] <- get_peak_std(var_by_rep_21_10,peak_place_and_value_rep_21_10,p_low_21_10_up)
all_var_at_peak[p_low_89_up,3] <- get_peak_std(var_by_rep_89,peak_place_and_value_rep_89,p_low_89_up)
colnames(all_var_at_peak)<-c(paste('peak_var',exp[1:3]))


## only for pic 14

all_peak_value_14 <- matrix(0,length(gene_pool_pic_14),2)
rownames(all_peak_value_14) <- gene_pool_pic_14
colnames(all_peak_value_14)<- exp[3:4]
all_peak_value_14[p_low_89_only_14_up,1] <- peak_place_and_value_rep_89[p_low_89_only_14_up,1]
all_peak_value_14[p_low_21_14_only_14_up,2] <- peak_place_and_value_rep_21_14[p_low_21_14_only_14_up,1]


all_peak_time_14 <- all_peak_value_14
all_peak_time_14[p_low_89_only_14_up,1] <- peak_place_and_value_rep_89[p_low_89_only_14_up,2]
all_peak_time_14[p_low_89_only_14_up,2] <- peak_place_and_value_rep_21_14[p_low_89_only_14_up,2]
colnames(all_peak_time_14)<-c(paste('peak_time',exp[3:4]))
