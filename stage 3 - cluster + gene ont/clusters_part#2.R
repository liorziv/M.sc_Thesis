library(broom)
choose_fit<-function(label_1,label_2,df){
  if(label_1 == 'KCl 14 Days' &  label_2 == 'KCl 10 Days'){
    fit1 <-  lm(`KCl 10 Days` ~0+`KCl 14 Days`, data = df)
  }
  if(label_1 == 'Picrotoxin 14 1 Days' &  label_2 == 'KCl 10 Days'){
    fit1 <-  lm(`KCl 10 Days` ~0+`Picrotoxin 14 1 Days`, data = df)
  }
  if(label_1 == 'Picrotoxin 14 2 Days' &  label_2 == 'KCl 10 Days'){
    fit1 <-  lm(`KCl 10 Days` ~0+`Picrotoxin 14 2 Days`, data = df)
  }
  if(label_1 == 'Picrotoxin 10 Days' &  label_2 == 'KCl 10 Days'){
    fit1 <-  lm(`KCl 10 Days` ~0+`Picrotoxin 10 Days`, data = df)
  }
  if(label_1 == 'KCl 10 Days' &  label_2 == 'Picrotoxin 10 Days'){
    fit1 <-  lm(`Picrotoxin 10 Days` ~0+`KCl 10 Days`, data = df)
  }
  if(label_1 == 'KCl 14 Days' &  label_2 == 'Picrotoxin 10 Days'){
    fit1 <-  lm(`Picrotoxin 10 Days` ~0+`KCl 14 Days`, data = df)
  }
  if(label_1 == 'Picrotoxin 14 1 Days' &  label_2 == 'Picrotoxin 10 Days'){
    fit1 <-  lm(`Picrotoxin 10 Days` ~0+`Picrotoxin 14 1 Days`, data = df)
  }
  if(label_1 == 'Picrotoxin 14 2 Days' &  label_2 == 'Picrotoxin 10 Days'){
    fit1 <-  lm(`Picrotoxin 10 Days` ~0+`Picrotoxin 14 2 Days`, data = df)
  }
  if(label_1 == 'KCl 10 Days' &  label_2 == 'Picrotoxin 14 1 Days'){
    fit1 <-  lm(`Picrotoxin 14 1 Days` ~0+`KCl 10 Days`, data = df)
  }
  if(label_1 == 'KCl 14 Days' &  label_2 == 'Picrotoxin 14 1 Days'){
    fit1 <-  lm(`Picrotoxin 14 1 Days` ~0+`KCl 14 Days`, data = df)
  }
  if(label_1 == 'Picrotoxin 14 2 Days' &  label_2 == 'Picrotoxin 14 1 Days'){
    fit1 <-  lm(`Picrotoxin 14 1 Days` ~0+`Picrotoxin 14 2 Days`, data = df)
  }
  if(label_1 == 'Picrotoxin 10 Days' &  label_2 == 'Picrotoxin 14 1 Days'){
    fit1 <-  lm(`Picrotoxin 14 1 Days` ~0+`Picrotoxin 10 Days`, data = df)
  }
  if(label_1 == 'KCl 10 Days' &  label_2 == 'Picrotoxin 14 2 Days'){
    fit1 <-  lm(`Picrotoxin 14 2 Days` ~0+`KCl 10 Days`, data = df)
  }
  if(label_1 == 'KCl 14 Days' &  label_2 == 'Picrotoxin 14 2 Days'){
    fit1 <-  lm(`Picrotoxin 14 2 Days` ~0+`KCl 14 Days`, data = df)
  }
  if(label_1 == 'Picrotoxin 14 1 Days' &  label_2 == 'Picrotoxin 14 2 Days'){
    fit1 <-  lm(`Picrotoxin 14 2 Days` ~0+`Picrotoxin 14 1 Days`, data = df)
  }
  if(label_1 == 'Picrotoxin 10 Days' &  label_2 == 'Picrotoxin 14 2 Days'){
    fit1 <-  lm(`Picrotoxin 14 2 Days` ~0+`Picrotoxin 10 Days`, data = df)
  }
  
  if(label_1 == 'Picrotoxin 14 2 Days' &  label_2 == 'KCl 14 Days'){
    fit1 <-  lm(`KCl 14 Days` ~0+`Picrotoxin 14 2 Days`, data = df)
  }
  if(label_1 == 'Picrotoxin 14 1 Days' &  label_2 == 'KCl 14 Days'){
    fit1 <-  lm(`KCl 14 Days` ~0+`Picrotoxin 14 1 Days`, data = df)
  }
  
  if(label_1 == 'KCl 10 Days' &  label_2 == 'KCl 14 Days'){
    fit1 <-  lm(`KCl 14 Days` ~0+`KCl 10 Days`, data = df)
  }
  if(label_1 == 'Picrotoxin 10 Days' &  label_2 == 'KCl 14 Days'){
    fit1 <-  lm(`KCl 14 Days` ~0+`Picrotoxin 10 Days`, data = df)
  }
  
  return(fit1)
}
align_genes<- function(gene_1,gene_2,peak_time){
  peak_diff <- peak_time[1]-peak_time[2]
  if(abs(peak_diff)<=3 & peak_diff != 0){
    if(peak_diff>0){
      gene_1 <- gene_1[(peak_diff+1):6]
      gene_2 <- gene_2[1:(6-peak_diff)]
    }
    if(peak_diff<0){
      peak_diff <- -peak_diff
      gene_1 <- gene_1[1:(6-peak_diff)]
      gene_2 <- gene_2[(peak_diff+1):6]
    }
  }
  return (list(gene_1,gene_2))
}
plot_hist <- function(dat){
  dat <- data.frame(x = round(as.numeric(dat),2), group = rep('sienna1',length(dat)))
  dat_count <- rbind(ddply(dat,.(x,group),summarise,count = length(group)))
  
  p <-ggplot(dat_count, aes(x, y= count,fill = group, colour=group)) + 
    geom_bar(colour = "black",position = "dodge",stat = "identity")+
    ggtitle(paste0("Fit hist pic 14 1 to pic 10: ",round(curr_center[1],2),',kcl to pic 10: ',
                   round(curr_center[2],2),',kcl 10 to pic 14 1: ',round(curr_center[3],2)))+xlab("Scaling factor")+
    ylab("# of genes")+  ylim(0,max(as.numeric(dat_count[,3])))
  print(p+ scale_x_continuous(breaks =  breaks , labels = breaks)) 
  
  dat_norm <- dat_count
  dat_norm[!is.na(dat_norm$count),3]<- dat_count[!is.na(dat_count$count),3]/sum(dat_count[!is.na(dat_count$count),3]) 
  breaks<- seq(min(dat[,1]),max(dat[,1]),10)
  # dat[,1]<- dat[,1]/sum(dat[,1])
  p <-ggplot(dat_norm, aes(x, y= count,fill = "sienna1", colour="sienna1")) + 
    geom_bar(colour = "sienna1",position = "dodge",stat = "identity",fill = "sienna1",width=1)+
    ggtitle(paste0("Fit hist pic 14 1 to pic 10: ",round(curr_center[1],2),',kcl to pic 10: ',
                   round(curr_center[2],2),',kcl 10 to pic 14 1: ',round(curr_center[3],2)))+xlab("Scaling factor")+
    ylab("# of genes")+  ylim(0,max(as.numeric(dat_norm[,3]))) 
  # +geom_vline(xintercept = 0.9, linetype="dotted", 
  #                                                                          color = "blue", size=1.5)+geom_vline(xintercept = 1.1, linetype="dotted", 
  #             color = "blue", size=1.5)
  print(p+ scale_x_continuous(breaks =  breaks , labels = breaks)) + theme(axis.text =  element_text( size = 24),axis.ticks.x =element_line(size =2),title = element_text(size = 30))
  
  
}
create_unique_mat<-function(mat){
  mat<-as.matrix(mat)
  unique_mat <- c()
  unique_names <- c()
  for (i in 1:nrow(mat)){
    if(length(which(unique_mat == mat[i])) == 0){
      unique_names<- c(unique_names,rownames(mat)[i])
      unique_mat<- c(unique,mat[i])
      
    }
  }
  
  unique_mat <- as.matrix(unique_mat)
  rownames(unique_mat)<- unique_names
  return(unique_mat)
}
adjust_p_value <- function(data_to_adjust){
  n <- length(data_to_adjust)
  print(n)
  p <- data_to_adjust 
  adjusted.p_real <- as.matrix(p.adjust(p, "BH")) #same as p.adjust(p, "BH")
  rownames(adjusted.p_real) <- names(data_to_adjust)
  
  return(adjusted.p_real)
}
gene_polyfit <-function(gene_1,gene_2){
    df<- data.frame(gene_1 = as.numeric(gene_1), gene_2 = as.numeric(gene_2))
    fit1 <- lm( gene_1 ~ 0+gene_2, data = df)
    corr<-cor.test(gene_1,gene_2)
    # round_slope <- round(fit1$coefficients[1],1)
    # r2<-calculte_r2(df,round_slope)
    # 
  return(list(corr,fit1))
}
fit_gene_list<- function(data_1,data_2,group,peak_time){
  i = 1
  fit_data <- c()
  
  for (gene in group){
    aligned_genes <- align_genes(data_1[gene,],data_2[gene,],peak_time[gene,])
    gene_1 <- aligned_genes[[1]]#data_1[gene,]
    gene_2 <- aligned_genes[[2]]#data_2[gene,]
    fit_curr_data<- gene_polyfit(as.numeric(gene_1),as.numeric(gene_2))
    fit_data <- rbind(fit_data,c(fit_curr_data[[1]]$estimate,
                      fit_curr_data[[1]]$p.value,
                      fit_curr_data[[2]]$coefficients[1],
                      summary(fit_curr_data[[2]])$adj.r.squared,
                      glance(fit_curr_data[[2]])$p.value)) 
    
  }
  fit_data <- as.matrix(fit_data)
  rownames(fit_data)<- group
  colnames(fit_data) <- c("corr","pval",'Fit','R-squered','P-value')
  return(fit_data)
}
elbow_method <- function(data,k){
  wss <- c()
  data <- data
  for (i in 2:k){
    wss<-c(wss,kmeans(data, i, nstart=50)$tot.withinss)
    print(i)
  
  }
  
  plot(2:k, wss,
       type="b", pch = 19, frame = FALSE, 
       xlab="Number of clusters K",
       ylab="Total within-clusters sum of squares")
  return(wss)
}
unit_clusters<- function(clus_num,clusters){
  names_lst <- c()
  for (idx in clus_num){
    names_lst <- c(names_lst,rownames(clusters)[which(as.numeric(clusters[,1]) == idx)])
  }
  return(names_lst)
}
write_clusters_raw_data <- function(clusters,raw_1,raw_2,label_1,label_2,range){
  res <- c()
  colnames(raw_1)<- paste(colnames(raw_1),label_1)
  colnames(raw_2)<- paste(colnames(raw_2),label_2)
  for ( i in range){
  
    names<- rownames(clusters[[i]])
    res <- rbind(res,cbind(raw_1[names,],raw_2[names,],clusters[[i]][,ncol(clusters[[i]])]))
    
  }
  write.csv(res,paste0(cluster_path,"raw_data ",label_1,' ',label_2,'.csv'))
}
get_clusters <- function(dat_mean_1,dat_mean_2,peak_dat_1,peak_dat_2,exp_idx,
                         label_1,label_2){
  
  # total_names <- intersect(p_low_89_up,p_low_222_up)
  slope_arr <- c()
  clusters_arr <- c()
  total_names <- c()
  cluster_regime <- c(0.65,0.8,1.19,1.4,1.6,1.8,2.19,3,4,5)
  curr_names_list <- rownames(fit_mat)[which(as.numeric(fit_mat[,exp_idx]) > 0 )]
  for(idx in 1:(length(cluster_regime))){
    # names_2 <- unit_clusters(clusters_222[[idx]],clusters)
    
    names_2 <- names(which(fit_mat[curr_names_list,exp_idx] <= cluster_regime[idx]))
    if(idx == length(cluster_regime)-1){
      names_2 <- curr_names_list
    }
    curr_names_list <- setdiff(curr_names_list,names_2)
    total_names[[idx]]<- names_2
    
    # names_2 <- intersect(names(which(tmp_km$cluster ==idx)),rownames(big_mat_fit_fitted_by_max)[which(big_mat_fit_fitted_by_max[,exp_idx] !=0)])
    if(length(names_2) >=1){
      df_2<- as.data.frame(cbind(peak_dat_1[names_2,1],peak_dat_2[names_2,1]))
      colnames(df_2) <- c(label_1,label_2)
      fit1 <- choose_fit(label_1,label_2,df_2) ##cluster genes
      slope_arr <-rbind(slope_arr,cbind(round(fit1$coefficients[1],2),
                                        round(summary(fit1)$adj.r.squared,2),
                                        format(glance(fit1)$p.value, scientific=TRUE),
                                        length(names_2)))
      options(digits = 2)
      ggplotRegression_narrow(df_2,fit1,paste0(cluster_path,"cluster_ ",i,' ' 
                                               ,idx,' ',round(fit1$coefficients[1],2))
                              ,max(df_2),0,TRUE,0)
      
      clusters_arr[[idx]] <-cbind(fit_mat[names_2,],rep(slope_arr[idx],length(names_2)))
      
      df_2_corrected<- convert_data_to_scatter_plot(as.matrix(dat_mean_1[names_2,2:6]),
                                                    as.matrix(fit1$coefficients[1]*dat_mean_2[names_2,2:6]))
      colnames(df_2_corrected) <- c(label_1,label_2)
      
      fit1_corrected <- choose_fit(label_1,label_2,df_2_corrected)
    }
    else if(length(names_2) == 1){
      
      slope_arr <-rbind(slope_arr,cbind(round(fit_mat[names_2,exp_idx],2),
                                        NA,
                                        NA,
                                        length(names_2)))
    }
    else{
      slope_arr <- rbind(slope_arr,c(NA,NA,NA,NA))
    }
    
    
  }
  
  colnames(slope_arr) <- c('slope','R2','p-value','size')
  rownames(slope_arr)<- 1:length(cluster_regime)
  write.csv(slope_arr,paste0(cluster_path,'slope_arr.csv'))
  
  return(clusters_arr)
}

cc
start <- Reduce(union,list(intersect(p_low_21_10_up,p_low_89_up),intersect(p_low_21_10_up,p_low_222_up),intersect(p_low_89_up,p_low_222_up)))
## creat/e initial fit matrix
big_mat_fit_fitted_by_max <- matrix(0,length(gene_pool),6) ## all fit - wirg p val - tommy
rownames(big_mat_fit_fitted_by_max) <- gene_pool
colnames(big_mat_fit_fitted_by_max) <- c('Fit pic 14 1 to pic 10',
                           'Fit kcl 10 to pic 10','Fit kcl 10 to pic 14 1',
                           'Peak time Pic 10',
                           'Peak time Pic 14 1',
                           'Peak time kcl 10')



big_mat_fit_fitted_by_max[intersect(p_low_21_10_up,p_low_89_up),c(1)] <-
  round(peak_place_and_value_rep_21_10[intersect(p_low_21_10_up,p_low_89_up),1]/
          peak_place_and_value_rep_89[intersect(p_low_21_10_up,p_low_89_up),1],2)
big_mat_fit_fitted_by_max[intersect(p_low_21_10_up,p_low_222_up),c(2)]<-
  round(peak_place_and_value_rep_21_10[intersect(p_low_21_10_up,p_low_222_up),1]/
          peak_place_and_value_rep_222[intersect(p_low_21_10_up,p_low_222_up),1],2)
big_mat_fit_fitted_by_max[intersect(p_low_222_up,p_low_89_up),c(3)]<-
  round(peak_place_and_value_rep_89[intersect(p_low_222_up,p_low_89_up),1]/
          peak_place_and_value_rep_222[intersect(p_low_222_up,p_low_89_up),1],2)

# remove zeros
# if there are more than two cor which passed than the thired fit will pass too
big_mat_fit_fitted_by_max <- big_mat_fit_fitted_by_max[which(rowSums(big_mat_fit_fitted_by_max)>0),]

names_to_add <- Reduce(intersect,list(which(big_mat_fit_fitted_by_max[,1] == 0),which(big_mat_fit_fitted_by_max[,2] != 0),

                                                                            which(big_mat_fit_fitted_by_max[,3] != 0)))
big_mat_fit_fitted_by_max[names_to_add,1] <-round(peak_place_and_value_rep_21_10[names_to_add,1]/peak_place_and_value_rep_89[names_to_add,1],2)

names_to_add <- Reduce(intersect,list(which(big_mat_fit_fitted_by_max[,2] == 0),which(big_mat_fit_fitted_by_max[,1] != 0),
                                      which(big_mat_fit_fitted_by_max[,3] != 0)))
big_mat_fit_fitted_by_max[names_to_add,2] <-round(peak_place_and_value_rep_21_10[names_to_add,1]/peak_place_and_value_rep_222[names_to_add,1],2)


names_to_add <- Reduce(intersect,list(which(big_mat_fit_fitted_by_max[,3] == 0),which(big_mat_fit_fitted_by_max[,1] != 0),
                                      which(big_mat_fit_fitted_by_max[,2] != 0)))

big_mat_fit_fitted_by_max[names_to_add,3] <-round(peak_place_and_value_rep_89[names_to_add,1]/peak_place_and_value_rep_222[names_to_add,1],2)


##add peak time
big_mat_fit_fitted_by_max[rownames(big_mat_fit_fitted_by_max),c(5)]<-peak_place_and_value_rep_89[rownames(big_mat_fit_fitted_by_max),2]
big_mat_fit_fitted_by_max[rownames(big_mat_fit_fitted_by_max),c(6)]<-peak_place_and_value_rep_222[rownames(big_mat_fit_fitted_by_max),2]
big_mat_fit_fitted_by_max[rownames(big_mat_fit_fitted_by_max),c(4)]<-peak_place_and_value_rep_21_10[rownames(big_mat_fit_fitted_by_max),2]

km_mat_ <- as.data.frame(cbind(round(big_mat_fit_fitted_by_max[,1],2)))
km_mat <- as.data.frame(km_mat_[which(km_mat_ >0),])
rownames(km_mat) <-rownames(km_mat_)[which(km_mat_ >0)]


tmp<-silhouette_score(df = km_mat,k = 2:100)
elbow_method(km_mat,100)

cluster_df_222 <- get_clusters(mean_fc_222,mean_fc_21_10,peak_place_and_value_rep_222,peak_place_and_value_rep_21_10,
                               2,exp[1],exp[2])

cluster_df_89 <- get_clusters(mean_fc_89,mean_fc_21_10,peak_place_and_value_rep_89,peak_place_and_value_rep_21_10,
                               1,exp[3],exp[2])

cluster_df_222_89 <- get_clusters(mean_fc_222,mean_fc_89,peak_place_and_value_rep_222,peak_place_and_value_rep_89,
                              3,exp[1],exp[3])

write_clusters_raw_data(cluster_df_222,data_21_10_esat_t_normalized,data_222_esat_t_normalized,exp[2],exp[1],1:9)
write_clusters_raw_data(cluster_df_89,data_21_10_esat_t_normalized,data_89_esat_t_normalized,exp[2],exp[3],1:9)
write_clusters_raw_data(cluster_df_222_89,data_89_esat_t_normalized,data_222_esat_t_normalized,exp[3],exp[1],1:9)



########## only pic 14 part #########################################################################

big_mat_fit_fitted_by_max_14 <- matrix(0,length(gene_pool_pic_14),3) ## all fit - wirg p val - tommy
rownames(big_mat_fit_fitted_by_max_14) <- gene_pool_pic_14
colnames(big_mat_fit_fitted_by_max_14) <- c('Fit pic 14 2 to pic 14 1',
                                         'Peak time Pic 14 1',
                                         'Peak time Pic 14 2')



big_mat_fit_fitted_by_max_14[intersect(p_low_21_14_only_14_up,p_low_89_only_14_up),c(1)] <-
  round(peak_place_and_value_rep_89[intersect(p_low_21_14_only_14_up,p_low_89_only_14_up),1]/
          peak_place_and_value_rep_21_14[intersect(p_low_21_14_only_14_up,p_low_89_only_14_up),1],2)

# remove zeros
# if there are more than two cor which passed than the thired fit will pass too
big_mat_fit_fitted_by_max_14 <- big_mat_fit_fitted_by_max_14[which(rowSums(big_mat_fit_fitted_by_max_14)>0),]


##add peak time
big_mat_fit_fitted_by_max_14[rownames(big_mat_fit_fitted_by_max_14),c(2)]<-peak_place_and_value_rep_89[rownames(big_mat_fit_fitted_by_max_14),2]
big_mat_fit_fitted_by_max_14[rownames(big_mat_fit_fitted_by_max_14),c(3)]<-peak_place_and_value_rep_21_14[rownames(big_mat_fit_fitted_by_max_14),2]

fit_mat <- big_mat_fit_fitted_by_max_14

cluster_df_pic_14_1_2 <- get_clusters(mean_fc_21_14,mean_fc_89,peak_place_and_value_rep_21_14,peak_place_and_value_rep_89,
                               1,exp[4],exp[3])

cluster_df_pic_14_1_2[[9]] <- c("Mir129-2",t(as.matrix(cluster_df_pic_14_1_2[[9]]))[1,])
write_clusters_raw_data(cluster_df_pic_14_1_2,data_89_esat_t_normalized,data_21_14_esat_t_normalized,exp[3],exp[4],1:8)
