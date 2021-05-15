# plot_gene_ont <-function(gene_ont_data){
#   
#   below_0.05 <- which(gene_ont_data$p.adjust < 0.05)
#   
#     data_names <- gene_ont_data$Description[below_0.05]
#   data <- data.frame(data_names,-log10(gene_ont_data$p.adjust[below_0.05]),gene_ont_data$Count[below_0.05])
#   rownames(data)<- data_names
#   colnames(data) <- c("group",'p.adjusted','Count')
#   # png(paste0(title,".png"), width = 8, height = 8, units = 'in', res = 300)
#   print({ggplot(as.data.frame(data), aes(x  = p.adjusted, y = group)) + geom_point(aes(size = Count)) +
#       scale_size_continuous(range = c(10, 15)) +geom_count(aes(color = Count,size =Count))+
#       geom_text(aes(label = Count))+ scale_x_discrete(limits=c(-log10(0.05), -log10(0.04),-log10(0.03),-log10(0.02),-log10(0.01),0) ,
#                                                      labels = round(c(-log10(0.05), -log10(0.04),-log10(0.03),-log10(0.02),-log10(0.01),0),2))+
#               guides(color = 'legend')+ggtitle(title)})
#         dev.off()
#   
#         
#   plot(tmp[,1],tmp[,2])
#   
# }

plot_gene_ont <-function(gene_ont_data,title){
  
  
  
  data_names <- gene_ont_data$Gene.Set.Name
  data <- data.frame(data_names,-log10(gene_ont_data$p.value),as.numeric(gene_ont_data$X..Genes.in.Overlap..k.))
  rownames(data)<- data_names
  colnames(data) <- c("group",'p.value','Count')
  data$group <- factor(data$group, levels = data$group[order(data$p.value)])
  # png(paste0(title,".png"), width = 8, height = 8, units = 'in', res = 300)
  print({ggplot(data, aes(x  = p.value, y = group)) + geom_point(stat = "identity",aes(size = Count)) +
      scale_size_continuous(range = c(4, 7)) +geom_count(aes(color = Count,size =Count))+
      geom_text(aes(label = Count))+ scale_x_discrete(limits=c(-log10(0.05), -log10(0.04),-log10(0.03),-log10(0.02),-log10(0.01),0) ,
                                                      labels = round(c(-log10(0.05), -log10(0.04),-log10(0.03),-log10(0.02),-log10(0.01),0),2))+
      guides(color = 'legend')+ggtitle(title)})
  dev.off()
  
  
  plot(tmp[,1],tmp[,2])
  
}
find_intersection<- function(struct_1,struct_2, label_1,label_2){
  res <- c()
  Flag = FALSE
  group_size<- c()
  res_mat <- list()
  for(j in 1:length(struct_1)){
    for(i in 1:length(struct_2)){
      
      res[[i]]<- intersect(rownames(struct_1[[j]]),rownames(struct_2[[i]]))
      group_size <- rbind(group_size,c(length(res[[i]]),length(rownames(struct_1[[j]])),length(rownames(struct_2[[i]])),
                          struct_1[[j]][1,ncol(struct_1[[j]])],struct_2[[i]][1,ncol(struct_2[[i]])]))
                          
      if(i >1){
        max_len<- max(max_len,length(res[[i]]))
      }
      else{
        max_len = length(res[[i]])
      }
      
    }
    res_mat_i <- matrix("empty",length(struct_1),max_len)
    
    
    for(k in 1:length(struct_2)){
      if(length(res[[k]])> 0){
        res_mat_i[k,1:length(res[[k]])] = res[[k]]
        
      }
    }
      res_mat_i <- cbind(group_size,res_mat_i)
      colnames(res_mat_i)<- c("intersect",paste(label_1,"size"),
                              paste(label_2,"size"),paste(label_1,"slope"),paste(label_2,"slope"),
                              paste("gene",1:max_len))
      rownames(res_mat_i)<- paste("cluster",1:length(struct_2))
      res_mat[[j]] <- as.matrix(res_mat_i)
      write.xlsx(res_mat_i, file=paste(label_1,"intersect",label_2,".xls"), append=Flag,sheetName= paste("cluster",j))
      Flag = TRUE
      group_size<- c()
      
    
    
  }
  
  
  return(res_mat)
}

cc
get_gene_ont_struct <- function(path,file_name){
  library(xlsx)
  Flag = FALSE
  Files=list.files(path=path, pattern=".xls")
  
  for (i in 1:length(Files)){
    gene_ont_matches <- read.csv(paste0(path,Files[i]),sep = '\t',nrows = 4,skip = 3)
    genes_ass <- read.csv(paste0(path,Files[i]),sep = '\t',nrows = gene_ont_matches[1,2],skip = 9)
    write.xlsx(genes_ass, file=paste0(path,file_name,"_gene_ont_terms.xls"), append=Flag,sheetName= Files[i])
  
    genes_ont_terms<- read.csv(paste0(path,Files[i]),sep = '\t',skip = 14+gene_ont_matches[1,2])
    write.xlsx(genes_ont_terms, file=paste0(path,file_name,"_gene_description.xls"),append=Flag ,sheetName=Files[i])
    Flag = TRUE    
  }
}
cc

library(org.Mm.eg.db)
path_89 <- 'C:/Users/Lior/Desktop/paper/new_org_code/stage 3 - cluster + gene ont/clusters  - pairwise/89/only 89/'
get_gene_ont_struct("overlap_89",path_89,1:9)
path_222 <- 'C:/Users/Lior/Desktop/paper/new_org_code/stage 3 - cluster + gene ont/clusters  - pairwise/222/only 222/'
get_gene_ont_struct("overlap_222",path_222,c(2:5,7,8))

path_222_89 <- 'C:/Users/Lior/Desktop/paper/new_org_code/stage 3 - cluster + gene ont/clusters  - pairwise/222_89/only 222_89/'
get_gene_ont_struct(path_222_89,'cluster_222_89')

path_89_21_14 <- 'C:/Users/Lior/Desktop/paper/new_org_code/stage 3 - cluster + gene ont/clusters  - pairwise/89_21_14/'
get_gene_ont_struct(path_89_21_14,'cluster_222_89')


plot_gene_ont(gene_TF_file," ")

for(i in 1:length(cluster_df_pic_14_1_2)){
   print(nrow(cluster_df_pic_14_1_2[[i]]))
  #print(cluster_df_pic_14_1_2[[i]][1,4])
}

int_89_222<- find_intersection(cluster_df_222,cluster_df_89,"kcl 10 fit pic 10","pic 14 1 fit pic 10")
int_89_222_222<- find_intersection(cluster_df_222,cluster_df_222_89,"kcl 10 fit pic 10","kcl 10 fit pic 14 1")
int_89_222_89<- find_intersection(cluster_df_89,cluster_df_222_89,"pic 14 1 fit pic 10","kcl 10 fit pic 14 1")
