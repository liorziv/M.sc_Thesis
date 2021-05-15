library("stats")
#############################################################################################################
# get a path which contains gene.txt files and create a normalized matrix (deseq2)
# return a pca plot with the normalized data (deseq2)
# files_path - the path to the gene.txt files
# file_names_ordered - all the file names (A1,A2 ...) according to the condition_names
# condition_names - conditions for the deseq2 normaliztion
# match - name of pca output file
#############################################################################################################
normalize_data_sep_esat <- function(files_path, file_names_ordered, condition_names, pca_file_name,data,title){
  
  library(DESeq2)
  #for the deseq2
  condition = data.frame(condition =file_names_ordered)
  colData = data.frame(condition = condition_names)
  rownames(colData) <- file_names_ordered
  
  #preper for normaliztion
  dds <- DESeqDataSetFromMatrix(as.matrix(data), colData = colData, design=   ~condition )
  dds$condition <- relevel(dds$condition, ref = "cont")
  vsd <- vst(dds, blind=FALSE)
  
  plt <- plotPCA(vsd, intgroup='condition')
  plot(plt)
  print(plt)
  dev.copy(png,paste0(files_path,pca_file_name))
  dev.off()
  
  #show libery depth - if smaller than 0.3 should be suspicious
  colnames(data) <- condition_names
  dds <- estimateSizeFactors(dds)
  print(sizeFactors(dds))
  
  
  dds<-DESeq(dds,minReplicatesForReplace = Inf)# betaPrior = FALSE)

  #get the normalized counts
  data_norm<- counts(dds, normalized=TRUE)
  colnames(data_norm) <- condition_names
 
  
  # plot_dendogram(data_norm,files_path,"dendogram_esat.png")
  return (data_norm)
  
}

#############################################################################################################
#loads the data from the aligned files
#files_path - the path of the wanted files
#files_names_ordered - the file names in a specific order
#############################################################################################################
load_data_esat <-function(files_path, file_names_ordered){
  #preper to read the files
  list_of_files = list.files(path = files_path , pattern = "\\.gene.txt$")
  file_name = c()
  data = vector()
  tmp_data = c()
 
  
  #read files
  for (csv_file in list_of_files){
    print(csv_file)
    tmp_name = c(strsplit(csv_file,"_")[[1]][2])
    
    
    file_name <- c(file_name, tmp_name)
    #read each file and append the expression coulmn into a big matrix
    tmp_data <- read.csv(file = paste0(files_path, csv_file),  header=T, sep="\t", stringsAsFactors=FALSE)
    data =  cbind(data,as.numeric(tmp_data$Exp1))
    
  }
  
  #adding the gene names and name of file
  colnames(data) <- c(file_name)
  rownames(data) <-c(tmp_data$Symbol)
  i_out <- setdiff(colnames(data),file_names_ordered)
  print(i_out)
  data <-data[,setdiff(colnames(data),i_out)]
  i1 <- match(file_names_ordered , colnames(data),nomatch = NA_integer_, incomparables = NULL)
  data <-data[,c(i1)] 
  
  
  #remove zero reads genes 
  keep <- rowSums(data) > 0
  filtered_data <- data[keep,]
  
  
  return (filtered_data)
}

#############################################################################################################
# filters the data according to pre defined read amount 
#############################################################################################################
filter_data_by_read_num <-function(data,threshold){
  
  summedGenes = apply(data,1,sum)
  print(quantile(summedGenes, probs = seq(from = 0.1, to = 1, by = 0.1)))
  #filtered_data = data[which(summedGenes > 10),]
  indexs_to_keep <- which(data >= threshold, arr.ind=TRUE)
  num_of_appear <-count_appearence(indexs_to_keep[,1])
  indexs_to_keep <- as.matrix(unique(indexs_to_keep[,1]))
  indexs_to_keep_f <-which(num_of_appear> floor((length(data[1,]))/2))# takes only if at least half of the data was above 10
  filtered_data = data[c(indexs_to_keep_f),]
  print(quantile(apply(filtered_data,1,sum), probs = seq(from = 0.1, to = 1, by = 0.1)))
  return (as.matrix(filtered_data))
}

#############################################################################################################
# count the number of times a gene appear in the filter
#for example a number of repeats the reads of the gene where above 10 
#############################################################################################################
count_appearence <- function(data){
  data_range = unique(data)
  count_vec  <- matrix(0, nrow = data_range, ncol = 1)
  for (i in data_range){
    count_vec[i] <-sum(data == i)
  }
  return(count_vec)
}

#############################################################################################################
#filters the data using 3 condition:
# significant p value < 0.05
# the peak is at the versed condition (120 for example)
#the peaks fold change is 1.5 
# same for minima (0.66)
#############################################################################################################
filter_data_by_pvalue_and_peak<-function(dds,condition_names,title){
  res <- results(dds) #get results
  up_list <- c()
  up_list_fc<-c()
  down_list <-c()
  down_list_fc<-c()
  cond_by_order <- c("cont","30_min","60_min","120_min","240_min","480_min")
  
  data_norm <- counts(dds, normalized=TRUE) #get normalized data
  colnames(data_norm) <- condition_names
  mean_norm_data <-   t(aggregate(t(data_norm), by = list(colnames(data_norm)),mean)) #save the data mean
  new_col_names <-(mean_norm_data[1,])
  mean_norm_data <- mean_norm_data[2:nrow(mean_norm_data),]
  colnames(mean_norm_data) <- new_col_names
  
  #remove control vs control
  comp_names<- list(resultsNames(dds))
  comp_names <- as.matrix(comp_names[[1]])[-c(which(as.matrix(comp_names[[1]]) == "Intercept"))]
  #save the peak/minima place for each gene
  peak_data <- colnames(mean_norm_data)[apply(mean_norm_data[,1:6],1,which.max)]
  minima_data <- colnames(mean_norm_data)[apply(mean_norm_data[,1:6],1,which.min)]
  #for each condition we take the up and down regulated genes
  for (comp in comp_names){
    #take the current condiotion from comp names
    cond<-paste0(strsplit(comp,"_")[[1]][2],"_",strsplit(comp,"_")[[1]][3])
    cont_idx <- which(cond_by_order == cond)
    #take the specific results for that condition
    cond_res <- results(dds, name=comp, independentFiltering = FALSE,  cooksCutoff = Inf )
    cond_fc <- as.matrix(cond_res$log2FoldChange)
    rownames(cond_fc) <- rownames(data_norm)
    
    gene_names <- rownames(data_norm)
    
    keep_cond_peak <- which(peak_data == cond) # keep the peak of that condition
    keep_cond_minima <- which(minima_data == cond) #same for min value
    keep_cond_pval <- (which(cond_res$padj <= 0.05))
    keep_pval_no_Na <- keep_cond_pval[!is.na(cond_res$padj[keep_cond_pval])]
    keep_cond_1.5 <- which(as.integer(mean_norm_data[,"cont"])*1.5 < as.integer(mean_norm_data[,c(cond)]))
    
    keep_cond_0.6 <- which(as.integer(mean_norm_data[,"cont"])*(2/3) > as.integer(mean_norm_data[,c(cond)]))
    up_idxs <- intersect(keep_cond_peak,intersect(keep_cond_1.5,keep_pval_no_Na))
    down_idxs <- intersect(keep_cond_minima,intersect(keep_cond_0.6,keep_pval_no_Na))
    up_list <- c(up_list,gene_names[up_idxs])
    down_list <-  c(down_list,gene_names[down_idxs])
    up_list_fc <- rbind(up_list_fc,cbind(cond_fc[up_idxs,],rep(cont_idx,length(up_idxs))))
    down_list_fc <-  rbind(down_list_fc,cbind(cond_fc[down_idxs,],rep(cont_idx,length(down_idxs))))
  }
  down_list= setdiff(down_list,up_list)
  down_list_fc <- down_list_fc[rownames(down_list_fc) %in% down_list,]
  total_fc <- as.data.frame(rbind(up_list_fc,down_list_fc))
  double_tagges_genes <- intersect(rownames(up_list_fc),rownames(down_list_fc))
  remove_idxs <- match(double_tagges_genes ,rownames(down_list_fc))
  if(!isEmpty(remove_idxs)){
    total_fc <- as.data.frame(rbind(up_list_fc,down_list_fc[-c(remove_idxs),]))
  }
  #scatter plot with up and down regulated density
  total_filtered_data <- mean_norm_data[rownames(total_fc),]
  library(ggplot2)
  png(paste0(title,".png"), width = 8, height = 8, units = 'in', res = 300)
  print({ggplot(total_fc, aes(x = V2, y = V1)) +geom_jitter(size = 1)
    scale_x_discrete(limits=c( "2", "3","4","5","6") ,labels = c("30_min","60_min","120_min","240_min","480_min"))})
  dev.off()
  
  return(list(up_list,down_list))
}

take_only_coding<-function(data){
  
  #remove non coding genes
  coding_protein_list = read.table(file ="/media/liorz/myVolume/DORON/summary of filters/Dorons_project/sources/list_protein_coding_GRCm38.txt",sep="\n")
  filtered_data <- data[rownames(data) %in% coding_protein_list[,1],]
  
  return(filtered_data)
  
}
