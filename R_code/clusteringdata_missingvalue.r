#load packages
library(stats)
library(Rtsne)
library(ADPclust)

#list all files in this directory
list<-list.files("cnnimpute/input/")


for(file in list){
  
  mat1<-read.csv(paste0("cnnimpute/input/", file),header = TRUE,row.names=1)
  
  # #remove cell that not expressed in all genes
  # mat1<- mat1[,-which(colSums(mat1)==0)]
  # #remove genes that not expressed in all cells
  # mat1<-mat1[-which(rowSums(mat1)==0),]
  
  
  #set seed for each step
  set.seed(1)
  #normalize data
  tsne_input <- log2(t(t(mat1)/colSums(mat1))*1000000+1)
  set.seed(1)
  
  #use Rrtsne to do dimension reduction
  if(ncol(mat1)<=30){
    tsne_output <- Rtsne(t(tsne_input), dims = 3, perplexity = 5, check_duplicates = FALSE, num_threads = 0)
	}else{
  	tsne_output <- Rtsne(t(tsne_input), dims = 3, perplexity = 30, check_duplicates = FALSE, num_threads = 0)
  }
  print("tsne finished")
  
  #cluster data by using adpclust and kmeans
  set.seed(1)
  adpOUTPUT <- adpclust(tsne_output$Y, htype = "amise",centroids="auto")
  set.seed(1)
  tsne_kmeansOUTPUT <- kmeans(tsne_output$Y, tsne_output$Y[adpOUTPUT$centers,], adpOUTPUT$nclust)
  
  #print("cluster result saved")
  #write.csv(file=paste0("/home/junie/data/",i),tsne_kmeansOUTPUT)  
  
  #label data with our new label:
  colname=colnames(mat1)
  for (i in 1:max(tsne_kmeansOUTPUT$cluster)) {
    loc1<-which(tsne_kmeansOUTPUT$cluster==i)
    colname[loc1]<-paste0(i,"_",seq(1,length(loc1)))
  }
  colnames(mat1)<-colname
  print("data labeled")
  #save labeled data
  write.csv(file=paste0("cnnimpute/tmp_file/",file),mat1) 
  
}

