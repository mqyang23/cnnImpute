
library(scImpute)
library(parallel)
library(doParallel)
library(kernlab)
#main func
out_dir="cnnimpute/tmp_file/"
#ncores=10
Kcluster=2
labeled=FALSE
drop_thre = 0.5
infile = "csv"
type = "count"
genelen = NULL
filetype = infile

#scimpute way to calculate dropout probabilty of the data

count<-read.csv("cnnimpute/tmp_file/data_original.csv",header = TRUE,row.names=1)
count = as.matrix(count)
count1 = count
print(paste("number of genes in raw count matrix", nrow(count)))
print(paste("number of cells in raw count matrix", ncol(count)))
totalCounts_by_cell = colSums(count)
saveRDS(totalCounts_by_cell, file = paste0(out_dir, "totalCounts_by_cell.rds"))
totalCounts_by_cell[totalCounts_by_cell == 0] = 1
count = sweep(count, MARGIN = 2, 10^6/totalCounts_by_cell, FUN = "*")
if (min(count) < 0) {
        stop("smallest read count cannot be negative!")
    }
count= log10(count + 1.01)
genenames = rownames(count)
cellnames = colnames(count)
#label=substring(colnames(count),2,2)
label=sub('.', '', gsub("\\_.*","",colnames(count)))
label=as.numeric(label)
cell_labels=label
drop_thre=0.5
ncores=as.numeric(detectCores())

count = as.matrix(count)
#rownames(count)<-seq(0,nrow(count)-1)
I = nrow(count)
J = ncol(count)
count_imp = count

count_hv = scImpute:::find_hv_genes(count, I, J)
print("searching candidate neighbors ... ")
neighbors_res = scImpute:::find_neighbors(count_hv = count_hv, labeled =TRUE, J = J,
                                          ncores = ncores, cell_labels = cell_labels)
dist_list = neighbors_res$dist_list
clust = neighbors_res$clust

# mixture model
nclust = sum(!is.na(unique(clust)))
cl = makeCluster(ncores, outfile="")
registerDoParallel(cl)

for(cc in 1:nclust){
  rowname<-seq(0,nrow(count)-1)
  print(paste("estimating dropout probability for type", cc, "..."))
  paste0(out_dir, "pars", cc, ".rds")
  scImpute:::get_mix_parameters(count = count[, which(clust == cc), drop = FALSE],
                                point = log10(1.01),
                                path = paste0(out_dir, "pars", cc, ".rds"), ncores = ncores)
  cells = which(clust == cc)
  if(length(cells) <= 1){ next }
  parslist = readRDS(paste0(out_dir, "pars", cc, ".rds"))
  print("searching for valid genes ...")
  valid_genes = scImpute:::find_va_genes(parslist, subcount = count[, cells])
  if(length(valid_genes) <= 10){ next }

  subcount = count[valid_genes, cells, drop = FALSE]
  subcount1 = count1[valid_genes, cells, drop = FALSE]
  rowname<-rowname[valid_genes]
  Ic = length(valid_genes)
  Jc = ncol(subcount)
  parslist = parslist[valid_genes, , drop = FALSE]

  droprate = t(sapply(1:Ic, function(i) {
    wt = scImpute:::calculate_weight(subcount[i, ], parslist[i, ])
    return(wt[, 1])
  }))
  #rownames(droprate)<-rownames(subcount)
  mucheck = sweep(subcount, MARGIN = 1, parslist[, "mu"], FUN = ">")
  droprate[mucheck & droprate > drop_thre] = 0
  
  #what I added
  #new
  loc<-intersect(which(droprate > 0.5), which(subcount1 != 0))
  if(length(loc) > 0) {
    droprate[loc] <- 0.4
  }
  
  # for (i in 1:nrow(droprate)) {
  #   for (j in 1:ncol(droprate)) {
  #     if((droprate[i,j] > 0.5) & (subcount1[i,j] != 0)){
  #       droprate[i,j] = 0.4
  #     }
  #   }
  # }
  
  #what I add, to map valid_genes matrix to origin matirx
  droprate1<-matrix(c(rep(1,nrow(count)*length(cells))),nrow = nrow(count))
  rownames(droprate1)<-seq(0,nrow(count)-1)
  #new
  loc1<-which(rownames(droprate1) %in% rowname)
  if(length(loc1) > 0) {
    droprate1[loc1,]<-droprate[which(rowname%in%rownames(droprate1)[loc1]),]
  }
  
  ###old one, not used
  # for (i in 1:nrow(droprate1)) {
  #   if (rownames(droprate1)[i] %in% rowname){
  #     loc<-which(rownames(droprate1)[i]==rowname)
  #     droprate1[i,]<-droprate[loc,]
  #   }
  # }
  
  rm(subcount)
  write.csv(droprate1,file=paste0("cnnimpute/tmp_file/droprate",cc,"_withlabel.csv"))
  write.table(rowname,file=paste0("cnnimpute/tmp_file/rowname",cc,"_withlabel.txt"))
}

