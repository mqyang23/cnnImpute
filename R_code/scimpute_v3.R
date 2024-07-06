#load functions
find_hv_genes<-function (count, I, J) 
{
  count_nzero = lapply(1:I, function(i) setdiff(count[i, ], 
                                                log10(1.01)))
  mu = sapply(count_nzero, mean)
  mu[is.na(mu)] = 0
  sd = sapply(count_nzero, sd)
  sd[is.na(sd)] = 0
  cv = sd/mu
  cv[is.na(cv)] = 0
  high_var_genes = which(mu >= 1 & cv >= quantile(cv, 0.25))
  if (length(high_var_genes) < 500) {
    high_var_genes = 1:I
  }
  count_hv = count[high_var_genes, ]
  return(count_hv)
}

find_neighbors<-function (count_hv, labeled, J, Kcluster = NULL, ncores, cell_labels = NULL) 
{
  if (labeled == TRUE) {
    if (class(cell_labels) == "character") {
      labels_uniq = unique(cell_labels)
      labels_mth = 1:length(labels_uniq)
      names(labels_mth) = labels_uniq
      clust = labels_mth[cell_labels]
    }
    else {
      clust = cell_labels
    }
    nclust = length(unique(clust))
    print("calculating cell distances ...")
    dist_list = lapply(1:nclust, function(ll) {
      cell_inds = which(clust == ll)
      count_hv_sub = count_hv[, cell_inds, drop = FALSE]
      if (length(cell_inds) < 1000) {
        var_thre = 0.4
        pca = prcomp(t(count_hv_sub))
        eigs = (pca$sdev)^2
        var_cum = cumsum(eigs)/sum(eigs)
        if (max(var_cum) <= var_thre) {
          npc = length(var_cum)
        }
        else {
          npc = which.max(var_cum > var_thre)
          if (labeled == FALSE) {
            npc = max(npc, Kcluster)
          }
        }
      }
      else {
        var_thre = 0.6
        pca = rpca(t(count_hv_sub), k = 1000, center = TRUE, 
                   scale = FALSE)
        eigs = (pca$sdev)^2
        var_cum = cumsum(eigs)/sum(eigs)
        if (max(var_cum) <= var_thre) {
          npc = length(var_cum)
        }
        else {
          npc = which.max(var_cum > var_thre)
          if (labeled == FALSE) {
            npc = max(npc, Kcluster)
          }
        }
      }
      if (npc < 3) {
        npc = 3
      }
      mat_pcs = t(pca$x[, 1:npc])
      dist_cells_list = mclapply(1:length(cell_inds), 
                                 function(id1) {
                                   d = sapply(1:id1, function(id2) {
                                     sse = sum((mat_pcs[, id1] - mat_pcs[, id2])^2)
                                     sqrt(sse)
                                   })
                                   return(c(d, rep(0, length(cell_inds) - id1)))
                                 }, mc.cores = ncores)
      dist_cells = matrix(0, nrow = length(cell_inds), 
                          ncol = length(cell_inds))
      for (cellid in 1:length(cell_inds)) {
        dist_cells[cellid, ] = dist_cells_list[[cellid]]
      }
      dist_cells = dist_cells + t(dist_cells)
      return(dist_cells)
    })
    return(list(dist_list = dist_list, clust = clust))
  }
  if (labeled == FALSE) {
    print("dimension reduction ...")
    if (J < 5000) {
      var_thre = 0.4
      pca = prcomp(t(count_hv))
      eigs = (pca$sdev)^2
      var_cum = cumsum(eigs)/sum(eigs)
      if (max(var_cum) <= var_thre) {
        npc = length(var_cum)
      }
      else {
        npc = which.max(var_cum > var_thre)
        if (labeled == FALSE) {
          npc = max(npc, Kcluster)
        }
      }
    }
    else {
      var_thre = 0.6
      pca = rpca(t(count_hv), k = 1000, center = TRUE, 
                 scale = FALSE)
      eigs = (pca$sdev)^2
      var_cum = cumsum(eigs)/sum(eigs)
      if (max(var_cum) <= var_thre) {
        npc = length(var_cum)
      }
      else {
        npc = which.max(var_cum > var_thre)
        if (labeled == FALSE) {
          npc = max(npc, Kcluster)
        }
      }
    }
    if (npc < 3) {
      npc = 3
    }
    mat_pcs = t(pca$x[, 1:npc])
    print("calculating cell distances ...")
    dist_cells_list = mclapply(1:J, function(id1) {
      d = sapply(1:id1, function(id2) {
        sse = sum((mat_pcs[, id1] - mat_pcs[, id2])^2)
        sqrt(sse)
      })
      return(c(d, rep(0, J - id1)))
    }, mc.cores = ncores)
    dist_cells = matrix(0, nrow = J, ncol = J)
    for (cellid in 1:J) {
      dist_cells[cellid, ] = dist_cells_list[[cellid]]
    }
    dist_cells = dist_cells + t(dist_cells)
    min_dist = sapply(1:J, function(i) {
      min(dist_cells[i, -i])
    })
    iqr = quantile(min_dist, 0.75) - quantile(min_dist, 
                                              0.25)
    outliers = which(min_dist > 1.5 * iqr + quantile(min_dist, 
                                                     0.75))
    non_out = setdiff(1:J, outliers)
    spec_res = specc(t(mat_pcs[, non_out]), centers = Kcluster, 
                     kernel = "rbfdot")
    print("cluster sizes:")
    print(spec_res@size)
    nbs = rep(NA, J)
    nbs[non_out] = spec_res
    return(list(dist_cells = dist_cells, clust = nbs))
  }
}

find_va_genes<-function (parslist, subcount) 
{
  point = log10(1.01)
  valid_genes = which((rowSums(subcount) > point * ncol(subcount)) & 
                        complete.cases(parslist))
  if (length(valid_genes) == 0) 
    return(valid_genes)
  mu = parslist[, "mu"]
  sgene1 = which(mu <= log10(1 + 1.01))
  dcheck1 = dgamma(mu + 1, shape = parslist[, "alpha"], rate = parslist[, 
                                                                        "beta"])
  dcheck2 = dnorm(mu + 1, mean = parslist[, "mu"], sd = parslist[, 
                                                                 "sigma"])
  sgene3 = which(dcheck1 >= dcheck2 & mu <= 1)
  sgene = union(sgene1, sgene3)
  valid_genes = setdiff(valid_genes, sgene)
  return(valid_genes)
}

get_mix<-function (xdata, point) 
{
  inits = rep(0, 5)
  inits[1] = sum(xdata == point)/length(xdata)
  if (inits[1] == 0) {
    inits[1] = 0.01
  }
  inits[2:3] = c(0.5, 1)
  xdata_rm = xdata[xdata > point]
  inits[4:5] = c(mean(xdata_rm), sd(xdata_rm))
  if (is.na(inits[5])) {
    inits[5] = 0
  }
  paramt = inits
  eps = 10
  iter = 0
  loglik_old = 0
  while (eps > 0.5) {
    wt = calculate_weight(xdata, paramt)
    paramt[1] = sum(wt[, 1])/nrow(wt)
    paramt[4] = sum(wt[, 2] * xdata)/sum(wt[, 2])
    paramt[5] = sqrt(sum(wt[, 2] * (xdata - paramt[4])^2)/sum(wt[, 
                                                                 2]))
    paramt[2:3] = update_gmm_pars(x = xdata, wt = wt[, 1])
    loglik = sum(log10(dmix(xdata, paramt)))
    eps = (loglik - loglik_old)^2
    loglik_old = loglik
    iter = iter + 1
    if (iter > 100) 
      break
  }
  return(paramt)
}


dmix<-function (x, pars) 
{
  pars[1] * dgamma(x, shape = pars[2], rate = pars[3]) + (1 - 
                                                            pars[1]) * dnorm(x, mean = pars[4], sd = pars[5])
}

### root-finding equation
fn = function(alpha, target){
  log(alpha) - digamma(alpha) - target
}


get_mix_parameters<-function (count, point = log10(1.01), path, ncores = 8) 
{
  count = as.matrix(count)
  null_genes = which(abs(rowSums(count) - point * ncol(count)) < 
                       1e-10)
  parslist = mclapply(1:nrow(count), function(ii) {
    if (ii%%2000 == 0) {
      gc()
      print(ii)
    }
    if (ii %in% null_genes) {
      return(rep(NA, 5))
    }
    xdata = count[ii, ]
    paramt = try(get_mix(xdata, point), silent = TRUE)
    if (class(paramt) == "try-error") {
      paramt = rep(NA, 5)
    }
    return(paramt)
  }, mc.cores = ncores)
  save(parslist, file = path)
  parslist = Reduce(rbind, parslist)
  colnames(parslist) = c("rate", "alpha", "beta", "mu", "sigma")
  saveRDS(parslist, file = path)
}

calculate_weight<-function (x, paramt) 
{
  pz1 = paramt[1] * dgamma(x, shape = paramt[2], rate = paramt[3])
  pz2 = (1 - paramt[1]) * dnorm(x, mean = paramt[4], sd = paramt[5])
  pz = pz1/(pz1 + pz2)
  pz[pz1 == 0] = 0
  return(cbind(pz, 1 - pz))
}

update_gmm_pars = function(x, wt){
  tp_s = sum(wt)
  tp_t = sum(wt * x)
  tp_u = sum(wt * log(x))
  tp_v = -tp_u / tp_s - log(tp_s / tp_t)
  if (tp_v <= 0){
    alpha = 20
  }else{
    alpha0 = (3 - tp_v + sqrt((tp_v - 3)^2 + 24 * tp_v)) / 12 / tp_v
    if (alpha0 >= 20){alpha = 20
    }else{
      alpha = uniroot(fn, c(0.9, 1.1) * alpha0, target = tp_v, 
                      extendInt = "yes")$root
    }
  }
  ## need to solve log(x) - digamma(x) = tp_v
  ## We use this approximation to compute the initial value
  beta = tp_s / tp_t * alpha
  return(c(alpha, beta))
}


library(parallel)
library(doParallel)
library(kernlab)
#main func
out_dir="cnnImpute/tmp_file/"
#ncores=10
Kcluster=2
labeled=FALSE
drop_thre = 0.5
infile = "csv"
type = "count"
genelen = NULL
filetype = infile

#calculate dropout probabilty of the data

count<-read.csv("cnnImpute/tmp_file/data_original.csv",header = TRUE,row.names=1)
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

count_hv = find_hv_genes(count, I, J)
print("searching candidate neighbors ... ")
neighbors_res = find_neighbors(count_hv = count_hv, labeled =TRUE, J = J,
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
  get_mix_parameters(count = count[, which(clust == cc), drop = FALSE],
                                point = log10(1.01),
                                path = paste0(out_dir, "pars", cc, ".rds"), ncores = ncores)
  cells = which(clust == cc)
  if(length(cells) <= 1){ next }
  parslist = readRDS(paste0(out_dir, "pars", cc, ".rds"))
  print("searching for valid genes ...")
  valid_genes = find_va_genes(parslist, subcount = count[, cells])
  if(length(valid_genes) <= 10){ next }

  subcount = count[valid_genes, cells, drop = FALSE]
  subcount1 = count1[valid_genes, cells, drop = FALSE]
  rowname<-rowname[valid_genes]
  Ic = length(valid_genes)
  Jc = ncol(subcount)
  parslist = parslist[valid_genes, , drop = FALSE]

  droprate = t(sapply(1:Ic, function(i) {
    wt = calculate_weight(subcount[i, ], parslist[i, ])
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
  write.csv(droprate1,file=paste0("cnnImpute/tmp_file/droprate",cc,"_withlabel.csv"))
  write.table(rowname,file=paste0("cnnImpute/tmp_file/rowname",cc,"_withlabel.txt"))
}

