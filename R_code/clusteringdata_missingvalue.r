FindH <- function(p, n, htype){
  if(!is.numeric(p)) stop('arg p must be numeric')
  if(length(p) != 1) stop('arg p must be scalar')    
  if(!is.numeric(n)) stop('arg n must be numeric')
  if(length(n) != 1) stop('arg n must be scalar')        
  if(!inherits(htype, 'character')) stop('arg htype not character')
  
  if(tolower(htype) == "amise"){
    return(AMISE(p, n))    
  }else if(tolower(htype) == "rot"){
    return(ROT(p, n))
  }else{
    stop("htype can only take two options 'amise' or 'rot'")
  }
}
AMISE <- function(x, y = NULL){
  if(is.null(y)){
    if(inherits(x, c("data.frame", "matrix"))){
      n <- nrow(x)
      p <- ncol(x)
    }
  }else{
    if(is.numeric(x) && x > 0 && is.numeric(y) && y > 0){
      p = x
      n = y
    }else
      stop("Wrong x, y.")
  }
  h <- (4 / (p + 2)) ^ (1 / (p + 4)) * n ^ (-1 / (p + 4))
  return(h)
}
FindFD <- function(distm, h, fdelta, f = NULL){
  # -------------------------------------------------------------------------    
  # Check arguments
  # -------------------------------------------------------------------------       
  if(!inherits(distm, 'dist')) stop("arg distm must inherit dist class. Got ", class(distm))
  
  n <- attr(distm, 'Size')
  distm <- as.matrix(distm)
  # ------------------------------------------------------------------------- 
  # Find f(x)
  # ------------------------------------------------------------------------- 
  if(is.null(f)){
    if(fdelta == "unorm"){
      f <- 1/(h * sqrt(2 * pi)) * rowSums(exp(-(distm/h)^2/2))
    }else if(fdelta == "weighted"){
      f <- rowSums(exp(-(distm/h)^2))        
    }else if(fdelta == "count"){
      f <- rowSums(distm < h) - 1   
    }else if(fdelta == "mnorm"){
      f <- rowSums(exp(-(distm / h) ^ 2 / 2))
    }else{
      stop("Wrong fdelta, try 'unorm', 'weighted', 'count' or 'mnorm' (recommended).")
    }       
  }
  
  # -------------------------------------------------------------------------    
  # Find delta
  # -------------------------------------------------------------------------       
  if(fdelta == "count"){
    f1 <- rank(f, ties.method = "first") # Break ties in f
    delta <- apply(distm / outer(f1, f1, FUN = ">"), 2, min, na.rm = TRUE)
    loc.max <- which.max(delta)
    delta[loc.max] <- max(delta[-loc.max]) # Equation in the Matlab code
  }else if(fdelta == "mnorm"){
    f.order <- order(f, decreasing = TRUE)
    delta <- rep(NA, n)
    delta[f.order[1]] <- Inf
    for(i in 2:length(f.order)){
      delta[f.order[i]] <- min(distm[f.order[i], f.order[1:(i - 1)]])
    }
    delta[f.order[1]] <- max(delta[-f.order[1]])
  }else{
    delta <- apply(distm / outer(f, f, FUN = ">"), 2, min, na.rm = TRUE)
    loc.max <- which.max(delta)
    delta[loc.max] <- max(delta[-loc.max]) # Equation in the Matlab code
  }
  
  return(list(f = f, delta = delta))
}

FindClustersAuto <- function(distm,
                             f, 
                             delta, 
                             ac = 1, 
                             nclust = 2:10, 
                             f.cut = c(0.1, 0.2, 0.3)){
  # -------------------------------------------------------------------------    
  # Check arguments
  # -------------------------------------------------------------------------           
  if(!inherits(distm, 'dist')) stop("arg distm must inherit dist class. Got ", class(distm))
  if(!is.numeric(f)) stop("arg f must be numeric. Got ", class(f))
  if(!is.numeric(delta)) stop("arg delta must be numeric. Got ", class(delta))
  if(attr(distm, 'Size') != length(f)) stop("length of f (", length(f),") not equal to number of observations in distm (", attr(distm, 'Size'), ")")
  if(attr(distm, 'Size') != length(delta)) stop("length of delta (", length(delta),") not equal to number of observations in distm (", attr(distm, 'Size'), ")")
  if(!all.equal(nclust, as.integer(nclust))) stop('arg nclust must all be integers. Got ', class(nclust))
  if(min(nclust) <= 0) stop('nclust must be positive integers')
  if(!is.numeric(f.cut)) stop('arg f.cut must be numeric. Got ', class(f.cut))
  if(length(f.cut) == 0) stop('arg f.cut is empty: ', f.cut)    
  if(min(f.cut) < 0) stop('arg f.cut must be between 0 - 1. Got', f.cut)
  if(max(f.cut) >= 1) stop('arg f.cut must be between 0 - 1. Got', f.cut)  
  if(length(ac) != 1) stop('arg ac must have length 1. Got', ac)
  
  # -------------------------------------------------------------------------
  # Find centers
  # -------------------------------------------------------------------------
  if(ac == 1){
    center.list <- FindCentersAutoV(f, delta, f.cut = f.cut, 
                                    nclust = nclust, rm.dup = FALSE)
  }else if(ac == 2){
    center.list <- FindCentersAutoD(f, delta, nclust = nclust)
  }else{
    stop("Wrong ac. Must be either 1 or 2. Got ", ac)
  }
  
  if(length(center.list) == 0){
    stop("Failed to find any centers")
  }
  
  # -------------------------------------------------------------------------
  # Cluster
  # -------------------------------------------------------------------------
  cluster.list <- lapply(center.list, 
                         function(x){
                           a <- FindClustersGivenCenters(distm, centers = x)
                           attributes(a) <- attributes(x)
                           return(a)
                         })
  sils.list <- lapply(cluster.list, 
                      function(x){
                        a <- FindSilhouette(distm, clusters = x)
                        attributes(a) <- attributes(x)
                        return(a)
                      })
  
  sils.vector <- unlist(sils.list)
  
  winner.i <- which.max(sils.vector)
  
  ans <- list()
  ans[['clusters']] <- cluster.list[[winner.i]]
  ans[['centers']] <- center.list[[winner.i]]
  ans[['silhouette']] <- sils.list[[winner.i]]
  ans[['nclust']] <- length(ans[['centers']])
  ans[['tested.sils']] <- sils.list
  return(ans)
}

FindCentersAutoV <- function(f, delta, f.cut = c(0.1, 0.2, 0.3), nclust, rm.dup = TRUE){
  # -------------------------------------------------------------------------
  # Check arguments
  # -------------------------------------------------------------------------    
  if(!is.numeric(f.cut)) stop('arg f.cut should inherit numeric. Got ', class(f.cut))
  if(length(f.cut) == 0) stop('arg f.cut is empty: ', f.cut)    
  if(min(f.cut) < 0) stop('arg f.cut must be between 0 - 1. Got', f.cut)
  if(max(f.cut) >= 1) stop('arg f.cut must be between 0 - 1. Got', f.cut)    
  if(!is.numeric(nclust)) stop('arg nclust should inherit numeric. Got ', class(nclust))
  if(!all.equal(nclust, as.integer(nclust))) stop('nclust must all be integers')
  if(min(nclust) <= 0) stop('nclust must be positive integers')    
  
  center.list <- list()
  
  for(i in seq_along(f.cut)){ # For each f.cuts
    ##f0 <- min(f) + f.cut[j] * (max(f) - min(f))
    f0 <- stats::quantile(f, probs = f.cut[i])
    delta1 <- delta
    delta1[f < f0] <- -Inf
    cts <- order(delta1, decreasing = TRUE)[1 : max(nclust)]
    for(j in seq_along(nclust)){ # For each nclust
      if(sum(f >= f0) < nclust[j]){ # Number of points that > f.cut is less than nclust. Stop
        stop("Only (", sum(f >= f0), ") points to the right of f.cut (", f0, "), but nclust = ", nclust[j])
      }
      centers <- cts[1 : nclust[j]]
      attributes(centers) <- list(f.cut = f.cut[i], 
                                  f.cut.value = f0,
                                  nclust = nclust[j])
      if(rm.dup){
        if(!IsDup(center.list, centers)){
          center.list <- c(center.list, list(centers))
        }
      }else{
        center.list <- c(center.list, list(centers))
      }
    }
  }
  return(center.list = center.list)
}

FindClustersGivenCenters <- function(distm, centers){
  if(!inherits(distm, 'dist')) stop("arg distm must be a dist class")
  if(!is.numeric(centers)) stop("arg centers must be numeric. Got ", class(centers))
  if(length(centers) == 0) stop("arg centers must have length > 1")
  if(!all(centers == floor(centers))) stop("arg centers must be a (vector of) integer(s). Got ", centers)
  if(anyDuplicated(centers)) stop("Duplications in centers not allowed")
  if(anyNA(centers)) stop("NA in centers not allowed")
  if(min(centers) <= 0) stop("min(centers) <= 0")
  if(max(centers) > attr(distm, 'Size')) stop("max(centers) larger than number of observations in distm")
  
  centers <- unique(centers)
  if(length(centers) <= 1) stop('length of unique(centers) must be greater than 1')
  distm <- as.matrix(distm)
  dist.to.centers <- distm[, centers]
  clusters <- apply(dist.to.centers, 1, FUN = which.min)
  return(clusters)
}

FindSilhouette <- function(distm, clusters){
  if(!inherits(distm, 'dist')) stop("arg distm must inherit dist class. Got: ", class(distm))
  if(!all(clusters == floor(clusters))) stop('arg clusters must all be integers')
  if(length(clusters) != attr(distm, 'Size')) stop('length of clusters', length(clusters),
                                                   'not equal to number of observations in distm', attr(distm, 'Size'))
  ans <- mean(cluster::silhouette(x = clusters, dist = distm)[,3])
  return(ans)
}

FindDistm <- function(x, normalize = FALSE, method = 'euclidean'){
  if(!inherits(x, 'data.frame') && !inherits(x, 'matrix')) stop('arg x must be data frame or matrix')
  if(nrow(x) == 0) stop('x is empty. Cannot calculate distance matrix.')
  if(!inherits(normalize, 'logical')) stop('arg normalize must be boolean')
  if(normalize){
    ## distm.std <- as.matrix(dist(scale(dat, center = FALSE, scale = TRUE),
    ##                             method = "euclidean"))
    sds <- apply(x, 2, stats::sd)
    distm <- stats::dist(scale(x, center = FALSE, scale = sds), method = method, upper = TRUE)
  }else{
    distm <- stats::dist(x, upper = TRUE, method = method)
  }
  return(distm)
}

adpclust <- function(x = NULL,
                     distm = NULL,
                     p = NULL,
                     centroids = 'auto',
                     h = NULL, 
                     htype = 'amise',
                     nclust = 2:10,
                     ac = 1,
                     f.cut = c(0.1, 0.2, 0.3),
                     fdelta = 'mnorm',
                     dmethod = 'euclidean',
                     draw = FALSE,
                     f = NULL
){
  # -------------------------------------------------------------------------
  # Check arguments
  # -------------------------------------------------------------------------    
  if(!centroids %in% c('user', 'auto')){
    stop('arg centroids must be one of c(\'user\', \'auto\') Got ', centroids)  
  }
  if(!is.null(distm) && !inherits(distm, 'dist')) stop("arg distm must inherit dist class. Got ", class(distm))
  if(!is.null(h)){
    if(!is.null(f)) stop('Both h and f are provided. If you provide the 
                         density f, then there is no need to specify h, 
                         since there is no need to use h to estimate f.')
    if(!is.numeric(h)) stop('arg h must be numeric. Got ', class(h))
    if(length(h) == 0) stop('arg h is empty: ', h)    
    if(min(h) <= 0) stop('arg h must be nonnegative. Got', h)
  }    
  if(!tolower(htype) %in% c('amise', 'rot')){
    stop('arg centroids must be one of c(\'amise\', \'rot\') Got ', htype)  
  }
  if(!all(nclust == floor(nclust))) stop('arg nclust must all be integers. Got ', nclust)
  if(min(nclust) <= 1) stop('arg nclust must be integers > 1. Got ', nclust)
  if(!ac %in% c(1,2)) stop('arg ac must be one of c(1,2). Got ', ac)
  if(!is.numeric(f.cut)) stop('arg f.cut must be numeric. Got ', class(f.cut))
  if(length(f.cut) == 0) stop('arg f.cut is empty: ', f.cut)
  if(min(f.cut) < 0) stop('arg f.cut must be between 0 - 1. Got', f.cut)
  if(max(f.cut) >= 1) stop('arg f.cut must be between 0 - 1. Got', f.cut)    
  if(!fdelta %in% c('mnorm', 'unorm', 'weighted', 'count')){
    stop('arg fdelta must be one of c(\'mnorm\', \'unorm\', \'weighted\', \'count\'). Got ', fdelta)        
  }
  if(!is.null(f)){
    if(!is.vector(f)) stop("f must be a vector, Got ", class(f))
    if(!is.null(x) && length(f) != nrow(x)){
      stop("Length of f and length of x don't match")
    } 
    if(!is.null(distm) && attr(distm, "Size") != length(f)){
      stop("Length of f and size of distm must match")
    }
  }
  
  if(is.null(x)){
    # Use distm
    if(is.null(distm)) stop("Must provide one of x or distm")
    if(!inherits(distm, 'dist')) stop("arg distm must inherit dist class. Got ", class(distm))
    if(is.null(p) && is.null(h) && is.null(f)){
      stop("Data x are not given. Densities values f are not given. 
           You must provide the bandwidth h for the model to calculate f, 
           or provide dimension p so we can estimate h.") 
    }
    }else{
      # Use x. Calculate distm.
      if(fdelta == "mnorm"){
        distm <- FindDistm(x, normalize = TRUE, method = "euclidean")
      }else{
        distm <- FindDistm(x, normalize = FALSE, method = dmethod)
      }
      p = ncol(x)
    }
  
  # -------------------------------------------------------------------------
  # Find bandwidth h
  # -------------------------------------------------------------------------    
  if(is.null(h) && is.null(f)){
    if(fdelta != "mnorm"){
      stop("Must give h unless fdelta == 'mnorm'")
    }
    h <- FindH(p, attr(distm, 'Size'), htype)
  }
  
  # -------------------------------------------------------------------------    
  # Clustering with the 'user' option
  # -------------------------------------------------------------------------    
  if(centroids == "user"){
    if(is.null(f) && (length(h) > 1)){
      stop("h must be a scalar when centroids == 'user'")
    }
    fd <- FindFD(distm = distm, h = h, fdelta = fdelta, f = f)
    ans <- FindClustersManual(distm, fd$f, fd$delta)
    ans[['h']] <- h
    ans[['f']] <- fd[['f']]
    ans[['delta']] <- fd[['delta']]
    ans[['selection.type']] <- 'user'
    class(ans) <- c("adpclust", "list")
    if(draw) plot.adpclust((ans))
    return(ans)
  }
  
  # -------------------------------------------------------------------------    
  # Clustering with the 'auto' option
  # -------------------------------------------------------------------------    
  if(centroids == "auto"){
    if(is.null(f)){
      if(length(h) > 1){
        h.seq <- h
      }else{
        h.seq <- seq(h / 3, h * 3, length.out = 10)
      }
      # Find f and delta for each h
      fd.list <- lapply(h.seq, function(h) FindFD(distm, h, fdelta))            
    }else{
      fd.list <- list(FindFD(distm = distm, h = h, fdelta = fdelta, f = f))
      h.seq <- list(NULL)
    }
    
    result.list <- lapply(fd.list, function(fd){
      FindClustersAuto(distm = distm, 
                       f = fd[['f']], 
                       delta = fd[['delta']], 
                       ac = ac,
                       nclust = nclust,
                       f.cut = f.cut)  
    })
    score.seq <- sapply(result.list, function(x) x$silhouette)
    iwinner <- which.max(score.seq)
    # Generate a list of all tested possibilities
    tested <- list()
    for(i in seq_along(h.seq)){
      for(one.sil in result.list[[i]]$tested.sils){
        one.tested <- list(f.cut = attr(one.sil, 'f.cut'),
                           f.cut.value = attr(one.sil, 'f.cut.value'),
                           nclust = attr(one.sil, 'nclust'),
                           h = h.seq[i],
                           sil = as.vector(one.sil))
        tested <- c(tested, list(one.tested))                
      }
    }
    
    ans <- result.list[[iwinner]]
    ans[['tested.sils']] <- NULL # Redundant. In 'tested'
    ans[['h']] <- h.seq[iwinner]
    fd <- fd.list[[iwinner]]
    ans[['f']] <- fd[['f']]
    ans[['delta']] <- fd[['delta']]
    ans[['selection.type']] <- 'auto'
    ans[['tested']] <- tested
    class(ans) <- c("adpclust", "list")
    if(draw) plot.adpclust((ans))
    return(ans)
  }
  stop('centroids not recognized') # should never reach here.
}





#load packages
library(R.utils)
library(stats)
library(Rtsne)


#list all files in this directory
#list<-list.files("cnnimpute/input/")

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)  
print(args)
mat1<-read.csv(args[1],header = TRUE,row.names=1)
  
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
file = sub("^.*/", "", args[1])
write.csv(file=paste0("cnnImpute/tmp_file/",file),mat1) 
  


