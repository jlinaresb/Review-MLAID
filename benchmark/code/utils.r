# Utils
# ===

select_k = function(res, maxK){
  
  # Extracted from: https://www.biostars.org/p/198789/
  Kvec = 2:maxK
  x1 = 0.1; x2 = 0.9 # threshold defining the intermediate sub-interval
  PAC = rep(NA,length(Kvec)) 
  names(PAC) = paste("K=",Kvec,sep="") # from 2 to maxK
  for(i in Kvec){
    M = res[[i]]$consensusMatrix
    Fn = ecdf(M[lower.tri(M)])
    PAC[i-1] = Fn(x2) - Fn(x1)
  }
  # The optimal K
  optK = Kvec[which.min(PAC)]
  
  return(optK)
  
}

mapping = function(rows){
  
  i = 1
  while(i^2 < rows){
    i = i + 1
  }
  
  return(i - 1)
}


select_k_icluster = function(ic, ids, plot = T){
  
  require(iClusterPlus)
  
  BIC = getBIC(ic)
  devR = getDevR(ic)
  minBICid = apply(BIC, 2, which.min)
  nK = length(ic)
  devRatMinBIC = rep(NA, nK)
  for (i in 1:nK) {
    devRatMinBIC[i] = devR[minBICid[i], i]
  }
  
  clusters = getClusters(ic)
  rownames(clusters)=ids
  colnames(clusters)=paste("K=",1:(length(ic)),sep="")
  k = which.max(devRatMinBIC)
  
  labels = clusters[,k]
  nclust = length(unique(labels))
  
  if (plot == T) {
    plot(1:(nK+1),c(0,devRatMinBIC),type="b",xlab="Number of clusters (K+1)", ylab="%Explained Variation")
  }
  
  return(list(nclust = nclust,
              labels = labels))
  
}


fast.cor.FS = function(data, thresh){
  
  stopifnot('target' %in% names(data))
  require(FCBF)
  
  y = as.factor(data$target)
  x = subset(data, select = - c(target))
  
  dis = discretize_exprs(t(x))
  # su_plot(dis, y)
  
  fcbf = fcbf(dis, y, verbose = T, minimum_su = thresh)
  xx = x[,fcbf$index]
  
  xx = as.data.frame(cbind(xx, target = y))
  
  return(xx)
  
}
