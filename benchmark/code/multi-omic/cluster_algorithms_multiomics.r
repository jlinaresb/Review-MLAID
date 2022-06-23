# Utils
# ===
require(Biobase)

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



# Consensus of k-Means cluster
# ===
conCluster.km = function(data1 ,data2, maxK = 10){
  
  require(ConsensusClusterPlus)
  stopifnot(colnames(data1) == colnames(data2))
  
  data_bind = rbind.data.frame(data1, data2)
  data_bind = scale(data_bind)
  
  # Fit model
  start = Sys.time()
  res = ConsensusClusterPlus(data_bind,
                             maxK = maxK,
                             reps = 1000,
                             pItem = 0.8,
                             pFeature = 1,
                             clusterAlg = 'km',
                             distance = 'euclidean',
                             seed = 1993,
                             plot = NULL)
  # res = calcICL(res)
  # res = calcICL(res,title="example")
  end = Sys.time()
  time = difftime(end, start, units = 'secs')
  
  nclust = select_k(res, maxK)
  labels = as.numeric(as.factor(res[[nclust]]$clrs[[1]]))
  names(labels) = colnames(data_bind)
  
  return(list(time = time,
              nclust = nclust,
              labels = labels,
              fit = res))
  
}


# Consensus of Hierarchical cluster
# ===
conCluster.h = function(data1, data2, maxK = 10){
  
  require(ConsensusClusterPlus)
  stopifnot(colnames(data1) == colnames(data2))
  
  data_bind = rbind.data.frame(data1, data2)
  data_bind = scale(data_bind)
  
  # Fit model
  start = Sys.time()
  res = ConsensusClusterPlus(data_bind,
                             maxK = maxK,
                             reps = 1000,
                             pItem = 0.8,
                             pFeature = 1,
                             clusterAlg = 'hc',
                             distance = 'pearson',
                             seed = 1993,
                             plot = NULL)
  # res = calcICL(res)
  # res = calcICL(res,title="example")
  end = Sys.time()
  time = difftime(end, start, units = 'secs')
  
  nclust = select_k(res, maxK)
  labels = as.numeric(as.factor(res[[nclust]]$clrs[[1]]))
  names(labels) = colnames(data_bind)
  
  return(list(time = time,
              nclust = nclust,
              labels = labels,
              fit = res))
  
}



# SNF
# ===
snf = function(data1, data2, K = 20, alpha = 0.5, iters = 30){
  
  data1 = t(data1)
  data2 = t(data2)
  
  require(SNFtool)
  d1 = standardNormalization(data1)
  d2 = standardNormalization(data2)
  
  start = Sys.time()
  Dist1 = SNFtool::dist2(as.matrix(d1), as.matrix(d1))
  Dist2 = SNFtool::dist2(as.matrix(d2), as.matrix(d2))
  
  W1 = affinityMatrix(Dist1, K, alpha)
  W2 = affinityMatrix(Dist2, K, alpha)
  
  W = SNF(list(W1, W2), K, t = iters)
  nclust = estimateNumberOfClustersGivenGraph(W, 2:10)[[1]]
  end = Sys.time()
  time = end - start
  
  return(list(time = time,
              nclust = nclust,
              labels = W))
}


# iClusterPlus
# ===
require(iClusterPlus)
icluster = function(data1, data2){
  
  # Fitting
  # ===
  cv = list()
  start = Sys.time()
  for(k in 1:10){
    cv[[k]] = tune.iClusterPlus(cpus = 15,
                                dt1 = t(exp),
                                dt2 = t(met),
                                type = c('gaussian','gaussian'),
                                K = k,
                                n.lambda = 55,
                                scale.lambda = c(1,1),
                                maxiter = 20)
  }
  end = Sys.time()
  time = difftime(end, start, units = 'secs')
  
  return(list(time = time,
              fit = cv))
  
}



# Run all

run_multiomic = function(data1, data2, outPath, file, return = F, save = T){
  
  
  print('Runing hierarchical consensus clustering')
  cc.h = conCluster.h(data1, data2, maxK = 10)
  
  print('Runing K-Means consensus clustering')
  cc.km = conCluster.km(data1, data2, maxK = 10)
  
  print('Runing SNF clustering')
  snf = snf(data1, data2)
  
  print('Runing iClusterPlust clustering')
  icluster = icluster(data1, data2)
  
  res = list(
    Consensus.H = cc.h,
    Consensus.KM = cc.km,
    SNF = snf,
    iCluster = icluster
  )
  
  if (save == T){
    saveRDS(res, file = paste0(outPath, file, '.rds'))
  }
  
  if (return == T) {
    return(res)
  }
  
}



