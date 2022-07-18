# Utils
# ===
require(Biobase)
source('/mnt/netapp2/Store_uni/home/ulc/co/jlb/git/Review-MLAID/benchmark/code/utils.r')


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


# SOM (?)
# ===
SOM = function(data1, data2){
  
  # Extracted from: https://www.shanelynn.ie/self-organising-maps-for-customer-segmentation-using-r/
  require(kohonen)
  stopifnot(colnames(data1) == colnames(data2))
  data_bind = rbind.data.frame(data1, data2)
  
  df = scale(t(data_bind))
  
  dim = mapping(nrow(df))
  
  start = Sys.time()
  grid = somgrid(xdim = dim, ydim = dim,
                 topo = "hexagonal")
  
  map = som(df, grid=grid)
  
  res = NbClust::NbClust(data = map$codes[[1]], method = 'ward.D', index = 'tau')
  nclust = length(unique(res$Best.partition))
  som_cluster = res$Best.partition
  cluster_assignment = som_cluster[map$unit.classif]
  # data$cluster = cluster_assignment
  
  labels = cluster_assignment
  names(labels) = colnames(data_bind)
  
  end = Sys.time()
  time = difftime(end, start, units = 'secs')
  
  nclust = length(unique(res$Best.partition))
  return(list(time = time,
              nclust = nclust,
              labels = labels,
              fit = res))
  
}



# NMF
# ===
NMF = function(data1, data2){
  
  require(NMF)
  stopifnot(colnames(data1) == colnames(data2))
  data_bind = rbind.data.frame(data1, data2)
  
  normalize <- function(x, na.rm = TRUE) {
    return((x- min(x)) /(max(x)-min(x)))
  }
  
  data = as.data.frame(t(apply(data_bind, 1, function(x) normalize(x))))
  
  # meth <- nmfAlgorithm(version='R')
  # meth <- c(names(meth), meth)
  
  start = Sys.time()
  res = nmf(data, 
            rank = c(2:10),
            method = 'brunet',
            nrun = 30,
            .options= 'vp15')
  
  end = Sys.time()
  time = difftime(end, start, units = 'secs')
  
  # plot(res)
  nclust = res$measures$rank[which.max(res$measures$silhouette.consensus)] # confirm!!! ¿max?
  
  h = coef(res$fit[[nclust-1]])
  labels = apply(h, 2, which.max)
  
  return(list(time = time,
              nclust = nclust,
              labels = labels,
              fit = res))
  
}

# Mclust
# ===
mclust = function(data1, data2){
  
  require(mclust)
  stopifnot(colnames(data1) == colnames(data2))
  data_bind = rbind.data.frame(data1, data2)
  
  data = t(data_bind)
  
  start = Sys.time()
  res = mclust::Mclust(data,
                       modelNames = 'EII',
                       na.rm = T)
  end = Sys.time()
  time = difftime(end, start, units = 'secs')
  nclust = res$G
  
  return(list(time = time,
              nclust = nclust,
              labels = res$classification,
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
                                dt1 = t(data1),
                                dt2 = t(data2),
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

COCA = function(data1, data2){
  
  require(coca)
  
  data = list(data1, data2)
  start = Sys.time()
  outputBuildMOC = coca::buildMOC(data,
                                  M = 2,
                                  maxK = 10,
                                  methods = 'hclust',
                                  distances = 'minkowski',
                                  fill = TRUE)
  
  moc = outputBuildMOC$moc
  res = coca::coca(moc, maxK = 10, hclustMethod = 'average')
  end = Sys.time()
  time = difftime(end, start, units = 'secs')
  
  return(list(time = time,
              fit = res))
  
  
}

# Run all

run_multiomic = function(data1, data2, outPath, file, return = F, save = T){
  
  
  print('Runing hierarchical consensus clustering')
  cc.h = conCluster.h(data1, data2, maxK = 10)
  
  print('Runing K-Means consensus clustering')
  cc.km = conCluster.km(data1, data2, maxK = 10)
  
  print('Runing SOM clustering')
  som = SOM(data1, data2)
  
  print('Runing NMF clustering')
  nmf = NMF(data1, data2)
  
  print('Runing mclust clustering')
  m = mclust(data1, data2)
  
  print('Runing SNF clustering')
  snf = snf(data1, data2)
  
  print('Runing iClusterPlust clustering')
  icluster = icluster(data1, data2)
  
  print('Runing COCA clustering')
  coca = COCA(data1, data2)
  
  res = list(
    Consensus.H = cc.h,
    Consensus.KM = cc.km,
    SOM = som,
    MClust = m,
    NMF = nmf,
    SNF = snf,
    iCluster = icluster,
    coca = coca
  )
  
  if (save == T){
    saveRDS(res, file = paste0(outPath, file, '.rds'))
  }
  
  if (return == T) {
    return(res)
  }
  
}



