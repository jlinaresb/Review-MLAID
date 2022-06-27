# Clustering algoirthms
# ====

# Utils
# ===
require(Biobase)
source('/mnt/netapp2/Store_uni/home/ulc/co/jlb/git/Review-MLAID/benchmark/code/utils.r')


# K-means clustering
# ===
KMEANS = function(data){
  
  require(NbClust)
  df = scale(data)
  
  start = Sys.time()
  res = NbClust(data = df, distance = 'euclidean',
                method = 'kmeans', index = 'tau', max.nc = 10)
  end = Sys.time()
  time = difftime(end, start, units = 'secs')
  
  nclust = length(unique(res$Best.partition))
  return(list(time = time,
              nclust = nclust,
              labels = res$Best.partition,
              fit = res))
  
}


# Hierarchical clustering
# ===
hierarchical = function(data){
  
  df = scale(data)
  
  start = Sys.time()
  res = NbClust::NbClust(data = df, method = 'ward.D', index = 'tau')
  end = Sys.time()
  time = difftime(end, start, units = 'secs')
  
  nclust = length(unique(res$Best.partition))
  return(list(time = time,
              nclust = nclust,
              labels = res$Best.partition,
              fit = res))
  
}


# SOM (?)
# ===
SOM = function(data){
  
  # Extracted from: https://www.shanelynn.ie/self-organising-maps-for-customer-segmentation-using-r/
  require(kohonen)
  df = scale(t(data))
  
  dim = mapping(nrow(df))
  
  start = Sys.time()
  grid = somgrid(xdim = dim, ydim = dim,
                          topo = "hexagonal")

  map = som(df, grid=grid)
  
  res = NbClust::NbClust(data = map$codes[[1]], method = 'ward.D', index = 'tau')
  nclust = length(unique(res$Best.partition))
  som_cluster = res$Best.partition
  cluster_assignment = som_cluster[map$unit.classif]
  data$cluster = cluster_assignment
  
  labels = data$cluster
  names(labels) = rownames(data)
  
  end = Sys.time()
  time = difftime(end, start, units = 'secs')
  
  nclust = length(unique(res$Best.partition))
  return(list(time = time,
              nclust = nclust,
              labels = labels,
              fit = res))
  
}


# Consensus of k-Means cluster
# ===
conCluster.km = function(data, maxK){
  
  require(ConsensusClusterPlus)
  
  data = scale(data)
  
  start = Sys.time()
  res = ConsensusClusterPlus(data,
                             maxK = maxK,
                             reps = 1000,
                             pItem = 0.8,
                             pFeature = 1,
                             clusterAlg = 'km',
                             distance = 'euclidean',
                             seed = 1993,
                             plot = NULL)

  end = Sys.time()
  time = difftime(end, start, units = 'secs')
  
  nclust = select_k(res, maxK)
  labels = as.numeric(as.factor(res[[nclust]]$clrs[[1]]))
  names(labels) = colnames(data)
  
  return(list(time = time,
              nclust = nclust,
              labels = labels,
              fit = res))
  
}



# Consensus of Hierarchical cluster
# ===
conCluster.h = function(data, maxK){
  
  require(ConsensusClusterPlus)

  data = scale(data)
  
  start = Sys.time()
  res = ConsensusClusterPlus(data,
                             maxK = maxK,
                             reps = 1000,
                             pItem = 0.8,
                             pFeature = 1,
                             clusterAlg = 'hc',
                             distance = 'pearson',
                             seed = 1993,
                             plot = NULL)

  end = Sys.time()
  time = difftime(end, start, units = 'secs')
  
  nclust = select_k(res, maxK)
  labels = as.numeric(as.factor(res[[nclust]]$clrs[[1]]))
  names(labels) = colnames(data)
  
  return(list(time = time,
              nclust = nclust,
              labels = labels,
              fit = res))
  
}



# NMF
# ===
NMF = function(data){
  
  require(NMF)
  
  normalize <- function(x, na.rm = TRUE) {
   return((x- min(x)) /(max(x)-min(x)))
  }
  
  data = as.data.frame(t(apply(data, 1, function(x) normalize(x))))

  # meth <- nmfAlgorithm(version='R')
  # meth <- c(names(meth), meth)
  
  data = t(data)
  
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
mclust = function(data){
  
  require(mclust)
  data = t(data)
  
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




# Run all
# ===
runall = function(data, file, outPath, return = F, save = T){
  
  # Arguments
  #   dim: for SOM map
 
  print(paste0('Runing hierarchical consensus clustering for: ', file))
  cc.h = conCluster.h(data, 10)

  print(paste0('Runing K-Means consensus clustering for: ', file))
  cc.km = conCluster.km(data, 10)
  
  print(paste0('Runing SOM clustering for: ', file))
  som = SOM(data)

  #print(paste0('Runing hierarchical clustering for: ', file))
  #h = hierarchical(data)  
  
  #print(paste0('Runing K-Means clustering for: ', file))
  # k = KMEANS(data)        
  
  print(paste0('Runing Mclust clustering for: ', file))
  m = mclust(data)
  
  print(paste0('Runing NMF clustering for: ', file))
  nmf = NMF(data)
  
  #print(paste0('Runing SNF clustering for: ', file))
  #snf = snf(data)           
  
  res = list(
            Consensus.H = cc.h,
            Consensus.KM = cc.km,
            SOM = som,
            MClust = m,
            NMF = nmf
            )

  if (save == T){
    saveRDS(res, file = paste0(outPath, file, '.rds'))
  }
   
  if (return == T) {
    return(res)
  }
  
}

