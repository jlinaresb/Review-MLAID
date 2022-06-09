# SNF
# ===
snf = function(data, K = 20, alpha = 0.5, iters = 30){
  
  require(SNFtool)
  df = standardNormalization(data)
  
  start = Sys.time()
  Dist = SNFtool::dist2(as.matrix(df), as.matrix(df))
  W = affinityMatrix(Dist, K, alpha)
  # W = SNF(list(W), K, t = iters)
  nclust = estimateNumberOfClustersGivenGraph(W, 2:20)[[1]]
  end = Sys.time()
  time = end - start
  
  return(list(time = time,
              nclust = nclust,
              labels = W))
  
}


# iClusterPlus
# ===
# icluster = function(){
#   
#   
#   
# }
