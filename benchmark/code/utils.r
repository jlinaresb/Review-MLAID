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
