grepCluster = function(data){
  
  # Arguments:
  # ===
  # data = data.frame
  #
  # Value:
  # ===
  # data.frame with clustering works 
  #   after quality analysis (disease in the abstract)

  # Check if 'cluster' or 'unsupervised' or 'stratify'
  x = grep('cluster', tolower(data$title))
  y = grep('unsupervised', tolower(data$title))
  z = grep('stratif', tolower(data$title))
  a = grep('subtype', tolower(data$title))
  data = data[unique(c(x, y, z, a)),]
  
  return(data)
}
