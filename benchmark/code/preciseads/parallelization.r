source('/mnt/netapp2/Store_uni/home/ulc/co/jlb/git/Review-MLAID/benchmark/code/preciseads/config_file.r')  
#source('~/git/Review-MLAID/benchmark/code/preciseads/config_file.r')

input.paths.ex = dir(path = paste0(input.dir.path, 'exp/'))
input.paths.met = dir(path = paste0(input.dir.path, 'met/'))
input.algs = dir(path = path.algs, pattern = pattern.algs)

# input.paths = c(input.paths.ex, input.paths.met)
# i = 1
# j = 1

for (i in 1:length(input.paths.ex)) {
  for (j in 1:length(input.algs)) {
    
    exec.dir = paste0(base.dir, 'benchmark/code/preciseads/Exec/')
    if (dir.exists(exec.dir) == FALSE) {
      dir.create(exec.dir)
      message('Creating directory --Exec-- ...')
    }
    
    jobname = input.algs[j]
    jobname = gsub('run_', '', jobname)
    jobname = gsub('.r', '', jobname, fixed = T)
    
    fstype = input.paths.ex[i]
    fstype = gsub('.rds', '', unlist(lapply(strsplit(fstype, '_'), '[[', 2)), fixed = T)
    
    if (cluster.type == 'multi') {
      expPath = input.paths.ex[i]
      metPath = input.paths.met[i]
    } else if (cluster.type == 'single'){
      if (jobname == 'expression') {
        expPath = input.paths.ex[i]
        metPath = NULL
      } else if (jobname == 'methylation'){
        expPath = NULL
        metPath = input.paths.met[i]
      }
    }
    
    filename = paste0(exec.dir, cluster.type, '_', jobname, '_', fstype, '.sh')
    name = gsub('.sh', '', basename(filename), fixed = T)
    
    sink(filename)
    
    cat("#!/bin/bash \n")
    
    cat(paste("#SBATCH", "-p", part, '\n'))
    cat(paste("#SBATCH", "--qos", qos, '\n'))
    cat(paste("#SBATCH", "-t", time, '\n'))
    cat(paste("#SBATCH", "--mem", mem, '\n'))
    cat(paste("#SBATCH", "-n", ntasks, '\n'))
    cat(paste("#SBATCH", "--mail-user=joselb1993@gmail.com", '\n'))
    cat(paste("#SBATCH", "--mail-type=END", '\n'))
    
    cat("module load cesga/2018 gcc/6.4.0 R/4.0.2", '\n')
    
    cat(paste("expPath=", input.dir.path, 'exp/', expPath, '\n', sep = ''))
    cat(paste("metPath=", input.dir.path, 'met/', metPath, '\n', sep = ''))
    cat(paste("name=", name, '\n', sep = ''))
    
    cat(paste('Rscript ', path.algs, input.algs[j],  " $expPath", " $metPath", " $name", sep=''))
    
    sink(file = NULL)
    
    system(paste('sbatch ', filename))
  }
}
