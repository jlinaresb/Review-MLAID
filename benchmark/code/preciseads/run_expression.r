cesga = T

if (cesga == T) {
  setwd('/mnt/netapp2/Store_uni/home/ulc/co/jlb/git/Review-MLAID')
} else{
  setwd('~/git/Review-MLAID/')
}

source('benchmark/code/utils.r')
source('benchmark/code/single-omic/cluster_algorithms.r')

exp = readRDS('benchmark/data/preciseads/expression_filtered_levene.rds')

runall(data = exp,
       outPath = 'benchmark/results/preciseads/',
       file = 'expression_preciseads',
       return = F,
       save = T)