cesga = T

if (cesga == T) {
  setwd('/mnt/netapp2/Store_uni/home/ulc/co/jlb/git/Review-MLAID')
} else{
  setwd('~/git/Review-MLAID/')
}

source('benchmark/code/utils.r')
source('benchmark/code/single-omic/cluster_algorithms.r')

met = readRDS('benchmark/data/preciseads/methylation_filtered_levene.rds')

runall(data = met,
       outPath = 'benchmark/results/preciseads/',
       file = 'methylation_preciseads',
       return = F,
       save = T)