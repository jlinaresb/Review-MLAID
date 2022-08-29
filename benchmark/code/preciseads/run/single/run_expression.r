cesga = T

if (cesga == T) {
  setwd('/mnt/netapp2/Store_uni/home/ulc/co/jlb/git/Review-MLAID')
} else{
  setwd('~/git/Review-MLAID/')
}

source('benchmark/code/utils.r')
source('benchmark/code/single-omic/cluster_algorithms.r')

args = commandArgs(trailingOnly = T)

exp = readRDS(args[1])


runall(data = exp,
       outPath = 'benchmark/results/preciseads/',
       file = args[3],
       return = F,
       save = T)