cesga = T

if (cesga == T) {
  setwd('/mnt/netapp2/Store_uni/home/ulc/co/jlb/git/Review-MLAID')
} else{
  setwd('~/git/Review-MLAID/')
}

source('benchmark/code/utils.r')
source('benchmark/code/single-omic/cluster_algorithms.r')

args = commandArgs(trailingOnly = T)

met = readRDS(args[2])

runall(data = met,
       outPath = 'benchmark/results/preciseads/',
       file = args[3],
       return = F,
       save = T)