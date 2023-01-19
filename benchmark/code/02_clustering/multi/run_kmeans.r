# Run

cesga = T

if (cesga == T) {
  setwd('/mnt/netapp2/Store_uni/home/ulc/co/jlb/git/Review-MLAID')
} else{
  setwd('~/git/Review-MLAID/')
}

source('benchmark/code/utils.r')
source('benchmark/code/multi-omic/cluster_algorithms_multiomics.r')

args = commandArgs(trailingOnly = T)

exp = readRDS(args[1])
met = readRDS(args[2])

conCluster.km(data1 = exp,
              data2 = met,
              maxK = 10,
              outPath = 'benchmark/results/preciseads/',
              file = args[3],
              return = F,
              save = T)