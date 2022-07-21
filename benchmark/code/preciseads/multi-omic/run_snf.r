# Run

cesga = F

if (cesga == T) {
  setwd('/mnt/netapp2/Store_uni/home/ulc/co/jlb/git/Review-MLAID')
} else{
  setwd('~/git/Review-MLAID/')
}

source('benchmark/code/utils.r')
source('benchmark/code/multi-omic/cluster_algorithms_multiomics.r')

exp = readRDS('benchmark/data/preciseads/expression_fcbf_noctrls.rds')
met = readRDS('benchmark/data/preciseads/methylation_fcbf_noctrls.rds')

res = snf(data1 = exp,
    data2 = met,
    K = 20,
    alpha = 0.5,
    iters = 30,
    outPath = 'benchmark/results/preciseads/',
    file = 'snf_preciseads_multiOmic',
    return = T,
    save = F)

res$nclust
table(res$labels, meta.noctrls$Clusters)
