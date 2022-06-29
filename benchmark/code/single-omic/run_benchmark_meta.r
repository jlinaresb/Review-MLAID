# Benchmark in transcriptomic data
# ===
require(phyloseq)
cesga = T

if (cesga == T) {
  setwd('/mnt/netapp2/Store_uni/home/ulc/co/jlb/git/Review-MLAID')
} else {
  setwd('~/git/Review-MLAID/')
}

source('benchmark/code/single-omic/cluster_algorithms.r')

path = 'benchmark/data/single-omic/'
files = list.files(path, pattern = 'meta')
files = files[grep('filtered', files)]

data = list()
for (i in seq_along(files)) {
  dat = readRDS(paste0(path, files[i]))
  data[[i]] = otu_table(dat)@.Data
}
names(data) = gsub('.rds', '', files)

for (i in seq_along(data)) {
  runall(data[[i]],
         file = paste0('metagenomic_', names(data)[i]),
         outPath = 'benchmark/results/single-omic/',
         return = F, save = T)
}

