# Real datasets
# ===

cesga = T

if (cesga == T) {
  setwd('/mnt/netapp2/Store_uni/home/ulc/co/jlb/git/Review-MLAID')
} else {
  setwd('~/git/Review-MLAID/')
}

source('benchmark/code/cluster_algorithms.r')

exp_path = 'benchmark/extdata/single-omic/transcriptomic/filtered/'             # Transcriptomic
epi_path = 'benchmark/extdata/single-omic/epigenomic/filtered/'                 # Epigenomic
meta_path = 'benchmark/extdata/single-omic/metagenomic/filtered/'               # Metagenomic

# Benchmark in metagenomic data
# ===
files = list.files(path = meta_path)
data= list()
for (i in seq_along(files)) {
  dat = readRDS(paste0(meta_path, files[i]))
  dat = subset(dat, select = -c(Condition))
  dat = dat[, which(colSums(dat) != 0)]
  data[[i]] = dat
}
names(data) = gsub('.rds', '', files)

for (i in seq_along(data)) {
  runall(data[[i]],
         file = paste0('metagenomic_', names(data)[i]),
         outPath = 'benchmark/results/',
         return = F, save = T)
}
