# Run GSE82221
# ===

cesga = T

if (cesga == T) {
  setwd('/mnt/netapp2/Store_uni/home/ulc/co/jlb/git/Review-MLAID')
} else {
  setwd('~/git/Review-MLAID/')
}


source('benchmark/code/multi-omic/cluster_algorithms_multiomics.r')



# Load data
# ===
path = 'benchmark/data/multi-omic/'
file1 = list.files(path, pattern = 'GSE82221_expression_filtered')
dat = readRDS(paste0(path, file1))
data1 = exprs(dat)

file2 = list.files(path, pattern = 'GSE82221_methylation_filtered')
dat = readRDS(paste0(path, file2))
data2 = exprs(dat)

# Run the models with multi-omic
# ===
run_multiomic(data1 = data1,
              data2 = data2,
              outPath = 'benchmark/results/multi-omic/',
              file = 'multiomic_GSE82221',
              return = F,
              save = T)



# Run models with single-omic
# ===
data = list(
  expression = data1,
  methylation = data2
)
source('benchmark/code/single-omic/cluster_algorithms.r')
for (i in seq_along(data)) {
  runall(data[[i]],
         file = paste0('GSE82221_', names(data)[i]),
         outPath = 'benchmark/results/multi-omic/',
         return = F, save = T)
}
