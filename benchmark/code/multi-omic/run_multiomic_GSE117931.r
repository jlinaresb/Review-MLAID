# Run GSE117931
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
file1 = list.files(path, pattern = 'GSE117931_expression_filtered')
dat = readRDS(paste0(path, file1))
data1 = exprs(dat)

file2 = list.files(path, patter = 'GSE117931_methylation_filtered')
dat = readRDS(paste0(path, file2))
data2 = exprs(dat)

# Run the models with multi-omic
# ===
run_multiomic(data1 = data1,
              data2 = data2,
              outPath = 'benchmark/results/multi-omic/',
              file = 'multiomic_GSE117931',
              return = F,
              save = T)


# Run models with single-omic
# ===
data = list(
  expression = data1,
  methylation = data2
)

for (i in seq_along(data)) {
  runall(data[[i]],
         file = paste0('GSE117931_', names(data)[i]),
         outPath = 'benchmark/results/multi-omic/',
         return = F, save = T)
}
