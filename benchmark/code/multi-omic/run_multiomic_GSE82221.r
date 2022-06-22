# Run GSE82221
# ===

cesga = F

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

# Run the models
# ===
run_multiomic(data1 = data1,
              data2 = data2,
              outPath = 'benchmark/results/multi-omic/',
              file = 'multiomic_GSE82221',
              return = F,
              save = T)




