# Run GSE82221
# ===

cesga = T

if (cesga == T) {
  setwd('/mnt/netapp2/Store_uni/home/ulc/co/jlb/git/Review-MLAID')
} else {
  setwd('~/git/Review-MLAID/')
}

source('benchmark/code/cluster_algorithms_multiomics.r')


# Load data
exp = readRDS('benchmark/extdata/multi-omic/GSE82221/filtered/expression.rds')
exp = exp[, complete.cases(exp)]
exp = as.data.frame(t(exp))

met = readRDS('benchmark/extdata/multi-omic/GSE82221/filtered/methylation.rds')
met = met[, complete.cases(met)]
met = as.data.frame(t(met))

# Remove these lines!!!
# exp = exp[1:100,]
# met = met[1:100,]
######

# Run the models
# ===
run_multiomic(data1 = exp,
              data2 = met,
              outPath = NULL,
              return = F,
              save = T)


