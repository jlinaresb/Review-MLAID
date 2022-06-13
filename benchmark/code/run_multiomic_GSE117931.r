# Run GSE117931
# ===

cesga = F

if (cesga == T) {
  setwd('/mnt/netapp2/Store_uni/home/ulc/co/jlb/git/Review-MLAID')
} else {
  setwd('~/git/Review-MLAID/')
}

source('benchmark/code/cluster_algorithms_multiomics.r')


# Load data
exp = read.table('benchmark/extdata/multi-omic/GSE117931/ADEx_data/expression.tsv',
                 header = T, row.names = 1)
exp = exp[complete.cases(exp),]

met = readRDS('benchmark/extdata/multi-omic/GSE117931/ADEx_data/methylation_M.rds')
met = met[complete.cases(met),]


# Remove these lines!!!
exp = exp[1:100,]
met = met[1:100,]
######

# Run the models
# ===
run_multiomic(data1 = exp,
              data2 = met,
              outPath = NULL,
              return = T,
              save = F)


