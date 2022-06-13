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
met = read.table('benchmark/extdata/multi-omic/GSE117931/ADEx_data/methylation.tsv',
                 header = T, row.names = 1)
met = met[complete.cases(met),]
require(lumi)
metm = beta2m(met)   #change B values to M values



# Remove these lines!!!
exp = exp[1:100,]
met = met[1:100,]
######

# Run the models
# ===
run_multiomic(data1 = exp,
              data2 = met,
              outPath = '',
              return = T,
              save = F)


