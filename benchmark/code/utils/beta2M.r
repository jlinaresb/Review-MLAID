setwd('~/git/Review-MLAID/')

# GSE117931
met = read.table('benchmark/extdata/multi-omic/GSE117931/ADEx_data/methylation.tsv',
                 header = T, row.names = 1)
met = met[complete.cases(met),]
require(lumi)
metm = beta2m(met)   #change B values to M values
saveRDS(metm, file = 'benchmark/extdata/multi-omic/GSE117931/ADEx_data/methylation_M.rds')

# GSE82221
met = read.table('benchmark/extdata/multi-omic/GSE82221/ADEx_data/methylation.tsv',
                 header = T, row.names = 1)
met = met[complete.cases(met),]
require(lumi)
metm = beta2m(met)   #change B values to M values
saveRDS(metm, file = 'benchmark/extdata/multi-omic/GSE82221/ADEx_data/methylation_M.rds')
