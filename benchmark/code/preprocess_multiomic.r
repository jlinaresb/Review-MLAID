# Preprocessing multi-omic
# ====
source('~/git/Review-MLAID/benchmark/code/preprocessing_functions.r')
setwd("~/git/Review-MLAID/benchmark/extdata/multi-omic/")

# Transcriptomics
# ===
# GSE117931
data = data.table::fread('GSE117931/ADEx_data/expression.tsv', header = T)
data = as.data.frame(data)
data = data[, -1]

pats = colnames(data)
meta = read.delim2('GSE117931/ADEx_data/metadata.tsv', header = T)
meta = meta[match(pats, meta$Sample), ]

data = as.data.frame(t(data))
data$Condition = meta$Condition

res = levene(data)
outpath = 'GSE117931/filtered/expression.rds'
saveRDS(res, file = outpath)

# GSE82221
data = data.table::fread('GSE82221/ADEx_data/expression.tsv', header = T)
data = as.data.frame(data)
data = data[, -1]

pats = colnames(data)
meta = read.delim2('GSE82221/ADEx_data/metadata.tsv', header = T)
meta = meta[match(pats, meta$Sample), ]

data = as.data.frame(t(data))
data$Condition = meta$Condition

res = levene(data)
outpath = 'GSE82221/filtered/expression.rds'
saveRDS(res, file = outpath)



# Epigenomics
# ===
# GSE117931
data = readRDS('GSE117931/ADEx_data/methylation_M.rds')

pats = colnames(data)
meta = read.delim2('GSE117931/ADEx_data/metadata.tsv', header = T)
meta = meta[match(pats, meta$Sample), ]

data = as.data.frame(t(data))
data$Condition = meta$Condition

y = data$Condition
x = subset(data, select = -c(Condition))
x = as.matrix(x)

data = na.delete(x)
data = remove.constant(data)
res = MostVariables(data, ntopGenes = 1000)
res = as.data.frame(res)

outpath = 'filtered/'
outpath = 'GSE117931/filtered/methylation.rds'
saveRDS(res, file = outpath)



# GSE82221
data = readRDS('GSE82221/ADEx_data/methylation_M.rds')

pats = colnames(data)
meta = read.delim2('GSE82221/ADEx_data/metadata.tsv', header = T)
meta = meta[match(pats, meta$Sample), ]

data = as.data.frame(t(data))
data$Condition = meta$Condition

y = data$Condition
x = subset(data, select = -c(Condition))
x = as.matrix(x)

data = na.delete(x)
data = remove.constant(data)
res = MostVariables(data, ntopGenes = 1000)
res = as.data.frame(res)

outpath = 'filtered/'
outpath = 'GSE82221/filtered/methylation.rds'
saveRDS(res, file = outpath)
