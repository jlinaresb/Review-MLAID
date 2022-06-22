# Preprocessing multi-omic
# ====
cesga = T

if (cesga == T) {
  setwd('/mnt/netapp2/Store_uni/home/ulc/co/jlb/git/Review-MLAID/')
} else {
  setwd('~/git/Review-MLAID/')
}

# Transcriptomics
# ===
# GSE117931
assay = fread('benchmark/extdata/multi-omic/GSE117931/expression.tsv',
              header = T, data.table = F)
rownames(assay) = assay[,1]; assay = assay[, -1]

meta = read.delim2(
  'benchmark/extdata/multi-omic/GSE117931/metadata.tsv',
  header = T)
meta = meta[match(colnames(assay), meta$Sample),]
rownames(meta) = meta$Sample

eset1 = ExpressionSet(
  assayData = as.matrix(assay),
  phenoData = AnnotatedDataFrame(meta)
)

eset1.f = genefilter::varFilter(eset1, var.cutoff = 0.7)

saveRDS(eset1, file = 'benchmark/data/multi-omic/eset_GSE117931_expression_complete.rds')
saveRDS(eset1.f, file = 'benchmark/data/multi-omic/eset_GSE117931_expression_filtered.rds')


# GSE82221
assay = fread('benchmark/extdata/multi-omic/GSE82221/expression.tsv',
              header = T, data.table = F)
rownames(assay) = assay[,1]; assay = assay[, -1]

meta = read.delim2(
  'benchmark/extdata/multi-omic/GSE82221/ADEx_data/metadata.tsv',
  header = T)
meta = meta[match(colnames(assay), meta$Sample),]
rownames(meta) = meta$Sample

eset1 = ExpressionSet(
  assayData = as.matrix(assay),
  phenoData = AnnotatedDataFrame(meta)
)

eset1.f = genefilter::varFilter(eset1, var.cutoff = 0.7)

saveRDS(eset1, file = 'benchmark/data/multi-omic/eset_GSE82221_expression_complete.rds')
saveRDS(eset1.f, file = 'benchmark/data/multi-omic/eset_GSE82221_expression_filtered.rds')



# Epigenomics
# ===
# GSE117931
assay = fread('benchmark/extdata/multi-omic/GSE117931/methylation.tsv',
              header = T, data.table = F)
rownames(assay) = assay[,1]; assay = assay[, -1]
assay = na.omit(assay)
assay = as.data.frame(t(apply(assay, 1, function(x) log2(x / (1-x)))))
meta = read.delim2(
  'benchmark/extdata/multi-omic/GSE117931/metadata.tsv',
  header = T)
meta = meta[match(colnames(assay), meta$Sample),]
rownames(meta) = meta$Sample

eset1 = ExpressionSet(
  assay = as.matrix(assay),
  phenoData = AnnotatedDataFrame(meta)
)

eset1.f = genefilter::varFilter(eset1, var.cutoff = 0.01)

saveRDS(eset1, file = 'benchmark/data/multi-omic/eset_GSE117931_methylation_complete.rds')
saveRDS(eset1.f, file = 'benchmark/data/multi-omic/eset_GSE117931_methylation_filtered.rds')



# GSE82221
assay = fread('benchmark/extdata/multi-omic/GSE82221/methylation.tsv',
              header = T, data.table = F)
rownames(assay) = assay[,1]; assay = assay[, -1]
assay = na.omit(assay)
assay = as.data.frame(t(apply(assay, 1, function(x) log2(x / (1-x)))))
meta = read.delim2(
  'benchmark/extdata/multi-omic/GSE82221/metadata.tsv',
  header = T)
meta = meta[match(colnames(assay), meta$Sample),]
rownames(meta) = meta$Sample

eset1 = ExpressionSet(
  assay = as.matrix(assay),
  phenoData = AnnotatedDataFrame(meta)
)

eset1.f = genefilter::varFilter(eset1, var.cutoff = 0.01)

saveRDS(eset1, file = 'benchmark/data/multi-omic/eset_GSE82221_methylation_complete.rds')
saveRDS(eset1.f, file = 'benchmark/data/multi-omic/eset_GSE82221_methylation_filtered.rds')

