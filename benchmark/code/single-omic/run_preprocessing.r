cesga = T

if (cesga == T) {
  setwd('/mnt/netapp2/Store_uni/home/ulc/co/jlb/git/Review-MLAID/')
} else {
  setwd('~/git/Review-MLAID/')
}

# Convert transcriptomic to eset object
require(SummarizedExperiment)
require(data.table)

# Transcriptomics
# ====
# GSE11907_T1D_GPL96
assay = fread('benchmark/extdata/single-omic/transcriptomic/GSE11907_T1D_GPL96.tsv',
              header = T, data.table = F)
rownames(assay) = assay[,1]; assay = assay[, -1]

meta = read.delim2(
  'benchmark/extdata/single-omic/transcriptomic/metadata.tsv',
  header = T)
meta = meta[match(colnames(assay), meta$Sample),]
rownames(meta) = meta$Sample

eset1 = ExpressionSet(
  assayData = as.matrix(assay),
  phenoData = AnnotatedDataFrame(meta)
)

eset1.f = genefilter::varFilter(eset1, var.cutoff = 0.7)

saveRDS(eset1, file = 'benchmark/data/single-omic/exp_GSE11907_T1D_GPL96_complete.rds')
saveRDS(eset1.f, file = 'benchmark/data/single-omic/exp_GSE11907_T1D_GPL96_filtered.rds')


# GSE45291_RA
assay = fread('benchmark/extdata/single-omic/transcriptomic/GSE45291_RA.tsv',
              header = T, data.table = F)
rownames(assay) = assay[,1]; assay = assay[, -1]

meta = read.delim2(
  'benchmark/extdata/single-omic/transcriptomic/metadata.tsv',
  header = T)
meta = meta[match(colnames(assay), meta$Sample),]
rownames(meta) = meta$Sample

eset2 = ExpressionSet(
  assayData = as.matrix(assay),
  phenoData = AnnotatedDataFrame(meta)
)

eset2.f = genefilter::varFilter(eset2, var.cutoff = 0.7)

saveRDS(eset2, file = 'benchmark/data/single-omic/exp_GSE45291_RA_complete.rds')
saveRDS(eset2.f, file = 'benchmark/data/single-omic/exp_GSE45291_RA_filtered.rds')



# GSE61635
assay = fread('benchmark/extdata/single-omic/transcriptomic/GSE61635.tsv',
              header = T, data.table = F)
rownames(assay) = assay[,1]; assay = assay[, -1]

meta = read.delim2(
  'benchmark/extdata/single-omic/transcriptomic/metadata.tsv',
  header = T)
meta = meta[match(colnames(assay), meta$Sample),]
rownames(meta) = meta$Sample

eset3 = ExpressionSet(
  assayData = as.matrix(assay),
  phenoData = AnnotatedDataFrame(meta)
)

eset3.f = genefilter::varFilter(eset3, var.cutoff = 0.7)

saveRDS(eset3, file = 'benchmark/data/single-omic/exp_GSE61635_complete.rds')
saveRDS(eset3.f, file = 'benchmark/data/single-omic/exp_GSE61635_filtered.rds')






# Epigenomics
# =====
# GSE42861
assay = fread('benchmark/extdata/single-omic/epigenomic/GSE42861.tsv',
              header = T, data.table = F)
rownames(assay) = assay[,1]; assay = assay[, -1]
assay = na.omit(assay)
assay = as.data.frame(t(apply(assay, 1, function(x) log2(x / (1-x)))))
meta = read.delim2(
  'benchmark/extdata/single-omic/epigenomic/metadata.tsv',
  header = T)
meta = meta[match(colnames(assay), meta$Sample),]
rownames(meta) = meta$Sample

eset1 = ExpressionSet(
  assay = as.matrix(assay),
  phenoData = AnnotatedDataFrame(meta)
)

eset1.f = genefilter::varFilter(eset1, var.cutoff = 0.01)

saveRDS(eset1, file = 'benchmark/data/single-omic/epi_GSE42861_complete.rds')
saveRDS(eset1.f, file = 'benchmark/data/single-omic/epi_GSE42861_filtered.rds')



# GSE56606_Monocytes
assay = fread('benchmark/extdata/single-omic/epigenomic/GSE56606_Monocytes.tsv',
              header = T, data.table = F)
rownames(assay) = assay[,1]; assay = assay[, -1]
assay = na.omit(assay)
assay = as.data.frame(t(apply(assay, 1, function(x) log2(x / (1-x)))))

meta = read.delim2(
  'benchmark/extdata/single-omic/epigenomic/metadata.tsv',
  header = T)
meta = meta[match(colnames(assay), meta$Sample),]
rownames(meta) = meta$Sample

eset2 = ExpressionSet(
  assay = as.matrix(assay),
  phenoData = AnnotatedDataFrame(meta)
)

eset2.f = genefilter::varFilter(eset2, var.cutoff = 0.01)

saveRDS(eset2, file = 'benchmark/data/single-omic/epi_GSE56606_Monocytes_complete.rds')
saveRDS(eset2.f, file = 'benchmark/data/single-omic/epi_GSE56606_Monocytes_filtered.rds')


# GSE59250_B_cells
assay = fread('benchmark/extdata/single-omic/epigenomic/GSE59250_B_cells.tsv',
              header = T, data.table = F)
rownames(assay) = assay[,1]; assay = assay[, -1]
assay = na.omit(assay)
assay = as.data.frame(t(apply(assay, 1, function(x) log2(x / (1-x)))))

meta = read.delim2(
  'benchmark/extdata/single-omic/epigenomic/metadata.tsv',
  header = T)
meta = meta[match(colnames(assay), meta$Sample),]
rownames(meta) = meta$Sample

eset3 = ExpressionSet(
  assay = as.matrix(assay),
  phenoData = AnnotatedDataFrame(meta)
)

eset3.f = genefilter::varFilter(eset3, var.cutoff = 0.01)

saveRDS(eset3, file = 'benchmark/data/single-omic/epi_GSE59250_B_cells_complete.rds')
saveRDS(eset3.f, file = 'benchmark/data/single-omic/epi_GSE59250_B_cells_filtered.rds')



# Metagenomics
# ======
require(phyloseq)
require(genefilter)
flist = filterfun(kOverA(5, 2e-05))

# Morgan
otuPath = 'benchmark/extdata/single-omic/metagenomic/morgan/otutable.txt'
clinPath = 'benchmark/extdata/single-omic/metagenomic/morgan/clinical.txt'
taxPath = 'benchmark/extdata/single-omic/metagenomic/morgan/TaxonomyTable.rds'

otu = read.delim2(otuPath, header = T, sep = '\t')
clin = read.delim2(clinPath, header = T, sep = '\t')
rownames(clin) = paste0('X', clin$X.SampleID)
tax = readRDS(taxPath)

TAXA = tax_table(as.matrix(tax))
otu = subset(otu, select = -c(X.OTU.ID))
rownames(otu) = rownames(tax)
otu = as.matrix(otu)
OTU = otu_table(otu, taxa_are_rows = T)

CLIN = sample_data(clin)

morgan = phyloseq(OTU, CLIN, TAXA)
morgan = tax_glom(morgan, 'Rank6')

morgan.f = filter_taxa(morgan, flist, prune=TRUE)

saveRDS(morgan, file = 'benchmark/data/single-omic/meta_morgan_complete.rds')
saveRDS(morgan.f, file = 'benchmark/data/single-omic/meta_morgan_filtered.rds')



# Gevers
otuPath = 'benchmark/extdata/single-omic/metagenomic/gevers/otutable.txt'
clinPath = 'benchmark/extdata/single-omic/metagenomic/gevers/clinical.txt'
taxPath = 'benchmark/extdata/single-omic/metagenomic/gevers/TaxonomyTable.rds'

otu = read.delim2(otuPath, header = T, sep = '\t')
clin = read.delim2(clinPath, header = T, sep = '\t')
rownames(clin) = clin$X.SampleID
tax = readRDS(taxPath)

TAXA = tax_table(as.matrix(tax))
otu = subset(otu, select = -c(X.OTU.ID))
rownames(otu) = rownames(tax)
otu = as.matrix(otu)
OTU = otu_table(otu, taxa_are_rows = T)

CLIN = sample_data(clin)

gevers = phyloseq(OTU, CLIN, TAXA)
gevers = tax_glom(gevers, 'Rank6')

gevers.f = filter_taxa(gevers, flist, prune=T)

saveRDS(gevers, file = 'benchmark/data/single-omic/meta_gevers_complete.rds')
saveRDS(gevers.f, file = 'benchmark/data/single-omic/meta_gevers_filtered.rds')


