require(phyloseq)

setwd('~/git/Review-MLAID/benchmark/')

# Morgan
# ===
otuPath = 'extdata/single-omic/metagenomic/morgan/otutable.txt'
clinPath = 'extdata/single-omic/metagenomic/morgan/clinical.txt'
taxPath = 'extdata/single-omic/metagenomic/morgan/TaxonomyTable.rds'

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

# Gevers
# ===
otuPath = 'extdata/single-omic/metagenomic/gevers/otutable.txt'
clinPath = 'extdata/single-omic/metagenomic/gevers/clinical.txt'
taxPath = 'extdata/single-omic/metagenomic/gevers/TaxonomyTable.rds'

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

rm(list = setdiff(ls(), c('morgan', 'gevers')))


# Convert to data.frame
# ===

morgan_ = cbind.data.frame(
  t(morgan@otu_table@.Data),
  Condition = get_variable(morgan, 'Var')
)

gevers_ = cbind.data.frame(
  t(gevers@otu_table@.Data),
  Condition = get_variable(gevers, 'Var')
)

saveRDS(morgan_, file = 'extdata/single-omic/metagenomic/filtered/morgan.rds')
saveRDS(gevers_, file = 'extdata/single-omic/metagenomic/filtered/gevers.rds')
