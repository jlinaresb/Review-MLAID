# Download from PubMed
# ====

require(easyPubMed)
require(data.table)

setwd('/home/joselinares/git/Review-MLAID/data/')

query = 'machine learning[Title/Abstract] AND type 1 diabetes'
out = batch_pubmed_download(pubmed_query_string = query,
                            dest_file_prefix = 'example_',
                            encoding = 'ASCII')

df = list()
for (i in seq_along(out)) {
  file = out[[i]]
  df[[i]] = table_articles_byAuth(pubmed_data = file,
                             included_authors = 'first',
                             max_chars = 0,
                             encoding = "ASCII")
}
res = rbindlist(df)
res = subset(res, select = -c(abstract, month, day, journal, firstname, address, email))

saveRDS(res, file = '~/git/Review-MLAID/data/pubmed.rds')