# Scopus search
require(data.table)
require(rscopus)
require(dplyr)

# Add API KEY
# ===
rscopus::set_api_key('52a7417360eefd56eda60627cc2e1620')

queryID = 'aid_clusters'

# Build the query
# ===
# Disease (systemic lupus erythematosus)
diseases = c('systemic lupus erythematosus', 'rheumatoid arthritis', 'sjögren syndrome', 'type 1 diabetes', 'systemic sclerosis', 'multiple sclerosis', 'inflammatory bowel disease')

# Omic data type
# Genome
genome = c('genome', 'genomic', 'SNP', 'GWAS', 'gene sequence', 'CNV', 'CNA', 'copy number variation', 'copy number aberration')
# Transcriptomic
transcriptome = c('transcriptome', 'transcriptomic', 'gene expression', 'miRNA', 'mRNA', 'RNASeq', 'RNA-Seq', 
                  'microarray', 'lncRNA', 'transcription factor', 'co-expression', 'exome')
# Epigenomic
epigenome = c('epigenome', 'epigenomic', 'methylation', 'acethylation', 'histone', 'chromatin')
# Proteomic
proteome = c('proteome', 'proteomic', 'protein expression', 'protein array', 'cytokine array', 'phosphoproteomic')
# Metabolomic
metabolome = c('metabolome', 'metabolomic', 'metabolites')
# Metagenomic
metagenome = c('metagenome', 'metagenomic', 'microbiota', '16S', 'microbiome')

# Approach
approach = c('unsupervised', 'cluster', 'machine learning')

omics = list(epigenome = epigenome,
             transcriptome = transcriptome,
             proteome = proteome,
             metabolome = metabolome,
             metagenome = metagenome,
             genome = genome)

queries = list()
query = list()
for (d in seq_along(diseases)) {
  for (i in seq_along(omics)) {
    y = paste0('TITLE-ABS-KEY(', diseases[d], ') & ', 'TITLE-ABS-KEY(', approach, ') & ')
    x = paste0('TITLE-ABS-KEY(', omics[[i]], ')')
    x = paste(x, collapse = " OR ")
    queries[[i]] = paste0(y, x)
  }
  query[[d]] = unlist(queries, use.names = F)
}
query = unlist(query)

# Omic queries
# ===
# Gene expression
queryRes = list()
completeArticle = list()
for (q in seq_along(query)) {
  
  # Run query 
  print(paste0('Search for ...', query[q]))
  completeArticle[[q]] <- scopus_search(
    query = query[q], 
    view = "COMPLETE", 
    count = 25)
  
  # Identify disease from query
  disease = strsplit(query[q], '&')[[1]][1]
  disease = gsub("[\\(\\)]", "", regmatches(disease, gregexpr("\\(.*?\\)", disease))[[1]])
  
  # Format results
  res = list()
  for (i in seq_along(completeArticle[[q]]$entries)) {
    
    pubmedID = completeArticle[[q]]$entries[[i]]$`pubmed-id`
    if (is.null(pubmedID)) {
      pubmedID = NA
    }
    
    DOI = completeArticle[[q]]$entries[[i]]$`prism:doi`
    if (is.null(DOI)){
      DOI = NA
    }
    
    type = completeArticle[[q]]$entries[[i]]$`prism:aggregationType`
    if (is.null(type)){
      type = NA
    }
    
    abstract = completeArticle[[q]]$entries[[i]]$`dc:description`
    if (is.null(abstract)){
      abstract = NA
    }
    
    omicType = strsplit(strsplit(query[q], ' & ')[[1]][3], ' OR ')[[1]][1]
    omicType = gsub("[\\(\\)]", "", regmatches(omicType, gregexpr("\\(.*?\\)", omicType))[[1]])
    
    nAuthors = as.numeric(completeArticle[[q]]$entries[[i]]$`author-count`$`@total`)
    if (length(nAuthors) == 0) {
      nAuthors = NA
      lastAuthor = NA
      firstAuthor = NA
    } else if (nAuthors == 0) {
      lastAuthor = NA
      firstAuthor = NA
    } else if (nAuthors > 0 & nAuthors <= 100){
      lastAuthor = completeArticle[[q]]$entries[[i]]$author[[nAuthors]]$authname
      firstAuthor = completeArticle[[q]]$entries[[i]]$author[[1]]$authname
    } else if (nAuthors > 100){
      lastAuthor = 'More than 100 authors'
      firstAuthor = completeArticle[[q]]$entries[[i]]$author[[1]]$authname
    }
    
    if (completeArticle[[q]]$total_results == 0) {
      res[[i]] = NULL
    } else{
      # Create data.frame with paper information
      res[[i]] = data.frame(
        disease = disease,
        pubmedID = pubmedID,
        URL = completeArticle[[q]]$entries[[i]]$`prism:url`,
        title = completeArticle[[q]]$entries[[i]]$`dc:title`,
        DOI = DOI,
        journal = completeArticle[[q]]$entries[[i]]$`prism:publicationName`,
        date = completeArticle[[q]]$entries[[i]]$`prism:coverDisplayDate`,
        type = type,
        nAuthors = nAuthors,
        firstAuthor = firstAuthor,
        lastAuthor = lastAuthor,
        abstract = abstract,
        omicType = omicType,
        citations = completeArticle[[q]]$entries[[i]]$`citedby-count`,
        OA = completeArticle[[q]]$entries[[i]]$openaccess,
        subtype = ifelse(is.null(completeArticle[[q]]$entries[[i]]$subtypeDescription),
                         NA, completeArticle[[q]]$entries[[i]]$subtypeDescription)
      )
    }
  }
  result = rbindlist(res)
  
  # save query results
  queryRes[[q]] = result
  
}
scopus = rbindlist(queryRes)

# Save it!
# ===
saveRDS(scopus, file = paste0('~/git/Review-MLAID/data/scopus_', queryID, '.rds'))
