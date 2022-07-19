# Filter articles
# ===
setwd('~/git/Review-MLAID/data/')
source('~/git/Review-MLAID/utils.r')

# Load search
# ===
scopus = readRDS('scopus_aid_clusters.rds')

# Filter by: article type
# ===
scopus = scopus[which(scopus$subtype == 'Article'), ]

# Delete articles without DOI
# ===
scopus = scopus[which(!is.na(scopus$DOI)),]
rownames(scopus) = 1:nrow(scopus)

# Label duplicated papers
# ===
scopus$duplicated = NA
scopus$duplicated[which(duplicated(scopus$DOI) == TRUE)] = 'yes'
table(scopus$duplicated)

# Remove duplicated papers
# ===
scopus.d = scopus[which(is.na(scopus$duplicated)),]

# Check if article is based on clustering
# ===
# scopus.c = grepCluster(scopus.d)
query1 = scopus.d[grep('disease classification', tolower(scopus.d$abstract)),]$DOI
query2 = scopus.d[grep('patient stratification', tolower(scopus.d$abstract)),]$DOI
query3 = scopus.d[grep('cluster patients', tolower(scopus.d$abstract)),]$DOI
query4 = scopus.d[grep('disease stratification', tolower(scopus.d$abstract)),]$DOI
byCohort = scopus.d[grep('cohort', tolower(scopus.d$abstract)),]$DOI

# news = setdiff(c(query1, query2, query3, query4), byCohort)
# news = news[c(1:6, 8:9)]
# 
# 
# View(scopus.d[match(news, scopus.d$DOI),])

# Export to excel file
# ===
require(xlsx)
# write.xlsx(scopus.c, file = "~/git/Review-MLAID/review-papers/SLE/aid_cluster.xlsx")
write.xlsx(byCohort, file = "~/git/Review-MLAID/review-papers/aid_cohort.xlsx")

