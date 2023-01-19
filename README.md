# Molecular-based clustering models in Autoimmune diseases. A systematic review and benchmarking


## Study Design
This is a systematic review of unsupervised models carried out in Autoimmune diseases (AID), which follows Preferred Reporting Items for Systematic Reviews [PRISMA-P](http://prisma-statement.org/Extensions/Protocols) guidelines.


## Motivation
AID are characterised by a clinical heterogenity, impeding a early diagnosis and a effective treatment. Previous studies [Barturen, G. et al. 2020](https://onlinelibrary.wiley.com/doi/10.1002/art.41610) have identified common clusters in seven AIDs. Based on this previous results, we aim to identify records that have used omic data to clusterized patients of several AIDs.


## Search Strategies
The literature search was carried out using Scopus API's through [rscopus](https://johnmuschelli.com/rscopus/) package. The search was updated from February 10th 2022 to November 3rd 2022.


## Diseases included:
This work is running in the following diseases:

* SLE = Systematic Lupus Erythematosus
* RA = Rheumatoid Arthritis
* SjS = Sjögren's Syndrome
* SSc = Systemic Sclerosis
* MSc = Multiple Sclerosis
* T1D = Type 1 Diabetes
* IBD = Inflammatory Bowel Disease

All except IBD and MSc are included in [ADEx](https://adex.genyo.es/).


### Query on Scopus
Query was constructed by three conditional aspects: disease, approach and omic type. Query was built on [search_scopus_v2.r](https://github.com/jlinaresb/Review-MLAID/blob/master/search_scopus_v2.r) script. 


#### Disease type:
``` r
'TITLE-ABS-KEY("systemic lupus erythematosus") OR TITLE-ABS-KEY("rheumatoid arthritis") OR TITLE-ABS-KEY("sjögren syndrome") OR TITLE-ABS-KEY("type 1 diabetes") OR TITLE-ABS-KEY("systemic sclerosis") OR TITLE-ABS-KEY("multiple sclerosis") OR TITLE-ABS-KEY("inflammatory bowel disease") '
```

#### Approach
``` r
'TITLE-ABS-KEY("machine learning") OR TITLE-ABS-KEY("cluster") OR TITLE-ABS-KEY("unsupervised')'
```

#### Omic type

1. Genomic
``` r
'TITLE-ABS-KEY("genome") OR TITLE-ABS-KEY("genomic") OR TITLE-ABS-KEY("SNP") OR TITLE-ABS-KEY("SNP") OR TITLE-ABS-KEY("GWAS") OR TITLE-ABS-KEY("gene secuence") OR TITLE-ABS-KEY("CNV") OR TITLE-ABS-KEY("CNA") OR TITLE-ABS-KEY("copy number variation") OR TITLE-ABS-KEY("copy number aberration")'
```

2. Transcriptomic
``` r
'TITLE-ABS-KEY("transcriptome") OR TITLE-ABS-KEY("transcriptomic") OR TITLE-ABS-KEY("gene expression") OR TITLE-ABS-KEY("miRNA") OR TITLE-ABS-KEY("mRNA") OR TITLE-ABS-KEY("RNASeq") OR TITLE-ABS-KEY("RNA-Seq") OR TITLE-ABS-KEY("microarray") OR TITLE-ABS-KEY("lncRNA") OR TITLE-ABS-KEY("transcription factor") OR TITLE-ABS-KEY("co-expression") OR TITLE-ABS-KEY("exome")'
```

3. Epigenomic
``` r
'TITLE-ABS-KEY("epigenome") OR TITLE-ABS-KEY("epigenomic") OR TITLE-ABS-KEY("methylation") OR TITLE-ABS-KEY("acethylation") OR TITLE-ABS-KEY("histone") OR TITLE-ABS-KEY("chromatin")'
```

4. Metabolomic
``` r
'TITLE-ABS-KEY("metabolome") OR TITLE-ABS-KEY("metabolomic") OR TITLE-ABS-KEY("metabolites")'
```

5. Proteomic
``` r
'TITLE-ABS-KEY("proteome") OR TITLE-ABS-KEY("proteomic") OR TITLE-ABS-KEY("protein expression") OR TITLE-ABS-KEY("protein array") OR TITLE-ABS-KEY("cytokine array") OR TITLE-ABS-KEY("phosphoproteomic")'
```

6. Metagenomic
``` r
'TITLE-ABS-KEY("metagenome") OR TITLE-ABS-KEY("metagenomic") OR TITLE-ABS-KEY("microbiota") OR TITLE-ABS-KEY("16S") OR TITLE-ABS-KEY("microbiome")'
```


## Quality Assessment of Selected Studies
First, based on Scopus annotation we select only articles types. Therefore, book chapters, conference papers, editorials, letters, notes, reviews and short surverys were removed. Those articles without avialble DOI and duplicated were also removed.

After automatic screening titles, we removed records without the word 'cluster', 'unsupervised', 'stratification' o 'subtype' in the title.
Methodological

##  flowchart and number of papers reviewd at each stage

![flowchart](systematic-review/figures/figure1.jpg?raw=true "flow chart")


## Inclusion criteria
Articles that did not meet the eligibility criteria were discarded. The criteria were composed by three main questions:
  * Patient cohort study
  * Has the study used any type of omic data?
  * Has the study used any unsupervised method to stratify patients?
  * Has the study used an accurate methodology and is it described with enough details to validate its findings?
If all these questions were answered with ‘yes’, the record was selected for synthesis.


## Exclusion Criteria
Articles that accomplish the following characteristics were excluded from the systematic review:

1. Original reseach (Article type)
2. DOI available
3. Models build with any type of omic data
4. Data availability?? 


## Articles classification 
Those articles that meet all eligibility criteria were further classify according this items:

  * Disease(s) includes
  * Tissue used
  * Number of patients
  * Validation cohort?
  * Data used to make the clusters
  * Time points
  * Cluster(s) algorithm(s)
  * Number of clusters identified
  * Is there available data to reproduce this study?
