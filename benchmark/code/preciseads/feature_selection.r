setwd('~/git/Review-MLAID/benchmark/extdata/PRECISEADS/')

# Load data
# ===
require(data.table)
exp = read.table('TPMs_filtered.tsv', header = T)

clinical = read.delim2('Metadata.tsv')
rownames(clinical) = make.names(rownames(clinical))

met = fread('MValues_filtered.tsv', data.table = F)
rownames(met) = met$V1; met = met[,-1]
names(met) = make.names(names(met))

stopifnot(names(met) == names(exp))
stopifnot(rownames(clinical) == names(met))

# Select features by levene test
# ===
clinical$cases = as.factor(ifelse(clinical$Diagnosis == 'CTRL', 'control', 'case'))
# In expression data
res.levene.exp = apply(exp,  1, function(x) car::leveneTest(x, clinical$cases)$`Pr(>F)`[1])
res.levene.adj.exp = p.adjust(res.levene.exp, method = 'bonferroni')

head(res.levene.adj.exp)
length(which(res.levene.adj.exp < 0.05))

# In mehtylation data
res.levene.met = apply(met, 1, function(x) car::leveneTest(x, clinical$cases)$`Pr(>F)`[1])
res.levene.adj.met = p.adjust(res.levene.met, method = 'bonferroni')

head(res.levene.adj.met)
length(which(res.levene.adj.met < 0.05))

