cesga = F

if (cesga == T) {
  setwd('/mnt/netapp2/Store_uni/home/ulc/co/jlb/git/Review-MLAID')
} else{
  setwd('~/git/Review-MLAID/')
}

outPath = 'benchmark/data/preciseads/'


# Load data
# ===
require(data.table)
# clinical
clinical = read.delim2('benchmark/extdata/PRECISEADS/Metadata.tsv')
rownames(clinical) = make.names(rownames(clinical))
clinical$cases = as.factor(ifelse(clinical$Diagnosis == 'CTRL', 'control', 'case'))
cases = which(clinical$cases == 'case')

# expression
exp = read.table('benchmark/extdata/PRECISEADS/TPMs_filtered.tsv', header = T)

# methylation
met = fread('benchmark/extdata/PRECISEADS/MValues_filtered.tsv', data.table = F)
rownames(met) = met$V1; met = met[,-1]
names(met) = make.names(names(met))

stopifnot(names(met) == names(exp))
stopifnot(rownames(clinical) == names(met))

# Fast correlation based feature selection (FCBF)
# ===
require(FCBF)
# In expression
print('FCBF in expression ...')
dis_exp = as.data.frame(discretize_exprs(exp))
fit.exp = fcbf(feature_table = dis_exp,
            target_vector = clinical$cases,
            minimum_su = 0.0005)

exp.fcbf = exp[fit.exp$index,]
exp.fcbf.noctrls = exp[fit.exp$index, cases]

saveRDS(exp.fcbf, file = paste0(outPath, 'expression_fcbf.rds'))
saveRDS(exp.fcbf.noctrls, file = paste0(outPath, 'expression_fcbf_noctrls.rds'))

# In methylation
print('FCBF  in methylation ...')
dis_met = as.data.frame(discretize_exprs(met))
fit.met = fcbf(feature_table = dis_met,
                target_vector = clinical$cases,
                minimum_su = 0.0005)

met.fcbf = met[fit.met$index,]
met.fcbf.noctrls = met[fit.met$index, cases]

saveRDS(met.fcbf, file = paste0(outPath, 'methylation_fcbf.rds'))
saveRDS(met.fcbf.noctrls, file = paste0(outPath, 'methylation_fcbf_noctrls.rds'))

# Select features by levene test
# ===
# In expression data
print('Levene & kruskal feature selection  in expression ...')
res.levene.exp = apply(exp,  1, function(x) car::leveneTest(x, clinical$cases)$`Pr(>F)`[1])
res.levene.adj.exp = p.adjust(res.levene.exp, method = 'bonferroni')

res.kruskal.exp = apply(exp, 1, function(x) kruskal.test(x, clinical$cases)$p.value)
res.kruskal.adj.exp = p.adjust(res.kruskal.exp, method = 'bonferroni')

l = which(res.levene.adj.exp < 0.05)
k = which(res.kruskal.adj.exp < 0.05)

print(paste0('Selected features: ', length(intersect(l, k))))

exp.filter = exp[intersect(l, k),]
exp.filter.noctrls = exp[intersect(l, k), cases]

saveRDS(exp.filter, file = paste0(outPath, 'expression_filter_levene_kruskal.rds'))
saveRDS(exp.filter.noctrls, file = paste0(outPath, 'expression_filter_levene_kruskal_noctrls.rds'))


# In mehtylation data
print('Levene & kruskal feature selection  in methylation ...')
res.levene.met = apply(met, 1, function(x) car::leveneTest(x, clinical$cases)$`Pr(>F)`[1])
res.levene.adj.met = p.adjust(res.levene.met, method = 'bonferroni')

res.kruskal.met = apply(met, 1, function(x) kruskal.test(x, clinical$cases)$p.value)
res.kruskal.adj.met = p.adjust(res.kruskal.met, method = 'bonferroni')

l = which(res.levene.adj.met < 0.05)
k = which(res.kruskal.adj.met < 0.05)

print(paste0('Selected features: ', length(intersect(l, k))))

met.f = met[which(res.levene.adj.met < 0.05),]
met.f.noctrls = met[which(res.levene.adj.met < 0.05), cases]

saveRDS(met.f, file = paste0(outPath, 'methylation_filter_levene_kruskal.rds'))
saveRDS(met.f.noctrls, file = paste0(outPath, 'methylation_filter_levene_kruskal_noctrls.rds'))

# saveRDS(exp.f, file = '~/git/Review-MLAID/benchmark/data/preciseads/expression_filtered_levene.rds')
# saveRDS(met.f, file = '~/git/Review-MLAID/benchmark/data/preciseads/methylation_filtered_levene.rds')
