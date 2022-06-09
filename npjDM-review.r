# Load data (pdf --> csv)
# ===
xx = read.csv('~/git/Review-MLAID/data/supp-npjDM.csv', header = T, skip = 12)

# Preprocess
# ===
y = xx[1:401,]
y = y[,1:7]

y2 = y[which(y$Paper != '' &  y$Paper != 'Paper'), ]
y3 = y2[match(setdiff(y2$Paper, c(1:16)), y2$Paper), ]
y4 = y3[which(y3$MultipleAIDsStudied != ''), ]
rownames(y4) = 1:nrow(y4)

d = y3[which(y3$MultipleAIDsStudied == ''), ]
d = d[which(d$Machine.LearningMethod == '' & d$Study.Size..N. == '' & d$Prediction.orClassificationTask == ''), ]$Paper
d = d[-c(2:4)]

kk = c(rep(d[1], 41), rep(d[2], 32), rep(d[3], 30), rep(d[4], 16),
       rep(d[5], 14), rep(d[6], 11), rep(d[7], 7), rep(d[8], 6),
       rep(d[9], 5), rep(d[10], 4), rep(d[11], 1), rep(d[12], 1))

y4$disease = kk

data = y4

# Descriptive analysis
# ===
table(data$disease)
table(data$Type.of.Data)

# Is there omic data??
# ===
omic = unique(c(grep('Gene', data$Type.of.Data),
grep('Trans', data$Type.of.Data),
grep('omic', data$Type.of.Data),
grep('Meta', data$Type.of.Data),
grep('Micro', data$Type.of.Data),
grep('GWAS', data$Type.of.Data),
grep('Expres', data$Type.of.Data),
grep('Exome', data$Type.of.Data),
grep('Proteome', data$Type.of.Data),
grep('Mass', data$Type.of.Data),
grep('SNP', data$Type.of.Data),
grep('RNA', data$Type.of.Data)))

data.omic = data[match(omic, rownames(data)),]

View(data.omic[,c(1,7,8)])
