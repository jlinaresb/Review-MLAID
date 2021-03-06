setwd('~/git/Review-MLAID/benchmark/results/single-omic/')

files = list.files()
time = list()
l = list()
res = list()
# f = 1
for (f in seq_along(files)) {
  data = readRDS(files[f])
  # data = data[-c(3)]
  
  for (i in seq_along(data)) {
    l[[i]] = data.frame(
      ids = names(data[[i]]$labels),
      labels = data[[i]]$labels
    )
    colnames(l[[i]])[2] = names(data)[i]
  }
  
  x = as.data.frame(l)
  
  res[[f]] = x[,-grep('id', colnames(x))]
  
} 

names(res) = files

# ARI calculation
# ===

require(mclust)
require(corrplot)
require(viridis)
# kk = res$metagenomic_morgan.rds
corrplots = list()
for (i in seq_along(res)) {
  kk = res[[i]]
  name = names(res)[i]
  
  n <- ncol(kk)
  algs = colnames(kk)
  
  ari <- function(i, j, data) {adjustedRandIndex(data[,i], data[,j])}
  corp <- Vectorize(ari, vectorize.args = list("i","j"))
  ari_res = outer(1:n, 1:n, corp, data=kk)
  ari_res[ari_res < 0] = 0
  
  rownames(ari_res) = algs
  colnames(ari_res) = algs
  
  corrplots[[i]] = corrplot(ari_res,
                            type = 'lower',
                            col = viridis(10),
                            tl.col = 'black',
                            title = name)
}
names(corrplots) = names(res)


# 
require(ggcorrplot)
require(viridis)
plots = list()
for (i in seq_along(corrplots)) {
  toPlot = corrplots[[i]]$corr
  
  toPlot[which(is.nan(toPlot))] = 0
  
  plots[[i]] = ggcorrplot(toPlot,
                          hc.order = F,
                          type = 'lower',
                          lab = T,
                          insig = 'blank',
                          title = strsplit(names(corrplots)[i], '_')[[1]][3],
                          show.legend = F,
                          show.diag = F, 
                          tl.cex = 7,
                          lab_size = 2, 
                          digits = 1,
                          colors = c(viridis(1),viridis(2),viridis(3)))
  
}
require(ggpubr)
finalPlot = ggarrange(plots[[6]], plots[[7]], plots[[8]],
          plots[[1]], plots[[2]], plots[[3]], 
          plots[[4]], plots[[5]],
          ncol = 3, nrow = 4, 
          common.legend = T)

ggsave(plot = finalPlot, 
       filename = 'ari_singleomic.pdf',
       device = 'pdf',
       path = '~/git/Review-MLAID/benchmark/plots/',
       height = 9,
       width = 7)
