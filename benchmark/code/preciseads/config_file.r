# config file
# ===

#base.dir = '~/git/Review-MLAID/'
base.dir = '/mnt/netapp2/Store_uni/home/ulc/co/jlb/git/Review-MLAID/'

cluster.type = 'multi'  #multi single
input.dir.path = paste0(base.dir, 'benchmark/data/preciseads/')
path.algs = paste0(base.dir, 'benchmark/code/preciseads/run/', cluster.type, '/')
pattern.algs = 'run_'


part = 'thin-shared'
qos = 'default'
time = '4-04:00:00'
mem = '100GB'
ntasks = '20'
