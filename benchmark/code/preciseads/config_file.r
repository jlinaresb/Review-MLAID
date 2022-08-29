# config file
# ===

#base.dir = '~/git/Review-MLAID/'
base.dir = '/mnt/netapp2/Store_uni/home/ulc/co/jlb/git/Review-MLAID/'

cluster.type = 'multi'
input.dir.path = paste0(base.dir, 'benchmark/data/preciseads/')
path.algs = paste0(base.dir, 'benchmark/code/preciseads/run/', cluster.type, '/')
pattern.algs = 'run_'


part = 'cola-corta'
qos = 'default'
time = '10:00:00'
mem = '64GB'
ntasks = '15'
