# config file
# ===

base.dir = '~/git/Review-MLAID/'

cluster.type = 'single'
input.dir.path = paste0(base.dir, 'benchmark/data/preciseads/')
path.algs = paste0(base.dir, 'benchmark/code/preciseads/run/', cluster.type, '/')
pattern.algs = 'run_'


part = 'cola-corta'
qos = 'default'
time = '10:00:00'
mem = '32GB'
ntasks = '20'
