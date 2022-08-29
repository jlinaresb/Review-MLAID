#!/bin/bash 
#SBATCH -p cola-corta 
#SBATCH --qos default 
#SBATCH -t 10:00:00 
#SBATCH --mem 32GB 
#SBATCH -n 20 
module load cesga/2018 gcc/6.4.0 R/4.0.2 
expPath=~/git/Review-MLAID/benchmark/data/preciseads/exp/expression_fcbf.rds
metPath=~/git/Review-MLAID/benchmark/data/preciseads/met/methylation_fcbf.rds
Rscript ~/git/Review-MLAID/benchmark/code/preciseads/run/multi/run_mclust.r $expPath $metPath