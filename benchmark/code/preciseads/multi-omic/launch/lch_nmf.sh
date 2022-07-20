#!/bin/bash

#SBATCH -p cola-corta
#SBATCH --qos default
#SBATCH -t 10:00:00
#SBATCH --mem=64GB
#SBATCH -n 15
#SBATCH --mail-user=joselb1993@gmail.com
#SBATCH --mail-type=BEGIN,END

module load cesga/2018 gcc/6.4.0 R/4.0.2
Rscript /mnt/netapp2/Store_uni/home/ulc/co/jlb/git/Review-MLAID/benchmark/code/preciseads/multi-omic/run_nmf.r