#!/bin/bash

#SBATCH -p amd-shared
#SBATCH --qos amd-shared
#SBATCH -t 10:00:00
#SBATCH --mem=200GB
#SBATCH -n 15
#SBATCH --mail-user=joselb1993@gmail.com
#SBATCH --mail-type=BEGIN,END

module load cesga/2018 gcc/6.4.0 R/4.0.2
Rscript /mnt/netapp2/Store_uni/home/ulc/co/jlb/git/Review-MLAID/benchmark/code/multi-omic/run_multiomic_GSE117931.r