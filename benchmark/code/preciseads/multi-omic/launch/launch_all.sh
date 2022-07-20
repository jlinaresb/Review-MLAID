#!/bin/bash

sbatch /home/joselinares/git/Review-MLAID/benchmark/code/preciseads/multi-omic/launch/lch_coca.sh
sbatch /home/joselinares/git/Review-MLAID/benchmark/code/preciseads/multi-omic/launch/lch_kmeans.sh
sbatch /home/joselinares/git/Review-MLAID/benchmark/code/preciseads/multi-omic/launch/lch_snf.sh
sbatch /home/joselinares/git/Review-MLAID/benchmark/code/preciseads/multi-omic/launch/lch_hierarchical.sh
sbatch /home/joselinares/git/Review-MLAID/benchmark/code/preciseads/multi-omic/launch/lch_mclust.sh
sbatch /home/joselinares/git/Review-MLAID/benchmark/code/preciseads/multi-omic/launch/lch_som.sh
sbatch /home/joselinares/git/Review-MLAID/benchmark/code/preciseads/multi-omic/launch/lch_icluster.sh
sbatch /home/joselinares/git/Review-MLAID/benchmark/code/preciseads/multi-omic/launch/lch_nmf.sh