#!/bin/bash

#SBATCH --partition=himem
#SBATCH -c 14
#SBATCH --job-name='get per cluster files'
#SBATCH --output='/gstore/scratch/u/novikovg/Astrocytes_meta/AD_MS_PD_cellbender_harmony_2024/harmony_integration/another_run_after_removing_high_mito_clusters/clustering/per_cluster_file_%j.out'
#SBATCH -e /gstore/scratch/u/novikovg/Astrocytes_meta/AD_MS_PD_cellbender_harmony_2024/harmony_integration/another_run_after_removing_high_mito_clusters/clustering/per_cluster_file_%j.err    # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=novikovg@gene.com

module load R
Rscript $1
