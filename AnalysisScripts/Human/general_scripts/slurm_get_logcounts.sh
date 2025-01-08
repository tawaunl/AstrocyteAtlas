#!/bin/bash

#SBATCH --partition=gpu
#SBATCH --gres=gpu:2
#SBATCH -c 28
#SBATCH --mem=400G
#SBATCH --job-name='Get var genes and run harmony with chan'
#SBATCH --output='/gstore/scratch/u/novikovg/Astrocytes_meta/AD_MS_PD_cellbender_harmony_2024/harmony_integration/get_logcounts_%j.out'
#SBATCH -e /gstore/scratch/u/novikovg/Astrocytes_meta/AD_MS_PD_cellbender_harmony_2024/harmony_integration/get_logcounts_%j.err    # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=novikovg@gene.com

module load R
Rscript /gstore/scratch/u/novikovg/Astrocytes_meta/AD_MS_PD_cellbender_harmony_2024/harmony_integration/get_logcounts_matrix.R
