#!/bin/bash

#SBATCH --partition=gpu
#SBATCH --gres=gpu:2
#SBATCH -c 28
#SBATCH --mem=400G
#SBATCH --job-name='Combine datasets v3'
#SBATCH --output='/gstore/scratch/u/novikovg/Astrocytes_meta/AD_MS_PD_cellbender_harmony_2024/harmony_integration/Combine_datasets_v3_%j.out'
#SBATCH -e /gstore/scratch/u/novikovg/Astrocytes_meta/AD_MS_PD_cellbender_harmony_2024/harmony_integration/Combine_datasets_v3_%j.err    # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=novikovg@gene.com

module load R
Rscript /gstore/scratch/u/novikovg/Astrocytes_meta/AD_MS_PD_cellbender_harmony_2024/harmony_integration/combine_all_datasets.v3.liddelow_cellranger_the_rest_cellbender.R
