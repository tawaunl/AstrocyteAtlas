#!/bin/bash
#SBATCH --job-name=CellBenderGetAstrocytes    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=lucast3@gene.com     # Where to send mail	
#SBATCH --ntasks=1                # Run on a single CPU
#SBATCH --partition="defq"         #Run on Himem or defq partition
#SBATCH --qos=long
#SBATCH --mem=150gb                     # Job memory request
#SBATCH --output=/gstore/data/project/neurodegen_meta/CEllBenderGetAstrocytes_out.json   # Standard output and error log
#SBATCH --cpus-per-task=4            # Number of CPU cores per task

pwd; hostname; date

ml R/prd

Rscript /gstore/data/project/neurodegen_meta/neurodegeneration_meta-analysis/Rscripts/LabelCellBenderCounts.R
 
