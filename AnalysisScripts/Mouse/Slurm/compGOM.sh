#!/bin/bash
#SBATCH --job-name=RCTD    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=lucast3@gene.com     # Where to send mail	
#SBATCH --ntasks=1                # Run on a single CPU
#SBATCH --partition="himem"         #Run on Himem or defq partition
#SBATCH --mem=450gb                     # Job memory request
#SBATCH --cpus-per-task=7            # Number of CPU cores per task
#SBATCH --output=/gstore/project/neurodegen_meta/CompGOM.json   # Standard output and error log

pwd; hostname; date

ml R/prd

Rscript /gstore/project/neurodegen_meta/neurodegeneration_meta-analysis/Rscripts/TopicModeling/writeTablesExtractBIC.R
