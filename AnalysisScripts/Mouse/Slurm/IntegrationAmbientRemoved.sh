#!/bin/bash
#SBATCH --job-name=AstrocyteIntegration_AmbientRemoved    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=lucast3@gene.com     # Where to send mail	
#SBATCH --ntasks=1                # Run on a single CPU
#SBATCH --partition="himem"         #Run on Himem or defq partition
#SBATCH --mem=450gb                     # Job memory request
#SBATCH --output=/gstore/project/neurodegen_meta/AstrocyteIntegration_outAmbientRemoval_new.json   # Standard output and error log
#SBATCH --cpus-per-task=10            # Number of CPU cores per task

pwd; hostname; date

ml R/prd

Rscript /gstore/data/project/neurodegen_meta/neurodegeneration_meta-analysis/Rscripts/CellBenderIntegration.R
 
