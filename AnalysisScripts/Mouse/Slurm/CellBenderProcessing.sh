#!/bin/bash
#SBATCH --job-name=CellBenderIntegratedProcessing   # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=lucast3@gene.com     # Where to send mail	
#SBATCH --ntasks=1                # Run on a single CPU
#SBATCH --partition="defq"         #Run on Himem or defq partition
#SBATCH --mem=220gb                     # Job memory request
#SBATCH --output=/gstore/project/neurodegen_meta/CellBenderIntegratedProcessing.json   # Standard output and error log
#SBATCH --cpus-per-task=7            # Number of CPU cores per task

pwd; hostname; date

ml R/prd

Rscript -e "rmarkdown::render('/gstore/project/neurodegen_meta/neurodegeneration_meta-analysis/Rscripts/FinalScripts/CellBenderIntegratedProcessing.Rmd')"
