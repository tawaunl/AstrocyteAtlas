#!/bin/bash
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=lucast3@gene.com     # Where to send mail	
#SBATCH --ntasks=1                # Run on a single CPU
#SBATCH --partition="himem"         #Run on Himem or defq partition
#SBATCH --mem=450gb                     # Job memory request
#SBATCH --cpus-per-task=3            # Number of CPU cores per task
#SBATCH --qos=long


pwd; hostname; date

ml R/prd

Rscript /gstore/project/neurodegen_meta/neurodegeneration_meta-analysis/Rscripts/TopicModeling/CountClust.R $1
