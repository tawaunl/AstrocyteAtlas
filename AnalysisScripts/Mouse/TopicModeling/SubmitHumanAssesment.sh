#!/bin/bash

#SBATCH --job-name=assess_models_human
#SBATCH --output=/gpfs/scratchfs01/site/u/lucast3/AstrocyteAtlas/AnalysisScripts/Mouse/TopicModeling/logs/assess_models_human_%j.log
#SBATCH --error=/gpfs/scratchfs01/site/u/lucast3/AstrocyteAtlas/AnalysisScripts/Mouse/TopicModeling/logs/assess_models_human_%j.err
#SBATCH --qos=1d
#SBATCH --partition=batch_cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=200G
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=lucast3@external.gene.com

# Load necessary modules
# Ensure the R environment is loaded. This might vary depending on your HPC setup.
ml CEDAR
ml R/cedar_r4.5_bioc3.21-release

# Run the script
echo "Starting Human assessment"

Rscript /gpfs/scratchfs01/site/u/lucast3/AstrocyteAtlas/AnalysisScripts/Mouse/TopicModeling/assess_models_human.R

echo "Human assessment finished."
