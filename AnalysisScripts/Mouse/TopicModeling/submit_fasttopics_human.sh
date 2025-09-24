#!/bin/bash
#SBATCH --job-name=fasttopics_array_human      # Job name
#SBATCH --array=4,6,8,10,12,14                     # The range of k values to run. Adjust as needed.
#SBATCH --output=/gpfs/scratchfs01/site/u/lucast3/AstrocyteAtlas/AnalysisScripts/Mouse/TopicModeling/logs/fasttopics_%A_%a.out  # Standard output file
#SBATCH --error=/gpfs/scratchfs01/site/u/lucast3/AstrocyteAtlas/AnalysisScripts/Mouse/TopicModeling/logs/fasttopics_%A_%a.err   # Standard error file
#SBATCH --mem-per-cpu=250G                # Memory per CPU. Adjust based on your data size.             
#SBATCH --cpus-per-task=2                # Number of CPUs per task. 
#SBATCH --partition=batch_cpu
#SBATCH --qos=1d

# --- Prepare Environment ---
# Ensure the R environment is loaded. This might vary depending on your HPC setup.

ml CEDAR
ml R/cedar_r4.5_bioc3.21-release

# --- Get the number of topics (k) from the SLURM_ARRAY_TASK_ID ---
# The %a in the array job name maps to the SLURM_ARRAY_TASK_ID.
K=${SLURM_ARRAY_TASK_ID}

# --- Run the R script with the current k value ---
# The Rscript command runs the R script from the command line.
echo "Starting job for k = $K"
Rscript /gpfs/scratchfs01/site/u/lucast3/AstrocyteAtlas/AnalysisScripts/Mouse/TopicModeling/run_fasttopics_human.R $K

echo "Job for k = $K finished."