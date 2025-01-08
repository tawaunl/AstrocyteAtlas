#!/bin/bash
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=lucast3@gene.com     # Where to send mail	
#SBATCH --ntasks=1                # Run on a single CPU
#SBATCH --mem=64gb                     # Job memory request
#SBATCH --cpus-per-task=2            # Number of CPU cores per task
#SBATCH -p gpu
#SBATCH --gres=gpu:2        # request 2 gpus

pwd; hostname; date

export PATH=~/cellranger-7.0.1:$PATH

ml Anaconda3
ml CUDA
source activate CellBender

source activate CellBender
source activate CellBender
conda info --envs

cellbender remove-background \
                 --cuda \
                 --input "$2" \
                 --output $3 \
                 --expected-cells $4 \
                 --total-droplets-included $5 \
                 --low-count-threshold 10 \
                 --fpr 0.01 0.05 0.1 0.15 \
                 --epochs 150
