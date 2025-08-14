#!/bin/bash
#SBATCH -n 8           # 1 core
#SBATCH -t 0-01:00:00   # 30 min
#SBATCH -J deduplicate # sensible name for the job
#SBATCH --output=opt.log   # Standard output and error log.
#SBATCH -p mit_normal # partition
#SBATCH --mem-per-cpu=4G # 8 gb

conda activate mechdata

python deduplicate.py --data_dir="results/" --output_file="results/templated_dedupe.txt"