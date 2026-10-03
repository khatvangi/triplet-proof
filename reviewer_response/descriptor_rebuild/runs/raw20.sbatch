#!/bin/bash
#SBATCH -J s3raw20
#SBATCH -c 32
#SBATCH -t 02:00:00
#SBATCH -o runs/%x_%A_%a.out
#SBATCH --array=0-1
# table S3 column: 20 standardized descriptors (pre-PCA), rebuilt and earlier matrices.
# same fixed seeds as runs/new.json etc., so all columns see identical nulls.
cd /storage/kiran-stuff/triplet-proof/reviewer_response/descriptor_rebuild
TAGS=(new_raw20 old_raw20); CSVS=(RAW20_new.csv RAW20_old.csv)
python stats.py ${TAGS[$SLURM_ARRAY_TASK_ID]} ${CSVS[$SLURM_ARRAY_TASK_ID]} d_ 1
