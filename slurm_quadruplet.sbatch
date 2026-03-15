#!/bin/bash
#SBATCH --job-name=quad_phase2
#SBATCH --output=logs/quad_phase2_%j.out
#SBATCH --error=logs/quad_phase2_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --partition=compute

source ~/.bashrc
conda activate triplet

cd /storage/kiran-stuff/triplet-proof

echo "=== Phase 2 Quadruplet Run ==="
date
python run_phase2_n2n4.py \
    --n-bases 4 \
    --n 1000000 \
    --workers 32 \
    --out results/phase2_quadruplet.auto.json
date
echo "=== Completed Quadruplet ==="

