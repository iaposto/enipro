#!/bin/bash
#SBATCH --job-name=admixture-1.3.0-case
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=30:00:00

# Load module
module load admixture/1.3.0

echo "[$(date)] Running ADMIXTURE K=2-30, 10-fold CV, 16 threads, seed=42"

for K in $(seq 2 30); do
    echo "[$(date)] K=${K}"
    admixture --cv=10 -j16 -s 42 graega_ldpruned.bed "$K" | tee "log${K}.out"
done
