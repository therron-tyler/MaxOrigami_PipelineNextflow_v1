#!/usr/bin/env bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 08:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tyler.therron@northwestern.edu
#SBATCH --output=%x.%j.out
#SBATCH --mem=40gb
#SBATCH --job-name=corigami_mm10
#SBATCH -N 1
#SBATCH -n 10

set -euo pipefail

module load samtools pigz bedtools
# mamba activate /projects/b1063/TylerT/C_oriGAMZzz

#export PYTHONNOUSERSITE=1
#unset PYTHONPATH


# Pass through any CLI args to the runner
./Corigami_predict_mm10.sh "$@"
