#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 40:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tyler.therron@northwestern.edu
#SBATCH --output=%x.%j.out
#SBATCH --mem=60gb
#SBATCH --job-name=maxATAC_Corigami_nf_test
#SBATCH -N 1
#SBATCH -n 10

module load nextflow/25.04.0

cd /home/ttm3567/rootdir_scratch/20250806_Corigami/maxATAC_Corigami_NextflowPipeline

nextflow run TAD_PredictionPipeline_NextflowMain_v1-1.nf \
  --input_pattern '../BAM/*.bam' \
  --chroms 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19' \
  --outdir BowdishATAC_Ly6cHigh_YoungOld_TNFKO_HiC_Prediction \
  --run_benchmark true
