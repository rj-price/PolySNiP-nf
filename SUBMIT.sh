#!/usr/bin/env bash
#SBATCH -J polysnip-nf
#SBATCH --partition=short
#SBATCH --mem=2G
#SBATCH --cpus-per-task=4

source activate polysnip-nf

nextflow run main.nf \
    --reads 'data/sample1_{R1,R2}*.fastq.gz' \
    --references 'refs/homoeologues.fasta' \
    --sgrna_seq 'GCAGAGAATACTTGCCGAGG' \
    --outdir 'sample1'