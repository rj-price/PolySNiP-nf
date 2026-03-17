#!/usr/bin/env bash
#SBATCH -J polysnip-loop
#SBATCH --partition=short
#SBATCH --mem=4G
#SBATCH --cpus-per-task=8

# Load environment
source activate polysnip-nf

# Define sgRNAs
GRNA1="GAAACATAATGCTTCCATGT"
GRNA2="GCAGAGAATACTTGCCGAGG"

# Reference file
REFS="refs/homoeologues.fasta"

# Loop through all R1 files in data/ to identify sample names
for r1_file in data/*_R1*.fastq.gz; do
    # Extract sample name (assuming format sample_name_R1...)
    # This removes the path and the suffix starting with _R1
    sample_name=$(basename "$r1_file" | sed 's/_R1.*//')
    
    echo "Processing sample: $sample_name"

    # Run for gRNA1
    echo "  Running gRNA1: $GRNA1"
    nextflow run main.nf \
        --reads "data/${sample_name}_{R1,R2}*.fastq.gz" \
        --references "$REFS" \
        --sgrna_seq "$GRNA1" \
        --outdir "${sample_name}_grna1"

    # Run for gRNA2
    echo "  Running gRNA2: $GRNA2"
    nextflow run main.nf \
        --reads "data/${sample_name}_{R1,R2}*.fastq.gz" \
        --references "$REFS" \
        --sgrna_seq "$GRNA2" \
        --outdir "${sample_name}_grna2"

done

echo "Loop completion: All samples processed."
