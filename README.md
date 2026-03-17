# PolySNiP-nf: Polyploid CRISPR Amplicon Genotyping

A reproducible Nextflow pipeline designed for accurate alignment and quantification of CRISPR-mediated mutations in polyploid genomes.

## Features
- **Comprehensive QC**: FastQC and MultiQC for read quality monitoring.
- **Adapter Trimming & Merging**: Uses `fastp` for quality filtering and merging paired-end reads.
- **Local Alignment**: Uses `BWA-MEM` for robust alignment against highly similar homoeologous sequences.
- **Sensitive Variant Calling**: Optimised `bcftools` parameters for high-coverage amplicon sequencing.
- **CRISPR Quantification**: Advanced SNP and INDEL detection within a custom window around the sgRNA cut site.
- **Automated Reporting**: Generates a self-contained HTML report per sample summarising all metrics and results.

---

## Setup & Installation

### Prerequisites
- [Conda](https://docs.conda.io/en/latest/miniconda.html) or Miniconda/Mamba.

### Initial Setup
Pull the repo and run the provided setup script to create the environment and set file permissions:

```bash
git clone https://github.com/rj-price/PolySNiP-nf.git && cd PolySNiP-nf
bash SETUP.sh
```

This script will:
1. Make all Python scripts in `bin/` executable.
2. Create the `polysnip-nf` Conda environment.
3. Prepare the `data/` and `refs/` directories.

---

## Usage

### 1. Prepare Inputs
- Place your raw paired-end FASTQ files in `data/`.
- Place your reference Multi-FASTA (containing homoeologue sequences) in `refs/`.

### 2. Run on Slurm (HPC)
A `SUBMIT.sh` script is provided for Slurm job submission. 

Set the `--reads`, `--references`, `--sgrna_seq` and `--outdir` parameters inside `SUBMIT.sh`, then run:

```bash
sbatch SUBMIT.sh
```

### 3. Manual Execution
To run locally or interactively:

```bash
nextflow run main.nf \
    --references 'refs/your_ref.fasta' \
    --sgrna_seq 'YOUR_GUIDE_SEQUENCE'
```

---

## Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--reads` | Glob pattern for raw FASTQ files | `data/*_{R1,R2}.fastq.gz` |
| `--references` | Multi-FASTA containing all homoeologues | `null` (Required) |
| `--sgrna_seq` | The target sgRNA sequence (guide) | `null` (Required) |
| `--quant_window`| Window size around cut site for quantification (bp) | `10` |
| `--min_edit_freq`| Min frequency (%) to report a variant in details table | `1.0` |
| `--min_variant_reads`| Min reads to support a variant call in the VCF | `10` |
| `--mapq_threshold`| MAPQ threshold for uniquely mapping reads | `20` |
| `--outdir` | Output directory for results | `results` |

---

## Outputs

All results are organised in the `--outdir` directory:

### **Primary Results**
- `*_Final_Report.html`: **The main output.** A comprehensive, self-contained report for each sample containing:
    - Pipeline parameters and execution time.
    - Alignment statistics.
    - Coverage depth histograms across homoeologues.
    - Editing frequency plots and summary tables.
    - Detailed mutation tables (filtered by frequency).

### **Data Tables & Plots**
- `quantification/`:
    - `*_Alleles_frequency_table.tsv`: Summary of modified vs. unmodified reads per homoeologue.
    - `*_Mutation_Details_table.tsv`: Detailed list of every SNP/INDEL found within the window.
- `reports/`:
    - `*_Editing_Frequencies.pdf`: Bar charts of editing efficiency.
    - `*_Editing_Frequencies.png`: PNG version used in the HTML report.

### **Intermediate Files & QC**
- `qc/`: Quality control reports from FastQC, Flagstat, and MultiQC.
- `alignments/`: Filtered and indexed BAM files (`*_filtered.bam`).
- `variants/`: Variant call files (`*.vcf.gz`) containing SNP/INDEL calls supported by >= 10 reads.
- `preprocessing/`: `fastp` reports and the merged FASTQ files.

---

## Methodology Notes
- **Cut Site**: Automatically calculated as 3bp upstream of the PAM (assumed 3' of the provided guide).
- **SNP Detection**: Unlike standard CIGAR-based tools, this pipeline performs a base-by-base comparison against the reference to ensure even single-nucleotide mismatches are quantified.
- **High-Coverage Sensitivity**: `bcftools` parameters are tuned to prevent sub-sampling in high-depth (e.g., >5000x) amplicon data.
