# PolySNiP-nf: Polyploid CRISPR Amplicon Genotyping

A reproducible Nextflow pipeline designed for accurate alignment and quantification of CRISPR-mediated mutations in polyploid genomes (e.g., octoploid strawberry).

## Features
- **Adapter Trimming & Merging**: Uses `fastp` for quality filtering and merging paired-end reads.
- **Local Alignment**: Uses `BWA-MEM` for robust alignment against highly similar homoeologous sequences.
- **MAPQ Filtering**: Isolates reads that map uniquely to specific homoeologues.
- **Polyploid-aware Variant Calling**: Leverages `bcftools` with adjustable ploidy.
- **CRISPR Quantification**: Custom Python scripts for per-homoeologue editing efficiency reporting and visualization.

## Quick Start

1. Install [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) and [Conda/Docker](https://docs.nextflow.io/en/latest/container.html).

2. Run the pipeline:

```bash
nextflow run main.nf \
    --reads 'data/*_{R1,R2}*.fastq.gz' \
    --references 'refs/homoeologues.fasta' \
    --sgrna_seq 'AGAGTTGATGCTTCTGGGAT' \
    --ploidy 8 \
    -profile conda
```

## Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--reads` | Glob pattern for raw Illumina FASTQ files | `data/*_{R1,R2}.fastq.gz` |
| `--references` | Multi-FASTA containing all homoeologue sequences | `null` (Required) |
| `--sgrna_seq` | The target sgRNA sequence (guide) | `null` (Required) |
| `--quant_window`| Quantification window around the cut site (bp) | `10` |
| `--ploidy` | Ploidy level for variant calling | `8` |
| `--mapq_threshold`| MAPQ threshold for uniquely mapping reads | `20` |
| `--outdir` | Output directory for results | `results` |

## Output Structure

- `alignments/`: Sorted, MAPQ-filtered `.bam` and `.bam.bai` files.
- `variants/`: Raw `.vcf.gz` files from `bcftools`.
- `quantification/`: Tabular summary of editing frequencies.
- `reports/`: Visualizations showing editing efficiency per homoeologue.
- `preprocessing/`: `fastp` reports and trimmed/merged reads.

## Dependencies
The pipeline manages dependencies via the `environment.yml` file using Conda or Docker. Key tools include:
- `fastp`
- `bwa`
- `samtools`
- `bcftools`
- `pysam`, `pandas`, `matplotlib`, `seaborn` (Python)
