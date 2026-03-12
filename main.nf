#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Validation of required parameters
if (params.references == null) { exit 1, "Please provide a reference FASTA file via --references" }
if (params.sgrna_seq == null) { exit 1, "Please provide an sgRNA sequence via --sgrna_seq" }

workflow {
    // 1. Create channels for input files
    reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    ref_ch = Channel.fromPath(params.references, checkIfExists: true)

    // 2. Process: Trim and Merge
    TRIM_AND_MERGE(reads_ch)

    // 3. Process: Index Reference
    INDEX_REF(ref_ch)

    // 4. Process: Align Reads
    ALIGN_READS(TRIM_AND_MERGE.out.merged_reads, INDEX_REF.out.index_files.collect())

    // 5. Process: Filter BAM
    FILTER_BAM(ALIGN_READS.out.bam)

    // 6. Process: Call Variants
    CALL_VARIANTS(FILTER_BAM.out.filtered_bam, ref_ch)

    // 7. Process: Parse and Plot
    PARSE_AND_PLOT(
        CALL_VARIANTS.out.vcf, 
        FILTER_BAM.out.filtered_bam, 
        ref_ch
    )
}

process TRIM_AND_MERGE {
    tag "$sample_id"
    cpus 4

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_merged.fastq.gz"), emit: merged_reads
    path "versions.yml"                                      , emit: versions

    script:
    """
    fastp \\
        --in1 ${reads[0]} \\
        --in2 ${reads[1]} \\
        --merge \\
        --merged_out ${sample_id}_merged.fastq.gz \\
        --html ${sample_id}_fastp.html \\
        --json ${sample_id}_fastp.json \\
        --thread $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
    END_VERSIONS
    """
}

process INDEX_REF {
    tag "${ref.baseName}"

    input:
    path ref

    output:
    path "${ref}*", emit: index_files

    script:
    """
    bwa index $ref
    """
}

process ALIGN_READS {
    tag "$sample_id"
    cpus 4

    input:
    tuple val(sample_id), path(merged_reads)
    path index_files

    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam"), emit: bam

    script:
    def ref_fasta = index_files.find { it.name.endsWith('.fa') || it.name.endsWith('.fasta') }
    if (!ref_fasta) {
        // Find by excluding bwa extensions
        ref_fasta = index_files.find { ! (it.name.endsWith('.amb') || it.name.endsWith('.ann') || it.name.endsWith('.bwt') || it.name.endsWith('.pac') || it.name.endsWith('.sa')) }
    }
    """
    bwa mem -t $task.cpus $ref_fasta $merged_reads | \\
        samtools sort -@ $task.cpus -o ${sample_id}_sorted.bam -
    """
}

process FILTER_BAM {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_filtered.bam"), path("${sample_id}_filtered.bam.bai"), emit: filtered_bam

    script:
    """
    samtools view -h -q ${params.mapq_threshold} -F 4 -o ${sample_id}_filtered.bam $bam
    samtools index ${sample_id}_filtered.bam
    """
}

process CALL_VARIANTS {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam), path(bai)
    path ref

    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), emit: vcf

    script:
    """
    bcftools mpileup -f $ref $bam | \\
        bcftools call -m --ploidy ${params.ploidy} -Oz -o ${sample_id}.vcf.gz
    bcftools index ${sample_id}.vcf.gz
    """
}

process PARSE_AND_PLOT {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(vcf)
    tuple val(sample_id), path(bam), path(bai)
    path ref

    output:
    path "${sample_id}_Alleles_frequency_table.tsv"      , emit: tsv
    path "${sample_id}_Editing_Frequencies.pdf"         , emit: pdf

    script:
    """
    parse_edits.py \\
        --vcf $vcf \\
        --bam $bam \\
        --ref $ref \\
        --sgrna "${params.sgrna_seq}" \\
        --window ${params.quant_window} \\
        --output ${sample_id}_Alleles_frequency_table.tsv

    plot_summaries.py \\
        --input ${sample_id}_Alleles_frequency_table.tsv \\
        --output ${sample_id}_Editing_Frequencies.pdf
    """
}
