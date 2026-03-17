#!/usr/bin/env nextflow

nextflow.enable.dsl=2

workflow {
    // 1. Create channels for input files
    reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    ref_ch = Channel.fromPath(params.references, checkIfExists: true)

    // 2. Process: FastQC pre-trimming
    FASTQC_PRE(reads_ch)

    // 3. Process: Trim and Merge
    TRIM_AND_MERGE(reads_ch)

    // 4. Process: FastQC post-trimming
    FASTQC_POST(TRIM_AND_MERGE.out.merged_reads)

    // 5. Process: Index Reference
    INDEX_REF(ref_ch)

    // 6. Process: Align Reads
    ALIGN_READS(
        TRIM_AND_MERGE.out.merged_reads, 
        INDEX_REF.out.ref.collect(),
        INDEX_REF.out.index_files.collect()
    )

    // 7. Process: Filter BAM
    FILTER_BAM(ALIGN_READS.out.bam)

    // 8. Process: Flagstat pre-filtering
    FLAGSTAT_PRE(ALIGN_READS.out.bam)

    // 9. Process: Flagstat post-filtering
    FLAGSTAT_POST(FILTER_BAM.out.filtered_bam.map { it -> [it[0], it[1]] })

    // 10. Process: Call Variants
    CALL_VARIANTS(FILTER_BAM.out.filtered_bam, INDEX_REF.out.ref.collect(), INDEX_REF.out.fai.collect())

    // 11. Process: Parse and Plot
    PARSE_AND_PLOT(
        CALL_VARIANTS.out.vcf, 
        FILTER_BAM.out.filtered_bam, 
        INDEX_REF.out.ref.collect()
    )

    // 12. Process: Coverage Depth
    GET_COVERAGE(FILTER_BAM.out.filtered_bam)

    // 13. Process: Final Report
    report_inputs_ch = PARSE_AND_PLOT.out.tsv
        .join(PARSE_AND_PLOT.out.details_tsv)
        .join(PARSE_AND_PLOT.out.pdf)
        .join(PARSE_AND_PLOT.out.png)
        .join(GET_COVERAGE.out.png)
        .join(FLAGSTAT_POST.out.txt)

    FINAL_REPORT(report_inputs_ch)

    // 14. Process: MultiQC
    MULTIQC(
        FASTQC_PRE.out.zip.collect().ifEmpty([]),
        FASTQC_POST.out.zip.collect().ifEmpty([]),
        TRIM_AND_MERGE.out.json.collect().ifEmpty([]),
        FLAGSTAT_PRE.out.txt.map{ it[1] }.collect().ifEmpty([]),
        FLAGSTAT_POST.out.txt.map{ it[1] }.collect().ifEmpty([])
    )

    // 15. Workflow completion message
    workflow.onComplete = {
        println "Pipeline completed at: ${workflow.complete}"
        println "Execution status: ${workflow.success ? 'OK' : 'failed'}"
        println "Duration: ${workflow.duration}"
    }
}

process FLAGSTAT_PRE {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_pre_filtering.flagstat"), emit: txt

    script:
    """
    samtools flagstat $bam > ${sample_id}_pre_filtering.flagstat
    """
}

process FLAGSTAT_POST {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_post_filtering.flagstat"), emit: txt

    script:
    """
    samtools flagstat $bam > ${sample_id}_post_filtering.flagstat
    """
}

process FASTQC_PRE {
    tag "$sample_id"
    cpus 2

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*.zip", emit: zip
    path "*.html", emit: html

    script:
    """
    fastqc -t $task.cpus $reads
    """
}

process FASTQC_POST {
    tag "$sample_id"
    cpus 2

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*.zip", emit: zip
    path "*.html", emit: html

    script:
    """
    fastqc -t $task.cpus $reads
    """
}

process TRIM_AND_MERGE {
    tag "$sample_id"
    cpus 4

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_merged.fastq.gz"), emit: merged_reads
    path "${sample_id}_fastp.json"                            , emit: json
    path "${sample_id}_fastp.html"                            , emit: html
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
    path ref, emit: ref
    path "${ref}.fai", emit: fai
    path "${ref}*", emit: index_files

    script:
    """
    bwa index $ref
    samtools faidx $ref
    """
}

process ALIGN_READS {
    tag "$sample_id"
    cpus 4

    input:
    tuple val(sample_id), path(merged_reads)
    path ref
    path index_files

    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam"), emit: bam

    script:
    """
    bwa mem -t $task.cpus $ref $merged_reads | \\
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
    path fai

    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), path("${sample_id}.vcf.gz.csi"), emit: vcf

    script:
    """
    bcftools mpileup -B -a AD,DP -d 10000 -L 10000 -f $ref $bam | \\
        bcftools call -mv | \\
        bcftools filter -i 'MAX(FMT/AD[0:1-]) >= ${params.min_variant_reads}' -Oz -o ${sample_id}.vcf.gz
    bcftools index ${sample_id}.vcf.gz
    """
}

process PARSE_AND_PLOT {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(vcf), path(vcf_idx)
    tuple val(sample_id), path(bam), path(bai)
    path ref

    output:
    tuple val(sample_id), path("${sample_id}_Alleles_frequency_table.tsv")      , emit: tsv
    tuple val(sample_id), path("${sample_id}_Mutation_Details_table.tsv")       , emit: details_tsv
    tuple val(sample_id), path("${sample_id}_Editing_Frequencies.pdf")         , emit: pdf
    tuple val(sample_id), path("${sample_id}_Editing_Frequencies.png")         , emit: png

    script:
    """
    parse_edits.py \\
        --vcf $vcf \\
        --bam $bam \\
        --ref $ref \\
        --sgrna "${params.sgrna_seq}" \\
        --window ${params.quant_window} \\
        --min_freq ${params.min_edit_freq} \\
        --output ${sample_id}_Alleles_frequency_table.tsv \\
        --output_details ${sample_id}_Mutation_Details_table.tsv

    plot_summaries.py \\
        --input ${sample_id}_Alleles_frequency_table.tsv \\
        --output ${sample_id}_Editing_Frequencies.pdf
    """
}

process GET_COVERAGE {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}_coverage.png"), emit: png

    script:
    """
    samtools depth -a $bam > ${sample_id}_depth.txt
    plot_coverage.py --input ${sample_id}_depth.txt --output ${sample_id}_coverage.png --sample $sample_id
    """
}

process FINAL_REPORT {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(summary_tsv), path(details_tsv), path(edit_pdf), path(edit_png), path(cov_png), path(flagstat)

    output:
    path "${sample_id}_Final_Report.html", emit: html

    script:
    """
    generate_report.py \\
        --sample $sample_id \\
        --sgrna "${params.sgrna_seq}" \\
        --window ${params.quant_window} \\
        --min_edit_freq ${params.min_edit_freq} \\
        --min_variant_reads ${params.min_variant_reads} \\
        --flagstat $flagstat \\
        --edit_plot $edit_png \\
        --summary_tsv $summary_tsv \\
        --details_tsv $details_tsv \\
        --cov_plot $cov_png \\
        --output ${sample_id}_Final_Report.html
    """
}

process MULTIQC {
    input:
    path 'fastqc_pre/*'
    path 'fastqc_post/*'
    path 'fastp/*'
    path 'flagstat_pre/*'
    path 'flagstat_post/*'

    output:
    path "multiqc_report.html", emit: report

    script:
    """
    multiqc .
    """
}
