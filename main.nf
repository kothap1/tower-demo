nextflow.enable.dsl = 2

process Hello {
    container 'ubuntu'
    script:
    """
    echo Hello
    """
}

process S3test {
    container 'ubuntu'
    input:
    tuple val(loc), path(redsheet)
    script:
    """
    echo ${loc}
    echo ${redsheet}
    cat ${redsheet}
    """
}
process parseManifests {
    container params.container_python
    publishDir params.publish_dir
    tag {"$samples"}

    input:
    tuple val(samples), file(redsheet), path(manifestdir), val(fastq_rootdir)

    output:
    path('*csv')

    script:
    """
    echo Parse manifests ${redsheet}
    ls -lrthL
    echo ===
    ls -lrth
    df -h
    parse_manifests.py \
        --processing per-lane \
        --redsheet ${redsheet} \
        --manifestdir ${manifestdir} \
        --fastq_rootdir ${fastq_rootdir}

    """
}
process mergeFastqs {
    // publishDir params.publish_dir
    container 'ubuntu'
    memory params.disk_merge_fastqs+' GB'
    tag {"$sampleID"}
    input:
    tuple val(sampleID), file("${sampleID}_R1_*.fastq.gz"), file("${sampleID}_R2_*.fastq.gz")
    output:
    tuple val(sampleID), path("*R1.merged.fastq.gz"),  path("*R2.merged.fastq.gz")
    script:
    """
    echo Merge fastqs ${sampleID}
    ls -lrthL
    echo ===
    ls -lrth
    df -h
    cat ${sampleID}_R1_*.fastq.gz > ${sampleID}_R1.merged.fastq.gz
    #rm -f `readlink ${sampleID}_R1_*.fastq.gz`
    cat ${sampleID}_R2_*.fastq.gz > ${sampleID}_R2.merged.fastq.gz
    #rm -f `readlink ${sampleID}_R2_*.fastq.gz`
    """
}
process fastp {
    container params.container_fastp
    publishDir params.publish_dir, pattern: "*.json"
    cpus params.fastp_threads
    memory params.disk_fastp+' GB'
    tag {"$sampleID"}
   
    input:
    tuple val(sampleID), file(read1), file(read2)

    output:
    tuple val(sampleID), file("${sampleID}_R1.trimmed.fq.gz"), file("${sampleID}_R2.trimmed.fq.gz"), emit: trimmed_fqs
    path("*json"), emit: fastp_qc

    script:
    """
    echo Fastp ${sampleID}
    ls -lrthL
    echo ===
    ls -lrth
    df -h
    fastp -i $read1 -I $read2 \
        -o ${sampleID}_R1.trimmed.fq.gz -O ${sampleID}_R2.trimmed.fq.gz \
        -h ${sampleID}.html -j ${sampleID}.json -R "fastp report for sample: ${sampleID}" \
        --thread ${params.fastp_threads} \
        --cut_window_size 5 \
        --qualified_quality_phred 20 \
        --length_required 30 \
        --detect_adapter_for_pe \
        --trim_poly_x \
        --trim_poly_g \
        --cut_tail
    rm -f `readlink $read1`
    rm -f `readlink $read2`
    """
}
process bwamem {
    container params.container_bwa_samtools
    // publishDir params.publish_dir
    cpus params.bwamem_threads
    memory params.disk_bwamem+' GB'
    tag {"$sampleID"}
    input:
    tuple val(sampleID), file(read1), file(read2) //, file(params.reference+'.amb'), file(params.reference+'.ann'), file(params.reference+'.bwt'), file(params.reference+'.pac'), file(params.reference+'.sa')
    output:
    tuple val(sampleID), file("*bam")
    script:
    """
    echo bwa-mem ${sampleID}
    ls -lrthL
    echo ===
    ls -lrth
    df -h
    bwa mem -t ${params.bwamem_threads} \
        ${params.reference} $read1 $read2 | samtools view -hb | samtools sort -o ${sampleID}.sorted.bam
    rm -f `readlink $read1`
    rm -f `readlink $read2`
    """
}
process sambamba_merge {
    container params.container_sambamba
    // publishDir params.publish_dir
    cpus params.sambamba_merge_threads
    memory params.disk_sambamba_merge+' GB'
    tag {"$sampleID"}
    input:
    tuple val(sampleID), file(bam)
    output:
    tuple val(sampleID), file("*.merged.bam"), file("*bai")
    script:
    """
    echo sambamba-merge ${sampleID}
    ls -lrthL
    echo ===
    ls -lrth
    df -h
    sambamba merge \
        --nthreads ${params.sambamba_merge_threads} \
        ${sampleID}.merged.bam $bam 
    rm -f `readlink $bam`
    """
}
process sambamba_markdup {
    container params.container_sambamba
    publishDir params.publish_dir, pattern: "*.{bam,bai}"
    cpus params.sambamba_markdup_threads
    memory params.disk_sambamba_markdup+' GB'
    tag {"$sampleID"}
    input:
    tuple val(sampleID), file(bam)
    output:
    tuple val(sampleID), file("*.markeddup.bam"), file("*bai")
    script:
    """
    echo sambamba-markdup
    ls -lrthL
    echo ===
    ls -lrth
    df -h
    sambamba markdup \
        --overflow-list-size 200000 \
        $bam ${sampleID}.markeddup.bam

    """
}
process picard_CollectInsertSizeMetrics {
    container params.container_picard
    publishDir params.publish_dir, pattern: "*.txt"
    memory params.disk_picard_CISM+' GB'
    tag {"$sampleID"}
    input:
    tuple val(sampleID), path(bam)
    output:
    tuple val(sampleID), path('*.{txt,pdf}')
    script:
    """
    echo Picard CISM ${sampleID}
    ls -lrthL
    echo ===
    ls -lrth
    df -h
    picard CollectInsertSizeMetrics  \
    -I $bam \
    -O ${sampleID}.insert_size_metrics.txt \
    -H ${sampleID}.insert_size_histogram.pdf \
    -VALIDATION_STRINGENCY SILENT
    """
}
process picard_CollectMultipleMetrics {
    container params.container_picard
    publishDir params.publish_dir, pattern: "*metrics"
    memory params.disk_picard_CMM+' GB'
    tag {"$sampleID"}
    input:
    tuple val(sampleID), path(bam), path(reference)
    output:
    tuple val(sampleID), path('*')
    script:
    // Read list of metrics-to-collect from the config and transform appropriately to use in the final command
    programs_list = []
    for (param in params.MULTI_METRICS_PROGRAMS.split(',')) { programs_list << "--PROGRAM $param" }
    PROGRAMS = programs_list.join(' ').trim()
    // reference = file(params.reference)
    """
    echo picard CMM ${sampleID}
    ls -lrthL
    echo ===
    ls -lrth
    df -h
    picard CollectMultipleMetrics  \
    -I $bam \
    -O ${sampleID} \
    -R ${reference} \
    --PROGRAM null \
    ${PROGRAMS} \
    -VALIDATION_STRINGENCY SILENT
    """
}
process mosdepth {
    container params.container_mosdepth
    publishDir params.publish_dir, pattern: "*.{html,xlsx}"
    cpus params.mosdepth_threads
    memory params.disk_mosdepth+' GB'
    tag {"$sampleID"}
    input:
    tuple val(sampleID), path(bam), path(bai)
    output:
    tuple val(sampleID), path("*.txt")
    script:
    """
    echo mosdepth ${sampleID}
    ls -lrthL
    echo ===
    ls -lrth
    df -h
    mosdepth --no-per-base -F 1796 -i 2 $sampleID $bam
    """
}
process multiqc {
    container params.container_multiqc
    publishDir params.publish_dir, pattern: "*.html"
    memory params.disk_multiqc+' GB'
    tag {"$redsheet_name"}
    afterScript "for sample in ${samplenames}; do cp multiqc_report.html multiqc_\${sample}.html; done; rm multiqc_report.html"
    input:
    val redsheet_name
    path(multiqc_config)
    path("fastp/*")
    val samplenames
    output:
    path("*html")
    script:
    samplenames = samplenames.join(' ').trim()
    """
    echo MultiQC ${redsheet_name}
    ls -lrthL
    echo ===
    ls -lrth
    df -h
    multiqc .
    """
}
process collateQC {
    container params.container_collateqc
    publishDir params.publish_dir
	tag {"$redsheet_name"}
    afterScript "for sample in ${samplenames}; do cp collated_qc.html batchqc_\${sample}.html; cp collated_qc.xlsx batchqc_\${sample}.xlsx; done; rm collated_qc*"

    input:
	val redsheet_name
    path(notebook)
    path(redsheet)
    path("QC_metrics/fastp/*")
	path("QC_metrics/mosdepth/*")
	path("QC_metrics/picard_CollectInsertSizeMetrics/*")
	path("QC_metrics/picard_CollectMultipleMetrics/*")
    val samplenames

    output:
	val true, emit: collation_completed
	file('*{html,xlsx}')

    script:
	samplenames = samplenames.join(' ').trim()
    """
    echo CollateQC ${redsheet}
    python --version
    echo ${redsheet_name}
    echo ${notebook}
    echo ${redsheet}
    echo ${samplenames}
	echo "_ --redsheet ${redsheet} --batch_path ./QC_metrics/ --xlsx collated_qc.xlsx" > .config_ipynb
	jupyter nbconvert --execute ${notebook} --to html --no-input --output collated_qc.html --output-dir .
    for sample in ${samplenames}; do cp collated_qc.html batchqc_\${sample}.html; cp collated_qc.xlsx batchqc_\${sample}.xlsx; done; rm collated_qc*
    ls QC_metrics
    ls -a
    ls QC_metrics/fastp/
    ls QC_metrics/mosdepth/
    ls QC_metrics/picard_CollectInsertSizeMetrics/
    ls QC_metrics/picard_CollectMultipleMetrics/
    df -h
    """
}
process generate_manifests {
    container params.container_collateqc
    publishDir params.publish_dir+'/manifests', pattern: "*.json"
    tag {"$redsheet"}

    input:
	val flag
    val user_id
    path(redsheet)
    path(manifestdir)
	each samplename

    output:
    path("*json")

    script:
    """
    echo Generate Manifests ${redsheet}
    ls -lrthL
    echo ===
    ls -lrth
    df -h
	generate_manifests_for_outputs.py \
	--redsheet ${redsheet} \
	--manifestdir ${manifestdir} \
	--samplename ${samplename} \
    --userid ${user_id}
    """
}

workflow {
    Hello()
    S3test([params.redsheet, file(params.redsheet)])
    parseManifests([params.redsheet, file(params.redsheet), file(params.manifestdir), params.fastq_rootdir])
    parseManifests.out
        .splitCsv(header: true, sep: '\t')
        .map { row-> tuple(row.sampleID, file(row.read1), file(row.read2)) }
        .groupTuple(by: [0])
        .set { ch_reads }
    ch_samples = mergeFastqs(ch_reads)
    fastp(ch_samples)
    bwamem(fastp.out.trimmed_fqs) //.map { row -> row +[file(params.reference+'.amb'), file(params.reference+'.ann'), file(params.reference+'.bwt'), file(params.reference+'.pac'), file(params.reference+'.sa')]})
    bwamem.out
        .map { it -> [ (it[0] =~ /^(.+?)_S\d+(?=_(?:L00[0-4]|\d{4,}))/)[0][0], it[1]] }
        .groupTuple(by: [0])
        .set { ch_merge_bams }
    sambamba_merge(ch_merge_bams)
    sambamba_markdup(sambamba_merge.out.map { it -> [it[0], it[1]] })
    ch_markdup_bam = sambamba_markdup.out
    picard_CollectInsertSizeMetrics(ch_markdup_bam.map {it -> [it[0], it[1]]})
    picard_CollectMultipleMetrics(ch_markdup_bam.map {it -> [it[0], it[1], file(params.reference)]})
    mosdepth(ch_markdup_bam)
    multiqc_config = file("$projectDir/multiqc_config.yml")
    multiqc (
        params.redsheet.split("\\.")[0],
        multiqc_config,
        fastp.out.fastp_qc.collect(),
        ch_markdup_bam.map {it -> it[0]}.toList()
    )
    generate_manifests(
        true,
        params.user_id,
        file(params.redsheet),
        file(params.manifestdir),
        ch_markdup_bam.map {it -> it[0]}.toList()
    )
}