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
    ls -lrth
    pwd
    ls -lrthL
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
    ls -lrth
    pwd
    ls -lrthL
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
    ls -lrth
    pwd
    ls -lrthL
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
    """
}
process bwamem {
    container params.container_bwa_samtools
    // publishDir params.publish_dir
    cpus params.bwamem_threads
    memory params.disk_bwamem+' GB'
    tag {"$sampleID"}
    input:
    tuple val(sampleID), file(read1), file(read2) file(reference), file(params.reference+'.bwt'), file(params.reference+'.pac'), file(params.reference+'.sa') , file(params.reference+'.amb'), file(params.reference+'.ann')
    output:
    tuple val(sampleID), file("*bam")
    script:
    """
    echo bwa-mem ${sampleID}
    ls -lrth
    pwd
    ls -lrthL
    bwa mem -t ${params.bwamem_threads} \
        $reference $read1 $read2 | samtools view -hb | samtools sort -o ${sampleID}.sorted.bam
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
    ls -lrth
    pwd
    ls -lrthL
    sambamba merge \
        --nthreads ${params.sambamba_merge_threads} \
        ${sampleID}.merged.bam $bam 
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
    ls -lrth
    pwd
    ls -lrthL
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
    ls -lrth
    pwd
    ls -lrthL
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
    ls -lrth
    pwd
    ls -lrthL
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
    ls -lrth
    pwd
    ls -lrthL
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
    ls -lrth
    pwd
    ls -lrthL
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
    ls -lrth
    pwd
    ls -lrthL
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
    ch_reads = tuple(params.sid,[file(params.r11), file(params.r12)],[file(params.r21), file(params.r22)])
    mergeFastqs(ch_reads)
    ch_samples = tuple(params.fsid, file(params.fr1), file(params.fr2))
    fastp(ch_samples)
    ch_samples = tuple(params.bsid, file(params.br1), file(params.br2), file(params.reference), file(params.reference+'.bwt'), file(params.reference+'.pac'), file(params.reference+'.sa'), file(params.reference+'.amb'), file(params.reference+'.ann'))
    bwamem(ch_samples)
    sambamba_merge(tuple(params.bmesid, [file(params.bme1), file(params.bme2), file(params.bme3), file(params.bme4)]))
    sambamba_markdup(tuple(params.bmdsid, file(params.bmd)))
    picard_CollectInsertSizeMetrics(tuple(params.psid, file(params.pb)))
    picard_CollectMultipleMetrics(tuple(params.psid, file(params.pb), file(params.reference)))
    mosdepth(tuple(params.psid, file(params.pb), file(params.pbi)))
}