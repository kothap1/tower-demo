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
workflow {
    Hello()
    S3test([params.redsheet, file(params.redsheet)])
    parseManifests([params.redsheet, file(params.redsheet), file(params.manifestdir), params.fastq_rootdir])
}