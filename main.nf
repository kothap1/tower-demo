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
    val loc
    path(redsheet)
    script:
    """
    echo ${loc}
    echo ${redsheet}
    cat ${redsheet}
    """
}
workflow {
    Hello()
    S3test([params.redsheet, file(params.redsheet)])
}