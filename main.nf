nextflow.enable.dsl = 2

process Hello {
    container 'ubuntu'
    script:
    """
    echo Hello
    """
}

workflow {
    Hello()
}