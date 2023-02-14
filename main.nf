nextflow.enable.dsl = 2

process Hello {
    script:
    """
    echo Hello
    """
}

workflow {
    Hello()
}