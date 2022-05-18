process R_GENERAL {
    tag "$meta.id"
    label 'process_high'

    container "chrischeshire/skinatlas-r:latest"

    input:
    tuple val(meta), path('input/*')

    output:
    tuple val(meta), path('*'), emit: files

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    Rscript $task.ext.script \
    --cores $task.cpus \
    --id $meta.id \
    --runtype nextflow \
    $args
    rm -r input
    """
}
