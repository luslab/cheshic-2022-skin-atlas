process R_SEURAT_CREATE_OBJECTS {
    tag "$meta.id"
    label 'process_min'

    container "chrischeshire/skinatlas-r:latest"

    input:
    tuple val(meta), path(folder)

    output:
    tuple val(meta), path('*.rds'), emit: rds

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"

    """
    Rscript $task.ext.script \
    --cores $task.cpus \
    --id $meta.id \
    --output $prefix \
    --folder $folder \
    --mincells $task.ext.mincells \
    --minfeatures $task.ext.minfeatures \
    $args
    """
}
