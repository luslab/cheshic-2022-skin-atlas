/*
 * Check input samplesheet and get read channels
 */

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    samplesheet
        .splitCsv ( header:true, sep:"," )
        .map { get_samplesheet_paths(it) }
        .set { folders }

    emit:
    folders // channel: [ val(meta), [ folders ] ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def get_samplesheet_paths(LinkedHashMap row) {
    def meta = [:]
    meta.id = row.id
    array = [ meta, row.folder ]

    return array
}
