/*
 * Check input samplesheet and get read channels
 */

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )

    SAMPLESHEET_CHECK.out.csv
        .splitCsv ( header:true, sep:"," )
        .map { get_samplesheet_paths(it) }
        .set { ch_folders }

    emit:
    folders = ch_folders // channel: [ val(meta), [ folders ] ]
    versions = SAMPLESHEET_CHECK.out.versions
}

def get_samplesheet_paths(LinkedHashMap row) {
    def meta = [:]
    meta.id = row.id
    array = [ meta, row.folder ]

    return array
}
