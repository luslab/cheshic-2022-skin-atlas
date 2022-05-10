#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     GENOME PARAMETER VALUES
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */

// if (!params.fasta) {
//     params.bowtie2   = WorkflowMain.getGenomeAttribute(params, 'bowtie2')
// } else {
//     params.bowtie2   = null
// }

// params.fasta     = WorkflowMain.getGenomeAttribute(params, 'fasta')
// params.gtf       = WorkflowMain.getGenomeAttribute(params, 'gtf')
// params.gene_bed  = WorkflowMain.getGenomeAttribute(params, 'bed12')
// params.blacklist = WorkflowMain.getGenomeAttribute(params, 'blacklist')

// /*
// ========================================================================================
//     SPIKEIN GENOME PARAMETER VALUES
// ========================================================================================
// */

// if (!params.spikein_fasta) {
//     params.spikein_bowtie2 = WorkflowMain.getGenomeAttributeSpikeIn(params, 'bowtie2')
// } else {
//     params.spikein_bowtie2 = null
// }
// params.spikein_fasta   = WorkflowMain.getGenomeAttributeSpikeIn(params, 'fasta')

// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     VALIDATE & PRINT PARAMETER SUMMARY
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */

// WorkflowMain.initialise(workflow, params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SKINATLAS } from './workflows/skinatlas'

workflow CRICK_SKINATLAS {
    //SKINATLAS ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    CRICK_SKINATLAS ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
