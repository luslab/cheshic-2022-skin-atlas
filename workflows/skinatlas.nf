/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

// Create summary input parameters map for reporting
def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters in specialised library
WorkflowSkinatlas.initialise(params, log)

// Check input path parameters to see if the files exist if they have been specified
checkPathParamList = [
    params.input
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters that cannot be checked in the groovy lib as we want a channel for them
if (params.input) { ch_input = file(params.input) } else { exit 1, "Input samplesheet not specified!" }

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

/*
========================================================================================
    INIALISE PARAMETERS AND OPTIONS
========================================================================================
*/

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

include { R_SEURAT_CREATE_OBJECTS } from "../modules/local/r_seurat_create_objects"

/*
========================================================================================
    IMPORT SUBWORKFLOWS
========================================================================================
*/

include { INPUT_CHECK } from "../subworkflows/local/input_check"

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow SKINATLAS {

    INPUT_CHECK (
        ch_input
    )
    //INPUT_CHECK.out.folders | view

    R_SEURAT_CREATE_OBJECTS (
        INPUT_CHECK.out.folders
    )
}

////////////////////////////////////////////////////
/* --              COMPLETION EMAIL            -- */
////////////////////////////////////////////////////

// workflow.onComplete {
//     NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
//     NfcoreTemplate.summary(workflow, params, log)
// }

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
