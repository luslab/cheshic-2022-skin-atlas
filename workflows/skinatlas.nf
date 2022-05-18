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

include { R_SEURAT_CREATE_OBJECTS       } from "../modules/local/r_seurat_create_objects"
include { R_GENERAL as R_PREPROCESSING  } from "../modules/local/r_general"
include { R_GENERAL as R_INTEGRATION    } from "../modules/local/r_general"
include { R_GENERAL as R_INTEGRATION_QC } from "../modules/local/r_general"

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

    // Check input
    INPUT_CHECK (
        ch_input
    )
    //INPUT_CHECK.out.data | view

    // Create intial seurat objects
    if(!params.rds_mode) {
        R_SEURAT_CREATE_OBJECTS (
            INPUT_CHECK.out.data
        )
    }
    //R_SEURAT_CREATE_OBJECTS.out.rds | view

    // Create merged list of seurat objects
    R_SEURAT_CREATE_OBJECTS.out.rds
        .map { row -> [row[1]] }
        .collect()
        .map { row -> [[id:params.project_name], row] }
        .set { ch_counts }   
    //ch_counts | view

    // Pre-process data and merge
    R_PREPROCESSING (
        ch_counts
    )
    //R_PREPROCESSING.out.files | view

    // Integration
    R_INTEGRATION (
        R_PREPROCESSING.out.files
    )
    //R_INTEGRATION.out.files | view

    // Integration QC
    R_INTEGRATION_QC (
        R_INTEGRATION_QC.out.files
    )
    //R_INTEGRATION_QC.out.files | view
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
