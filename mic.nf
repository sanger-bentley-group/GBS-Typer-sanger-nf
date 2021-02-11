/*
 * Nextflow pipeline for GBS MIC prediction.
 *
 */

// Enable DSL 2
nextflow.enable.dsl=2

// Import modules
include {printHelp} from './modules/mic_help.nf'
include {target2mic} from './modules/target2mic.nf'


// Help message
if (params.help){
    printHelp()
    exit 0
}

// Check if output specified
if (params.output == ""){
    println("Please specify and output prefix with --output.")
    println("Print help with --help")
    System.exit(1)
}


// Create results directory if it doesn't already exist
results_dir = file(params.results_dir)
results_dir.mkdir()

// Output files
params.alleles_variants_out = "${params.output}_drug_cat_alleles_variants.txt"
params.existing_pbp_alleles_out = "${params.output}_existing_pbp_alleles.txt"
params.target2mic_out = "${params.output}_mic_predictions.txt"


// Main Workflow
workflow {

     main:
        target2mic(
                file("${results_dir}/${params.alleles_variants_out}", checkIfExists: true),
                file("${results_dir}/${params.existing_pbp_alleles_out}", checkIfExists: true),
                file("${results_dir}/${params.target2mic_out}"))
}
