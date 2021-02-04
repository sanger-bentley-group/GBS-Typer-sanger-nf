/*
 * Nextflow pipeline for sero and resistance typing Group B Strep
 *
 */

// Enable DSL 2
nextflow.enable.dsl=2

// Import modules
include {printHelp} from './modules/help.nf'
include {serotyping} from './modules/serotyping.nf'
include {srst2_for_res_typing; split_target_RES_seq_from_sam_file; split_target_RES_sequences; freebayes} from './modules/res_alignments.nf'
include {res_typer} from './modules/res_typer.nf'
include {surface_typer} from './modules/surface_typer.nf'
include {srst2_for_mlst; get_mlst_allele_and_pileup} from './modules/mlst.nf'
include {get_pbp_genes; get_pbp_alleles} from './modules/pbp_typer.nf'
include {combine_results; finalise_surface_typer_results; finalise_pbp_existing_allele_results} from './modules/combine.nf'


// Help message
if (params.help){
    printHelp()
    exit 0
}

// Check if reads specified
if (params.run_sero_res | params.run_mlst | params.run_surfacetyper){
    if (params.reads == ""){
        println("Please specify reads with --reads.")
        println("Print help with --help")
        System.exit(1)
    }
    // Create read pairs channel
    Channel.fromFilePairs( params.reads, checkIfExists: true )
        .set { read_pairs_ch }
}

// Check if output specified
if (params.output == ""){
    println("Please specify and output prefix with --output.")
    println("Print help with --help")
    System.exit(1)
}

if (!params.run_sero_res && !params.run_surfacetyper && !params.run_mlst && !params.run_pbptyper){
    println("Please specify one or more pipelines to run.")
    println("Print help with --help")
    System.exit(1)
}

if (params.other_res_dbs.toString() != 'none'){

    // Check other resistance databases exist
    other_res_db_list = params.other_res_dbs.toString().tokenize(' ')
    for (db in other_res_db_list) {
        other_db = file(db, checkIfExists: true)
    }

    // Check number of resistance databases matches number of minimum coverage and max divergence parameters
    if (params.other_res_dbs.toString().tokenize(' ').size() != params.other_res_min_coverage.toString().tokenize(' ').size()){
        println("Number of --other_res_min_coverage values is not equal to the number of --other_res_dbs files")
        println("Please specify an equal number of values.")
        System.exit(1)
    }

    if (params.other_res_dbs.toString().tokenize(' ').size() != params.other_res_max_divergence.toString().tokenize(' ').size()){
        println("Number of --other_res_max_divergence values is not equal to the number of --other_res_dbs files")
        println("Please specify an equal number of values.")
        System.exit(1)
    }
}

// Check parameters are within range
if (params.gbs_res_min_coverage < 0 | params.gbs_res_min_coverage > 100){
    println("--gbs_res_min_coverage value not in range. Please specify a value between 0 and 100.")
    System.exit(1)
}

if (params.gbs_res_max_divergence < 0 | params.gbs_res_max_divergence > 100){
    println("--gbs_res_max_divergence value not in range. Please specify a value between 0 and 100.")
    System.exit(1)
}

other_res_min_coverage_list = params.other_res_min_coverage.toString().tokenize(' ')
for (other_res_min_coverage in other_res_min_coverage_list){
    if (other_res_min_coverage.toDouble() < 0 | other_res_min_coverage.toDouble() > 100){
        println("--other_res_min_coverage value(s) not in range. Please specify a value between 0 and 100.")
        System.exit(1)
    }
}

other_res_max_divergence_list = params.other_res_max_divergence.toString().tokenize(' ')
for (other_res_max_divergence in other_res_max_divergence_list){
    if (other_res_max_divergence.toDouble() < 0 | other_res_max_divergence.toDouble() > 100){
        println("--other_res_max_divergence value(s) not in range. Please specify a value between 0 and 100.")
        System.exit(1)
    }
}

if (params.restyper_min_read_depth < 0){
    println("--restyper_min_read_depth value not in range. Please specify a value of 0 or above.")
    System.exit(1)
}

if (params.serotyper_min_read_depth < 0){
    println("--serotyper_min_read_depth value not in range. Please specify a value of 0 or above.")
    System.exit(1)
}

if (params.mlst_min_coverage < 0 | params.mlst_min_coverage > 100){
    println("--mlst_min_coverage value not in range. Please specify a value between 0 and 100.")
    System.exit(1)
}

if (params.mlst_min_read_depth < 0){
    println("--mlst_min_read_depth value not in range. Please specify a value of 0 or above.")
    System.exit(1)
}

if (params.surfacetyper_min_coverage < 0 | params.surfacetyper_min_coverage > 100){
    println("--surfacetyper_min_coverage value not in range. Please specify a value between 0 and 100.")
    System.exit(1)
}

if (params.surfacetyper_max_divergence < 0 | params.surfacetyper_max_divergence > 100){
    println("--surfacetyper_max_divergence value not in range. Please specify a value between 0 and 100.")
    System.exit(1)
}

if (params.surfacetyper_min_read_depth < 0){
    println("--surfacetyper_min_read_depth value not in range. Please specify a value of 0 or above.")
    System.exit(1)
}

// Create tmp directory if it doesn't already exist
tmp_dir = file(params.tmp_dir)
tmp_dir.mkdir()

// Create results directory if it doesn't already exist
results_dir = file(params.results_dir)
results_dir.mkdir()

// Output files
params.sero_res_incidence_out = "${params.output}_serotype_res_incidence.txt"
params.variants_out =  "${params.output}_gbs_res_variants.txt"
params.alleles_variants_out = "${params.output}_drug_cat_alleles_variants.txt"
params.existing_pbp_alleles_out = "${params.output}_existing_pbp_alleles.txt"
params.surface_protein_incidence_out = "${params.output}_surface_protein_incidence.txt"
params.surface_protein_variants_out = "${params.output}_surface_protein_variants.txt"


// Resistance mapping with the GBS resistance database
workflow GBS_RES {

    take:
        reads

    main:

        gbs_res_typer_db = file(params.gbs_res_typer_db, checkIfExists: true)
        gbs_res_targets_db = file(params.gbs_res_targets_db, checkIfExists: true)

        // Split GBS target sequences from GBS resistance database into separate FASTA files per sequence
        split_target_RES_sequences(gbs_res_typer_db, gbs_res_targets_db)

        // Get path and name of GBS resistance database
        db_path = gbs_res_typer_db.getParent()
        db_name = gbs_res_typer_db.getName()

        // Map genomes to GBS resistance database using SRST2
        srst2_for_res_typing(reads, db_name, db_path, tmp_dir, 'GBS_RES', params.gbs_res_min_coverage, params.gbs_res_max_divergence)

        // Split sam file for each GBS target sequence
        split_target_RES_seq_from_sam_file(srst2_for_res_typing.out.bam_files, gbs_res_targets_db)

        // Get consensus sequence using freebayes
        freebayes(split_target_RES_seq_from_sam_file.out, split_target_RES_sequences.out, tmp_dir)

    emit:
        freebayes.out.id
}

// Resistance mapping with the other resistance databases
workflow OTHER_RES {

    take:
        reads

    main:
        // Map genomes to resistance database using SRST2
        srst2_for_res_typing(reads, params.other_res_dbs, file('.'), tmp_dir, 'OTHER_RES', params.other_res_min_coverage, params.other_res_max_divergence)

    emit:
        srst2_for_res_typing.out.id
}

// MLST pipeline
workflow MLST {

    take:
        reads

    main:
        // Run SRST2 MLST
        srst2_for_mlst(reads, file(params.mlst_allele_db, checkIfExists: true), file(params.mlst_definitions_db, checkIfExists: true), params.mlst_min_coverage)

        // Get new consensus allele and pileup data
        get_mlst_allele_and_pileup(srst2_for_mlst.out, params.mlst_min_read_depth, file(params.mlst_allele_db, checkIfExists: true))

    emit:
        get_mlst_allele_and_pileup.out
}

// PBP-1A allele typing pipeline
workflow PBP1A {

    take:
        pbp_typer_output

    main:
        // Run
        get_pbp_alleles(pbp_typer_output, 'GBS1A-1', file(params.gbs_blactam_1A_db, checkIfExists: true))

        // Output new PBP alleles to results directory
        get_pbp_alleles.out.new_pbp.subscribe { it ->
            it.copyTo(file("${results_dir}"))
        }

        // Combine existing PBP alleles results in one file
        finalise_pbp_existing_allele_results(get_pbp_alleles.out.existing_pbp)

    emit:
        // Emit existing PBP alleles for collection
        finalise_pbp_existing_allele_results.out
}

// PBP-2B allele typing pipeline
workflow PBP2B {

    take:
        pbp_typer_output

    main:
        // Run
        get_pbp_alleles(pbp_typer_output, 'GBS2B-1', file(params.gbs_blactam_2B_db, checkIfExists: true))

        // Output new PBP alleles to results directory
        get_pbp_alleles.out.new_pbp.subscribe { it ->
            it.copyTo(file("${results_dir}"))
        }

        // Combine existing PBP alleles results in one file
        finalise_pbp_existing_allele_results(get_pbp_alleles.out.existing_pbp)

    emit:
        // Emit existing PBP alleles for collection
        finalise_pbp_existing_allele_results.out
}

// PBP-2X allele typing pipeline
workflow PBP2X {

    take:
        pbp_typer_output

    main:
        // Run
        get_pbp_alleles(pbp_typer_output, 'GBS2X-1', file(params.gbs_blactam_2X_db, checkIfExists: true))

        // Output new PBP alleles to results directory
        get_pbp_alleles.out.new_pbp.subscribe { it ->
            it.copyTo(file("${results_dir}"))
        }

        // Combine existing PBP alleles results in one file
        finalise_pbp_existing_allele_results(get_pbp_alleles.out.existing_pbp)

    emit:
        // Emit existing PBP alleles for collection
        finalise_pbp_existing_allele_results.out
}

// Main Workflow
workflow {

    // Create new genotypes file
    Channel.fromPath( params.output )
        .set { output_ch }

    main:

        if (params.run_sero_res){

            // Serotyping Process
            serotyping(read_pairs_ch, file(params.serotyping_db, checkIfExists: true), params.serotyper_min_read_depth)

            // Resistance Mapping Workflows
            GBS_RES(read_pairs_ch)

            if (params.other_res_dbs != 'none'){
                OTHER_RES(read_pairs_ch)
                id_ch = GBS_RES.out.join(OTHER_RES.out)
            } else {
                id_ch = GBS_RES.out
            }

            // Once GBS or both resistance workflows are complete, trigger resistance typing
            res_typer(id_ch, tmp_dir, params.restyper_min_read_depth)

            // Combine serotype and resistance type results for each sample
            sero_res_ch = serotyping.out.join(res_typer.out)

            combine_results(sero_res_ch)

            // Combine samples and output results files
            combine_results.out.sero_res_incidence
                .collectFile(name: file("${results_dir}/${params.sero_res_incidence_out}"), keepHeader: true)

            combine_results.out.res_alleles_variants
                .collectFile(name: file("${results_dir}/${params.alleles_variants_out}"), keepHeader: true)

            combine_results.out.res_variants
                .collectFile(name: file("${results_dir}/${params.variants_out}"), keepHeader: true)

        }

        // MLST
        if (params.run_mlst){

            MLST(read_pairs_ch)
            MLST.out.subscribe { it ->
                it.copyTo(file("${results_dir}"))
            }

        }

        // Surface Typing Process
        if (params.run_surfacetyper){

            surface_typer(read_pairs_ch, file(params.gbs_surface_typer_db, checkIfExists: true),
                params.surfacetyper_min_read_depth, params.surfacetyper_min_coverage,
                params.surfacetyper_max_divergence)

            finalise_surface_typer_results(surface_typer.out)

            // Combine results for surface typing
            finalise_surface_typer_results.out.surface_protein_incidence
                .collectFile(name: file("${results_dir}/${params.surface_protein_incidence_out}"), keepHeader: true)
            finalise_surface_typer_results.out.surface_protein_variants
                .collectFile(name: file("${results_dir}/${params.surface_protein_variants_out}"), keepHeader: true)

        }

        // PBP Typer
        if (params.run_pbptyper){

            // Check if contigs specified
            if (params.contigs == ""){
                println("Please specify contigs with --contigs.")
                println("Print help with --contigs")
                System.exit(1)
            }

            contig_paths = Channel
                .fromPath(params.contigs, checkIfExists: true)
                .map { file -> tuple(file.baseName, file) }

            get_pbp_genes(contig_paths, file(params.gbs_blactam_db, checkIfExists: true), params.pbp_frac_align_threshold, params.pbp_frac_identity_threshold)

            // Get PBP existing and new alleles
            PBP1A(get_pbp_genes.out)
            PBP2B(get_pbp_genes.out)
            PBP2X(get_pbp_genes.out)

            PBP1A.out
                .collectFile(name: file("${results_dir}/${params.existing_pbp_alleles_out}"), keepHeader: true)
            PBP2B.out
                .collectFile(name: file("${results_dir}/${params.existing_pbp_alleles_out}"), keepHeader: true)
            PBP2X.out
                .collectFile(name: file("${results_dir}/${params.existing_pbp_alleles_out}"), keepHeader: true)
        }
}
