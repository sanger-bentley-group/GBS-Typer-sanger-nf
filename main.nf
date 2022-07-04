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
include {finalise_sero_res_results; finalise_surface_typer_results; finalise_pbp_existing_allele_results; combine_results; include_qc} from './modules/combine.nf'
include {add_version} from './modules/version.nf'
include {pf_data_symlink} from './modules/pf_data_symlink.nf'
include {get_lane_list} from './modules/get_lane_list.nf'

// Import GBS QC modules
include {get_file_destinations} from './GBS_QC_nf/modules/get_file_destinations.nf'
include {relative_abundance} from './GBS_QC_nf/modules/relative_abundance.nf'
include {number_of_contigs} from './GBS_QC_nf/modules/number_of_contigs.nf'
include {collate_qc_data} from './GBS_QC_nf/modules/collate_qc_data.nf'
include {contig_gc_content} from './GBS_QC_nf/modules/contig_gc_content.nf'
include {genome_length} from './GBS_QC_nf/modules/genome_length.nf'
include {get_qc_stats_from_pf} from './GBS_QC_nf/modules/get_qc_stats_from_pf.nf'
include {depth_of_coverage} from './GBS_QC_nf/modules/depth_of_coverage.nf'
include {breadth_of_coverage} from './GBS_QC_nf/modules/breadth_of_coverage.nf'
include {get_proportion_HET_SNPs} from './GBS_QC_nf/modules/get_proportion_HET_SNPs.nf'
include {HET_SNPs} from './GBS_QC_nf/modules/HET_SNPs.nf'


// Help message
if (params.help){
    printHelp()
    exit 0
}

// Check if QC specified
if (params.run_qc){
    if (params.lanes == ""){
        println("Please specify a file with a list of lane ids.")
        System.exit(1)
    }
    // Create lanes channel
    lanes_ch = Channel.fromPath( params.lanes, checkIfExists: true )
} else {
    if (params.reads == ""){
        println("Please specify reads with --reads.")
        System.exit(1)
    }
    // Create read pairs channel
    Channel.fromFilePairs( params.reads, checkIfExists: true )
        .set { read_pairs_ch }
}

// Check if output specified
if (params.output == ""){
    println("Please specify and output prefix with --output.")
    System.exit(1)
}

if (!params.run_sero_res && !params.run_surfacetyper && !params.run_mlst && !params.run_pbptyper){
    println("Please specify one or more pipelines to run.")
    System.exit(1)
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

// Create results directory if it doesn't already exist
results_dir = file(params.results_dir)
results_dir.mkdir()

// Output files
params.sero_res_incidence_out = "${params.output}_serotype_res_incidence.txt"
params.variants_out =  "${params.output}_gbs_res_variants.txt"
params.alleles_variants_out = "${params.output}_drug_cat_alleles_variants.txt"
params.res_accessions_out = "${params.output}_resfinder_accessions.txt"
params.existing_pbp_alleles_out = "${params.output}_existing_pbp_alleles.txt"
params.surface_protein_incidence_out = "${params.output}_surface_protein_incidence.txt"
params.surface_protein_variants_out = "${params.output}_surface_protein_variants.txt"
params.existing_mlst_alleles_out = "${params.output}_existing_sequence_types.txt"
params.new_mlst_alleles_status = "${params.output}_new_mlst_alleles.log"
params.gbs_typer_report_no_qc = "${params.output}_gbs_typer_report_no_qc.txt"
params.gbs_typer_report = "${params.output}_gbs_typer_report.txt"

// Workflow for reads QC
workflow reads_qc {
    take:
    file_dest_ch
    headers_ch
    lanes_ch

    main:
    relative_abundance(file_dest_ch, headers_ch, lanes_ch)
    qc_report = relative_abundance.out

    emit:
    qc_report
}

// Workflow for assemblies QC
workflow assemblies_qc {
    take:
    file_dest_ch
    qc_stats_ch
    het_stats_ch
    headers_ch
    lanes_ch

    main:
    number_of_contigs(file_dest_ch, headers_ch, lanes_ch)
    contig_gc_content(file_dest_ch, headers_ch, lanes_ch)
    genome_length(file_dest_ch, headers_ch, lanes_ch)
    depth_of_coverage(qc_stats_ch, headers_ch, lanes_ch)
    breadth_of_coverage(qc_stats_ch, headers_ch, lanes_ch)
    HET_SNPs(het_stats_ch, headers_ch, lanes_ch)

    number_of_contigs.out
    .combine(contig_gc_content.out)
    .combine(genome_length.out)
    .combine(depth_of_coverage.out)
    .combine(breadth_of_coverage.out)
    .combine(HET_SNPs.out)
    .set { qc_report }

    emit:
    qc_report
}

// Resistance mapping with the GBS resistance database
workflow GBS_RES {

    take:
        reads

    main:

        gbs_res_typer_db = file(params.gbs_res_typer_db, checkIfExists: true)
        gbs_res_targets_db = file(params.gbs_res_targets_db, checkIfExists: true)

        // Split GBS target sequences from GBS resistance database into separate FASTA files per sequence
        split_target_RES_sequences(gbs_res_typer_db, gbs_res_targets_db)

        // Map genomes to GBS resistance database using SRST2
        srst2_for_res_typing(reads, gbs_res_typer_db, params.gbs_res_min_coverage, params.gbs_res_max_divergence)
        fullgenes = srst2_for_res_typing.out.fullgenes

        // Split sam file for each GBS target sequence
        split_target_RES_seq_from_sam_file(srst2_for_res_typing.out.bam_files, gbs_res_targets_db)

        // Get consensus sequence using freebayes
        freebayes(split_target_RES_seq_from_sam_file.out, split_target_RES_sequences.out)
        consensus = freebayes.out.consensus

    emit:
        fullgenes
        consensus
}

// Resistance mapping with the other resistance databases
workflow OTHER_RES {

    take:
        reads

    main:
        other_res_db = file(params.other_res_db, checkIfExists: true)
        // Map genomes to resistance database using SRST2
        srst2_for_res_typing(reads, other_res_db, params.other_res_min_coverage, params.other_res_max_divergence)
        fullgenes = srst2_for_res_typing.out.fullgenes

    emit:
        fullgenes
}

// MLST pipeline
workflow MLST {

    take:
        reads

    main:
        // Run SRST2 MLST
        srst2_for_mlst(reads, params.mlst_min_coverage)

        // Get new consensus allele and pileup data
        get_mlst_allele_and_pileup(srst2_for_mlst.out.bam_and_srst2_results, params.mlst_min_read_depth)

        // Collect outputs
        new_alleles = get_mlst_allele_and_pileup.out.new_alleles
        pileup = get_mlst_allele_and_pileup.out.pileup
        existing_alleles = get_mlst_allele_and_pileup.out.existing_alleles
        status = get_mlst_allele_and_pileup.out.new_alleles_status
        srst2_results = srst2_for_mlst.out.srst2_results

    emit:
        new_alleles
        pileup
        existing_alleles
        status
        srst2_results
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
        finalise_pbp_existing_allele_results(get_pbp_alleles.out.existing_pbp, file(params.config, checkIfExists: true))

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
        finalise_pbp_existing_allele_results(get_pbp_alleles.out.existing_pbp, file(params.config, checkIfExists: true))

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
        finalise_pbp_existing_allele_results(get_pbp_alleles.out.existing_pbp, file(params.config, checkIfExists: true))

    emit:
        // Emit existing PBP alleles for collection
        finalise_pbp_existing_allele_results.out
}

// Main Workflow
workflow {

    main:

        if (params.run_qc){

            // Set project directory
            projectDir=params.qc_project_dir

            // Run assembly QC
            get_file_destinations(lanes_ch)
            get_qc_stats_from_pf(lanes_ch)
            get_proportion_HET_SNPs(lanes_ch)

            headers_ch = Channel.fromPath( params.config, checkIfExists: true )

            assemblies_qc(get_file_destinations.out, get_qc_stats_from_pf.out, get_proportion_HET_SNPs.out, headers_ch, lanes_ch)

            // Run reads QC
            reads_qc(get_file_destinations.out, headers_ch, lanes_ch)

            // Collate QC reports
            collate_qc_data(reads_qc.out.qc_report, assemblies_qc.out.qc_report)

            // Get list of lane ids that did not FAIL
            get_lane_list(collate_qc_data.out.complete)
            not_failed_lanes_ch=get_lane_list.out

            // Create a read pairs ch
            pf_data_symlink(not_failed_lanes_ch)
            Channel.fromFilePairs( pf_data_symlink.out, checkIfExists: true )
                .set { read_pairs_ch }
        }

        if (params.run_sero_res){

            // Serotyping Process
            serotyping(read_pairs_ch, params.serotyper_min_read_depth)

            // Resistance Mapping Workflows
            GBS_RES(read_pairs_ch)
            OTHER_RES(read_pairs_ch)

            // Once GBS or both resistance workflows are complete, trigger resistance typing
            GBS_RES.out.fullgenes
            .join(GBS_RES.out.consensus)
            .join(OTHER_RES.out.fullgenes)
            .set { res_files_ch }

            res_typer(res_files_ch, params.restyper_min_read_depth, file(params.config, checkIfExists: true))

            // Combine serotype and resistance type results for each sample
            sero_res_ch = serotyping.out.join(res_typer.out.res_out)

            finalise_sero_res_results(sero_res_ch, file(params.config, checkIfExists: true))

            // Combine samples and output results files
            finalise_sero_res_results.out.sero_res_incidence
                .collectFile(name: file("${results_dir}/${params.sero_res_incidence_out}"), keepHeader: true)

            finalise_sero_res_results.out.res_alleles_variants
                .collectFile(name: file("${results_dir}/${params.alleles_variants_out}"), keepHeader: true)

            finalise_sero_res_results.out.res_variants
                .collectFile(name: file("${results_dir}/${params.variants_out}"), keepHeader: true)

            res_typer.out.res_accessions
                .collectFile(name: file("${results_dir}/${params.res_accessions_out}"))
        }

        // MLST
        if (params.run_mlst){

            MLST(read_pairs_ch)
            MLST.out.new_alleles.subscribe { it ->
                it.copyTo(file("${results_dir}"))
            }
            MLST.out.pileup.subscribe { it ->
                it.copyTo(file("${results_dir}"))
            }
            MLST.out.existing_alleles
                .collectFile(name: file("${results_dir}/${params.existing_mlst_alleles_out}"), keepHeader: true, sort: true)
            MLST.out.status
                .collectFile(name: file("${results_dir}/${params.new_mlst_alleles_status}"), keepHeader: false, sort: true)

        }

        // Surface Typing Process
        if (params.run_surfacetyper){

            surface_typer(read_pairs_ch, file(params.gbs_surface_typer_db, checkIfExists: true),
                params.surfacetyper_min_read_depth, params.surfacetyper_min_coverage,
                params.surfacetyper_max_divergence)

            finalise_surface_typer_results(surface_typer.out, file(params.config, checkIfExists: true))

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
            .concat(PBP2B.out, PBP2X.out)
            .set { PBP_all }

            PBP_all
                .collectFile(name: file("${results_dir}/${params.existing_pbp_alleles_out}"), keepHeader: true, sort: true)
        }

        // Combine serotype, resistance, allelic profile, surface typer and GBS resistance variants
        if (params.run_sero_res & params.run_surfacetyper & params.run_mlst){

            // Combine serotype and resistance type results for each sample
            combined_ch = serotyping.out
                .join(res_typer.out.res_out)
                .join(surface_typer.out)
                .join(MLST.out.srst2_results)

            combine_results(combined_ch, file(params.config, checkIfExists: true))

            combine_results.out
                .collectFile(name: file("${results_dir}/${params.gbs_typer_report_no_qc}"), keepHeader: true, sort: true)
                .set { in_silico_ch }

            // Include QC in final output
            if (params.run_qc) {
                include_qc(collate_qc_data.out.summary, collate_qc_data.out.complete, in_silico_ch, headers_ch)
                report_ch = include_qc.out
            } else {
                report_ch = in_silico_ch
            }

            // Get version of pipeline
            add_version(report_ch)

            add_version.out
            .subscribe { it ->
                it.copyTo(file("${results_dir}/${params.gbs_typer_report}"))
            }
        }
}
