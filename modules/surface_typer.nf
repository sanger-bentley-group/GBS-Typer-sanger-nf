process surface_typer {

    input:
    tuple val(pair_id), file(reads)
    file(surface_protein_db)
    val(min_read_depth) // Minimum read depth threshold
    val(min_coverage) // Minimum coverage threshold
    val(max_divergence) // Maximum allowed divergence threshold

    output:
    tuple val(pair_id), file(inc_output_file), file(variants_output_file)

    script:
    inc_output_file="${pair_id}_surface_protein_incidence_sample.txt"
    variants_output_file="${pair_id}_surface_protein_variants_sample.txt"
    """
    set +e
    srst2 --samtools_args '\\-A' --input_pe ${reads[0]} ${reads[1]} --output ${pair_id}_SURFACE --log --save_scores --min_coverage ${min_coverage} --max_divergence ${max_divergence} --gene_db ${surface_protein_db}
    process_surface_typer_results.py --srst2_gbs_fullgenes ${pair_id}_SURFACE --surface_db ${surface_protein_db} --output_prefix ${pair_id} --min_read_depth ${min_read_depth}

    touch ${inc_output_file}
    touch ${variants_output_file}
    """
}
