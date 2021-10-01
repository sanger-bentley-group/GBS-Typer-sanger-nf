process res_typer {

    input:
    tuple val(pair_id), path(gbs_fullgenes), path(gbs_consensus), path(other_fullgenes)
    val(min_read_depth) // Minimum read depth threshold

    output:
    tuple val(pair_id), file(inc_output_file), file(alleles_output_file), file(variants_output_file)

    script:
    inc_output_file="${pair_id}_res_incidence.txt"
    alleles_output_file="${pair_id}_res_alleles_variants.txt"
    variants_output_file="${pair_id}_res_gbs_variants.txt"
    """
    set +e
    process_res_typer_results.py --srst2_gbs_fullgenes ${gbs_fullgenes} --srst2_gbs_consensus ${gbs_consensus} --srst2_other_fullgenes ${other_fullgenes} --min_read_depth ${min_read_depth} --output_prefix ${pair_id}

    touch ${inc_output_file}
    touch ${alleles_output_file}
    touch ${variants_output_file}
    """
}
