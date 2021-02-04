process res_typer {

    input:
    val(pair_id) // ID
    path(res_dir) // Directory of resistance mapping results
    val(min_read_depth) // Minimum read depth threshold

    output:
    tuple val(pair_id), file("${pair_id}_res_incidence.txt"), file("${pair_id}_res_alleles_variants.txt"), file("${pair_id}_res_gbs_variants.txt")

    """
    process_res_typer_results.py --srst2_gbs_fullgenes ${res_dir}/${pair_id}/${pair_id}_GBS_RES_*__fullgenes__*__results.txt --srst2_gbs_consensus ${res_dir}/${pair_id}/${pair_id}_consensus_seq.fna --srst2_other_fullgenes ${res_dir}/${pair_id}/${pair_id}_OTHER_RES_*__fullgenes__*__results.txt --min_read_depth ${min_read_depth} --output_prefix ${pair_id}
    rm -f ${res_dir}/${pair_id}/${pair_id}*
    """
}
