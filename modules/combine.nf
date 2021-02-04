process combine_results {

    input:
    // ID, serotyping results, resistance incidence, resistance alleles and resistance variants results,
    // surface protein incidence, surface protein variants
    tuple val(pair_id), file(sero_results), file(res_incidence), file(res_alleles), file(res_variants)

    output:
    path("${pair_id}_sero_res_incidence.txt"), emit: sero_res_incidence
    path("${pair_id}_id_alleles_variants.txt"), emit: res_alleles_variants
    path("${pair_id}_id_variants.txt"), emit: res_variants

    """
    combine_results.py sero_res\
        -i ${pair_id} -s "${sero_results}" -r "${res_incidence}" -a "${res_alleles}" -v "${res_variants}" -o ${pair_id}
    """
}

process finalise_surface_typer_results {

    input:
    // ID, surface protein incidence, surface protein variants
    tuple val(pair_id), file(surface_protein_incidence_in), file(surface_protein_variants_in)

    output:
    path("${pair_id}_surface_protein_incidence.txt"), emit: surface_protein_incidence
    path("${pair_id}_surface_protein_variants.txt"), emit: surface_protein_variants

    """
    combine_results.py surface_typer \
        -i ${pair_id} --surface_incidence_results "${surface_protein_incidence_in}" \
        --surface_variants_results "${surface_protein_variants_in}" \
        -o ${pair_id}
    """
}

process finalise_pbp_existing_allele_results {

    input:
    // ID, existing PBP allele file
    tuple val(pair_id), file(pbp_existing_allele)

    output:
    path("${pair_id}_existing_PBP_allele.txt")

    """
    combine_results.py pbp_typer \
        -i ${pair_id} --pbp_existing_allele_results "${pbp_existing_allele}" \
        -o ${pair_id}
    """
}
