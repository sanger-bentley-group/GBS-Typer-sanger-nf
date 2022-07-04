process combine_results {
    input:
    // ID, serotyping results, resistance incidence, resistance alleles, resistance variants, surface protein incidence, surface protein variants, MLST allelic frequency
    tuple val(pair_id), file(sero_results), file(res_incidence), file(res_alleles), file(res_variants), file(surface_protein_incidence), file(surface_protein_variants), file(mlst_allelic_frequency)
    path config

    output:
    path("${pair_id}_id_combined_output.txt") optional true

    """
    combine_results.py combine_all \
        -i ${pair_id} \
        -t ${config} \
        -s "${sero_results}" \
        -r "${res_incidence}" \
        -v "${res_variants}" \
        -m "${mlst_allelic_frequency}" \
        -x "${surface_protein_incidence}" \
        -o ${pair_id}
    """
}

process finalise_sero_res_results {

    input:
    // ID, serotyping results, resistance incidence, resistance alleles and resistance variants results,
    // surface protein incidence, surface protein variants
    tuple val(pair_id), file(sero_results), file(res_incidence), file(res_alleles), file(res_variants)
    path config

    output:
    path("${pair_id}_sero_res_incidence.txt"), emit: sero_res_incidence
    path("${pair_id}_id_alleles_variants.txt"), emit: res_alleles_variants
    path("${pair_id}_id_variants.txt"), emit: res_variants

    """
    combine_results.py sero_res \
        -i ${pair_id} \
        -t ${config} \
        -s "${sero_results}" \
        -r "${res_incidence}" \
        -a "${res_alleles}" \
        -v "${res_variants}" \
        -o ${pair_id}
    """
}

process finalise_surface_typer_results {

    input:
    // ID, surface protein incidence, surface protein variants
    tuple val(pair_id), file(surface_protein_incidence_in), file(surface_protein_variants_in)
    path config

    output:
    path("${pair_id}_surface_protein_incidence.txt"), emit: surface_protein_incidence
    path("${pair_id}_surface_protein_variants.txt"), emit: surface_protein_variants

    """
    combine_results.py surface_typer \
        -i ${pair_id} \
        -t ${config} \
        --surface_incidence_results "${surface_protein_incidence_in}" \
        --surface_variants_results "${surface_protein_variants_in}" \
        -o ${pair_id}
    """
}

process finalise_pbp_existing_allele_results {

    input:
    // ID, existing PBP allele file
    tuple val(pair_id), file(pbp_existing_allele)
    path config

    output:
    path("${pair_id}_existing_PBP_allele.txt")

    """
    combine_results.py pbp_typer \
        -i ${pair_id} \
        -t ${config} \
        --pbp_existing_allele_results "${pbp_existing_allele}" \
        -o ${pair_id}
    """
}

process include_qc {
    input:
    path qc_summary
    path qc_complete
    path in_silico
    path config

    output:
    path("${output}")

    script:
    output=qc_and_in_silico.txt
    """
    include_qc.py -h ${config} \
        -s ${qc_summary} \
        -c ${qc_complete} \
        -i ${in_silico} \
        -o ${output}
    """
}
