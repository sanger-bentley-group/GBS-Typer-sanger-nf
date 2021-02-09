process target2mic {

    input:
    tuple val(pair_id), file(sero_results), file(res_incidence), file(res_alleles), file(res_variants)


    output:
    tuple val(pair_id), file("${pair_id}_mic_predictions.txt")

    """
    get_target2MIC.py --res_file "${res_alleles}" --output_prefix "${pair_id}"
    """
}