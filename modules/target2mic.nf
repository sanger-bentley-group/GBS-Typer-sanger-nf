process target2mic {

    input:
    file(res_alleles)
    file(pbp_alleles)
    file(output_file)

    output:
    file(output_file)

    """
    get_target2MIC.py --res_file "${res_alleles}" --pbp_file "${pbp_alleles}" --output "${output_file}"
    """
}
