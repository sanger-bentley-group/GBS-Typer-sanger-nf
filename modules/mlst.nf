process srst2_for_mlst {

    input:
    tuple val(pair_id), file(reads) // ID and paired read files
    file(mlst_alleles) // Path to MLST alleles file
    file(mlst_definitions) // Path to MLST definition file
    val(min_coverage) // String of minimum coverage parameter(s) for SRST2

    //publishDir "./${tmp_dir}/${pair_id}", mode: 'move', overwrite: true, pattern: "${pair_id}_${db_name}_*__fullgenes__*__results.txt"

    output:
    tuple val(pair_id), file("${pair_id}*.bam"), file("${pair_id}__mlst__${mlst_name}__results.txt"), emit: bam_and_srst2_results
    tuple val(pair_id), file("${pair_id}__mlst__${mlst_name}__results.txt"), emit: srst2_results

    script:
    mlst_name=mlst_alleles.getSimpleName()

    """
    set +e
    srst2 --samtools_args '\\-A' --mlst_delimiter '_' --input_pe ${reads[0]} ${reads[1]} --output ${pair_id} --save_scores --mlst_db ${mlst_alleles} --mlst_definitions ${mlst_definitions} --min_coverage ${min_coverage}
    touch ${pair_id}__mlst__${mlst_name}__results.txt
    """
}

process get_mlst_allele_and_pileup {

    input:
    tuple val(pair_id), file(bam_file), file(results_file)
    val(min_read_depth)
    file(mlst_alleles)

    output:
    path("${pair_id}_new_mlst_alleles.fasta"), emit: new_alleles, optional: true
    path("${pair_id}_new_mlst_pileup.txt"), emit: pileup, optional: true
    path("${pair_id}_existing_mlst_alleles.txt"), emit: existing_alleles, optional: true
    path("${pair_id}_new_mlst_alleles.log"), emit: new_alleles_status

    """
    # Get alleles from mismatches in SRST2 MLST results file
    samtools index ${bam_file}
    get_alleles_from_srst2_mlst.py --mlst_results_file ${results_file} --min_read_depth ${min_read_depth} --output_prefix ${pair_id}
    if [[ -f ${pair_id}_new_mlst_alleles.txt ]]
    then
        num_alleles=\$(cat ${pair_id}_new_mlst_alleles.txt | wc -l)
    else
        num_alleles=0
    fi

    # Get consensus allele and variant pileup for each allele
    if [[ \${num_alleles} -gt 1 ]]
    then
        echo "${pair_id}: New MLST alleles found." > ${pair_id}_new_mlst_alleles.log
        for ((i=0;i<\${num_alleles};i++)); # Skip first line of alleles file
        do
            target_allele=\$(sed -n \"\${i+1}p\" ${pair_id}_mlst_alleles.txt)

            mlst_bam=\"\${target_allele}_${bam_file}\"
            mlst_vcf=\"\$(basename \${mlst_bam} .bam).vcf\"
            samtools view -b ${bam_file} \${target_allele} > \${mlst_bam}
            samtools index \${mlst_bam} \$(basename \${mlst_bam} .bam).bai
            echo \${target_allele} > \${target_allele}.txt
            get_targets_from_db.py -f ${mlst_alleles} -t \${target_allele}.txt -o CHECK_MLST_
            freebayes -q 20 -p 1 -f CHECK_MLST_\${target_allele}_ref.fna \${mlst_bam} -v \${mlst_vcf}
            bgzip \${mlst_vcf}
            tabix -p vcf \${mlst_vcf}.gz

            cat CHECK_MLST_\${target_allele}_ref.fna | vcf-consensus \${mlst_vcf}.gz >> ${pair_id}_new_mlst_alleles.fasta
            echo '\nNew MLST Allele Pileup:' >> ${pair_id}_new_mlst_pileup.txt
            samtools mpileup -f ${mlst_alleles} ${bam_file} -r \${target_allele} >> ${pair_id}_new_mlst_pileup.txt
        done
    elif [[ \${num_alleles} -eq 1 ]]
    then
        cat ${pair_id}_new_mlst_alleles.txt > ${pair_id}_new_mlst_alleles.log
    else
        echo "${pair_id}: No new MLST alleles found." > ${pair_id}_new_mlst_alleles.log
    fi

    """

}
