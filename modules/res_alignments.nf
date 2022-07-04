process srst2_for_res_typing {

    container 'bluemoon222/gbs-typer-sanger-nf:0.0.7'

    input:
    tuple val(pair_id), file(reads) // ID and paired read files
    path db // File of resistance database file(s)
    val(min_coverage) // String of minimum coverage parameter(s) for SRST2
    val(max_divergence) // String of maximum coverage parameter(s) for SRST2

    output:
    tuple val(pair_id), file("${pair_id}*.bam"), emit: bam_files
    tuple val(pair_id), file("${pair_id}__fullgenes__${db_name}__results.txt"), emit: fullgenes

    script:
    db_name=db.getSimpleName()
    """
    set +e
    srst2 --samtools_args '\\-A' --input_pe ${reads[0]} ${reads[1]} --output ${pair_id} --log --save_scores --min_coverage ${min_coverage} --max_divergence ${max_divergence} --gene_db ${db}
    touch ${pair_id}__fullgenes__${db_name}__results.txt
    """
}

process split_target_RES_sequences {

    container 'bluemoon222/gbs-typer-sanger-nf:0.0.7'

    input:
    file(fasta_file) // FASTA file of GBS target sequences
    file(targets_file) // Text file of GBS targets of interest

    output:
    file("CHECK_*")

    """
    get_targets_from_db.py -f ${fasta_file} -t ${targets_file} -o CHECK_
    """
}

process split_target_RES_seq_from_sam_file {

    container 'bluemoon222/gbs-typer-sanger-nf:0.0.7'

    input:
    tuple val(pair_id), file(bam_file) // ID and corresponding BAM file from mapping
    file(targets_file) // Text file of GBS targets of interest

    output:
    val(pair_id)
    file("*_*_${pair_id}*.bam")
    file("*_*_${pair_id}*.bai")

    """
    set +e
    samtools view -h ${bam_file} > \$(basename ${bam_file} .bam).sam
    get_targets_from_samfile.py -s \$(basename ${bam_file} .bam).sam -t ${targets_file} -i ${pair_id} -o CHECK_
    for check_sam_file in CHECK_*_${pair_id}*.sam; do
        samtools view -bS \${check_sam_file} > \$(basename \${check_sam_file} .sam).bam
        samtools index \$(basename \${check_sam_file} .sam).bam \$(basename \${check_sam_file} .sam).bai
    done
    touch dummy_dummy_${pair_id}_dummy.bam
    touch dummy_dummy_${pair_id}_dummy.bai
    """
}

process freebayes {

    container 'bluemoon222/gbs-typer-sanger-nf:0.0.7'

    input:
    val(pair_id) // ID
    file(target_bam) // BAM file from a mapped target sequence of interest
    file(target_bai) // Corresponding BAM index file
    file(target_ref) // FASTA file of target sequence

    output:
    tuple val(pair_id), file("${pair_id}_consensus_seq.fna"), emit: consensus

    """
    set +e
    for check_bam_file in CHECK_*_${pair_id}*.bam; do
        target=\$(echo \${check_bam_file} | sed 's/CHECK_//g' | sed 's/_${pair_id}.*//g')
        freebayes -q 20 -p 1 -f CHECK_\${target}_ref.fna \${check_bam_file} -v CHECK_\${target}_${pair_id}_seq.vcf
        bgzip CHECK_\${target}_${pair_id}_seq.vcf
        tabix -p vcf CHECK_\${target}_${pair_id}_seq.vcf.gz
        cat CHECK_\${target}_ref.fna | vcf-consensus CHECK_\${target}_${pair_id}_seq.vcf.gz >> ${pair_id}_consensus_seq.fna
    done
    touch ${pair_id}_consensus_seq.fna
    """
}
