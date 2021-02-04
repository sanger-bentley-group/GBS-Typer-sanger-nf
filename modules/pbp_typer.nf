process get_pbp_genes {
    input:
    tuple val(pair_id), file(contigs)
    file(blactam_ref)
    val(frac_align_len_threshold)
    val(frac_identity_len_threshold)

    output:
    tuple val(pair_id), file("${pair_id}_*bed"), file(contigs), optional: true

    """
    # Build a blast reference database from the assmeblies
    makeblastdb -in ${contigs} -dbtype nucl -out ${pair_id}_contig_blast_db

    # Blast the blactam database against the blast reference database
    blastn -db ${pair_id}_contig_blast_db -query ${blactam_ref} -outfmt 6 -word_size 7 -out ${pair_id}_blast_blactam.out

    # Get BED file of PBP fragments
    get_pbp_genes_from_contigs.py --blast_out_file ${pair_id}_blast_blactam.out --query_fasta ${blactam_ref} --frac_align_len_threshold ${frac_align_len_threshold} --frac_identity_threshold ${frac_identity_len_threshold} --output_prefix ${pair_id}_
    """
}

process get_pbp_alleles {
    input:
    tuple val(pair_id), file(bed_file), file(contigs)
    val(pbp_type)
    file(gbs_blactam_db)

    output:
    path "${pair_id}_${pbp_type}_PBP_new_allele.faa", optional: true, emit: new_pbp
    tuple val(pair_id), file("${pair_id}_${pbp_type}_PBP_existing_allele.txt"), optional: true, emit: existing_pbp

    """
    # Get PBP alleles
    if [ -f ${pair_id}_${pbp_type}.bed ]; then
        # Convert BED file to fasta
        bedtools getfasta -s -fi ${contigs} -bed ${pair_id}_${pbp_type}.bed -fo ${pair_id}_${pbp_type}_contig_fragments.fasta

        # Translate fasta file into separate faa files using the translator in python
        translate_pbp_genes.py --blactam_fasta ${pair_id}_${pbp_type}_contig_fragments.fasta --output_file ${pair_id}_${pbp_type}.faa

        # Build blast databases of PBP alleles and blast amino acids against database
        makeblastdb -in ${gbs_blactam_db} -dbtype prot -out ${gbs_blactam_db}
        blastp -db ${gbs_blactam_db} -query ${pair_id}_${pbp_type}.faa -outfmt 6 -out ${pair_id}_blast_${pbp_type}.out

        # Get identical or imperfect hits
        get_pbp_alleles.py --blast_out_file ${pair_id}_blast_${pbp_type}.out --query_fasta ${pair_id}_${pbp_type}.faa --output_prefix ${pair_id}_${pbp_type}_PBP
    fi
    """
}
