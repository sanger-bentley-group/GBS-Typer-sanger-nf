process pf_data_symlink {

    input:
    path lanes_file

    output:
    file "${read_files}"

    script:
    pf_version=params.pf_version
    read_files="read_files.txt"

    """
    module load pf/${pf_version}
    pf data -t file -i ${lanes_file} > dir_list.txt
    num=\$(cat dir_list.txt | wc -l)
    for ((i=1;i<=\${num};i++))
    do
        dir=\$(sed -n "\${i}p" dir_list.txt)
        lane=\$(echo \${dir} | rev | cut -d/ -f1 | rev)
        echo -e "\${lane}\t\${dir}/\${lane}_1.fastq.gz" >> ${read_files}
        echo -e "\${lane}\t\${dir}/\${lane}_2.fastq.gz" >> ${read_files}
    done
    """
}
