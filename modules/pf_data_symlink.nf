process pf_data_symlink {

    input:
    path lanes_file

    output:
    path "${data_dir}/*/*_{1,2}.fastq.gz"

    script:
    pf_version=params.pf_version
    data_dir="data_dir"

    """
    #!/bin/bash

    module load pf/${pf_version}
    pf data -t file -i ${lanes_file} -l ${data_dir}
    """
}
