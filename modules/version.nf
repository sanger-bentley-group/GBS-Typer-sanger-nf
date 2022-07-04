process add_version {
    input:
    path tab_file

    output:
    path("${output}")

    script:
    version=params.version
    output="output.txt"
    """
    echo "version" > ${output}

    if [ -z ${version} ];
    then
        use_version=\$(echo \$(git -C $baseDir describe --tags))
    else
        use_version=\$(echo ${version})
    fi

    head -1 ${tab_file} | sed -e 's/\$/\tversion/' > ${output}
    sed '1d' ${tab_file} | sed -e "s/\$/\t\${use_version}/" >> ${output}
    """
}
