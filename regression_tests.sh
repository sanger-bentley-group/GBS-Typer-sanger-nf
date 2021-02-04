#!/usr/bin/env bash

in_dir=regression_test_data/input_data
out_dir=regression_test_data/output_data

SUDO_OPT=
while getopts ":s" opt; do
  case ${opt} in
    s ) SUDO_OPT="sudo "
      ;;
    \? ) echo "Usage: $0 [-s to run with sudo]"
         exit 1
      ;;
  esac
done

echo "Starting regression tests..."
echo ""

# Run pipeline on test input data
${SUDO_OPT}nextflow -log nextflow_test.log run main.nf --reads "$in_dir/test_{1,2}.fastq.gz" --results_dir "$out_dir" --run_surfacetyper --run_pbptyper --contigs "$in_dir/test.fa" --output 'test'
cat nextflow_test.log
echo ""

function file_diff {
    local test=${1}
    local reference=${2}

    diff ${test} ${reference} > /dev/null 2>&1
    status=$?

    # If files different then warning
    if [[ ${status} -gt 0 ]]; then
        if [[ ${status} -eq 1 ]]; then
            echo ""
            echo "The contents of ${test} is not expected."
        elif [[ ${status} -eq 2 ]]; then
            echo ""
            echo "Unable to perform differences check."
            if [[ ! -f $test ]]; then
                echo "Expected output ${test} is missing."
            fi
        fi
        return 1
    fi
}

error_status=0
# Check for test_drug_cat_alleles.txt output
file_diff "${out_dir}/test_drug_cat_alleles_variants.txt" "${out_dir}/reference_drug_cat_alleles_variants.txt"
out=$?
error_status=$(($error_status | $out))

# Check for test_gbs_res_variants.txt output
file_diff "${out_dir}/test_gbs_res_variants.txt" "${out_dir}/reference_gbs_res_variants.txt"
out=$?
error_status=$(($error_status | $out))

# Check for test_serotype_res_incidence.txt output
file_diff "${out_dir}/test_serotype_res_incidence.txt" "${out_dir}/reference_serotype_res_incidence.txt"
out=$?
error_status=$(($error_status | $out))

# Check for test_surface_protein_incidence.txt output
file_diff "${out_dir}/test_surface_protein_incidence.txt" "${out_dir}/reference_surface_protein_incidence.txt"
out=$?
error_status=$(($error_status | $out))

# Check for test_surface_protein_variants.txt output
file_diff "${out_dir}/test_surface_protein_variants.txt" "${out_dir}/reference_surface_protein_variants.txt"
out=$?
error_status=$(($error_status | $out))

# Check for test_existing_pbp_alleles.txt output
file_diff "${out_dir}/test_existing_pbp_alleles.txt" "${out_dir}/reference_existing_pbp_alleles.txt"
out=$?
error_status=$(($error_status | $out))

# Error if any output files missing or not expected
if [[ ${error_status} -eq 1 ]]; then
    echo ""
    echo "Test failed. Outputs listed may be missing or their contents not expected."
    exit 1
else
    echo "Test passed. All outputs expected."
fi
