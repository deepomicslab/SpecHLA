#!/bin/env bash

set -euo pipefail

output_dir="./output"
whole_typing_test_output_dir="../whole/output"
whole_typing_sample_name="HG00733"

freq_list_file="${output_dir}/freq.list"
whole_typing_sample_dir="${whole_typing_test_output_dir}/${whole_typing_sample_name}"
whole_typing_genotype_file="${whole_typing_sample_dir}/hla.result.txt"

if test ! -d $whole_typing_test_output_dir; then
    # Run the test and generate the output.
    (cd $(dirname $whole_typing_test_output_dir) && \
        . $(basname $whole_typing_test_output_dir))
fi

mkdir -p "$output_dir"

ls "${whole_typing_sample_dir}/"*_freq.txt > "$freq_list_file"

# Call HLA LOH using ficticious purity & ploidy.
perl ../../script/cal.hla.copy.pl \
    -purity 0.9 \
    -ploidy 2 \
    -S test \
    -F "$freq_list_file" \
    -T "$whole_typing_genotype_file" \
    -O "$output_dir"
