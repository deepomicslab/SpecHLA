#!/bin/env bash

# Run all the test scripts.

set -euo pipefail

run_test () {
    local test_script_path=$1

    test_dir="$(dirname $test_script_path)"
    script="$(basename $test_script_path)"

    (cd $test_dir && ./$script)
}

run_test exon/test_exon.sh
run_test hybrid/test_hybrid.sh
run_test pacbio/test_pacbio.sh
run_test whole/test_whole.sh

# Ran last, depends on test output for whole.
run_test whole/test_hla_loh.sh
