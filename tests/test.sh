#!/bin/bash

set -beEu -o pipefail

cd ~/github/viral_usher
# Build a local docker image
time docker build --platform linux/amd64 -t angiehinrichs/viral_usher:development .

function run_test {
    workdir=$(mktemp -d -t viral_usher_test.XXXXXX)
    echo ========================================= $1 =========================================
    bash $1 $workdir
    rm -rf $workdir
    echo ""
}

# Theo's quick file & URL no GenBank test:
run_test tests/test_no_genbank.sh

# Quick test of simple build from scratch, with segment 4 of epizootic hemorrhagic disease virus
# which is small and doesn't have very many sequences:
run_test tests/test_ehdv_4.sh

# Same EHDV segment 4, but use viral_usher_trees as a base
run_test tests/test_ehdv_4_viral_usher_trees.sh

# Now build from scratch with extra_fasta (no metadata)
run_test tests/test_ehdv_4_extra_fasta.sh

# Build from scratch with extra_fasta and extra_metadata
run_test tests/test_ehdv_4_extra_fasta_metadata.sh

# Use viral_usher_trees as a base, add extra_fasta and extra_metadata
run_test tests/test_ehdv_4_extra_fasta_metadata_viral_usher_trees.sh

# Build from scratch: PRRSV2 with nextclade annotations, no extra fasta
run_test tests/test_prrsv2_nextclade.sh

# Build from viral_usher_trees: PRRSV2 with nextclade annotations and extra fasta
run_test tests/test_prrsv2_nextclade_extra_fasta.sh

# Build from viral_usher_trees: PRRSV2 with nextclade annotations and extra fasta, no genbank
run_test tests/test_prrsv2_nextclade_extra_fasta_no_genbank.sh

echo Success
