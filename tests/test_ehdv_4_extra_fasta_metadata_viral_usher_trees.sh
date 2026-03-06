# EHDV segment 4 with extra fasta and metadata, using viral_usher_trees as a base
set -beEu -o pipefail

workdir=$1

testdir=$(realpath $(dirname "${BASH_SOURCE[0]}"))

mkdir -p $workdir
workdir=$(realpath $workdir)

time viral_usher init \
    -r NC_013399.1 \
    -s "epizootic hemorrhagic disease virus" \
    -t 3431252 \
    -f tests/test_data/ehdv_extra.fa \
    -m tests/test_data/ehdv_extra.metadata.tsv \
    -w $workdir \
    --use_viral_usher_trees \
    -c $workdir/config.toml

# I think --update is forced by --use_viral_usher_trees but use it anyway
time viral_usher build \
    -d angiehinrichs/viral_usher:development \
    -c $workdir/config.toml \
    -u

# TODO: make sure the generated metadata has the expected additional columns and
# metadata rows for user-added sequences.

# Check outputs exist
cd $workdir
source $testdir/check_output_files.sh

