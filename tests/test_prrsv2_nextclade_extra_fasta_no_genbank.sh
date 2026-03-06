# Build from scratch: PRRSV2 with nextclade annotations, using viral_usher_trees, with extra fasta
set -beEu -o pipefail

workdir=$1

testdir=$(realpath $(dirname "${BASH_SOURCE[0]}"))

mkdir -p $workdir
workdir=$(realpath $workdir)

time viral_usher init \
    -r NC_038291.1 \
    -s "PRRSV2" \
    -t 2499685 \
    -x community/isuvdl/mazeller/prrsv2/orf5/yimim2023 \
    --use_viral_usher_trees \
    -f tests/test_data/prrsv2_extra.fa \
    -w $workdir \
    -c $workdir/config.toml

# Just one update --no_genbank run so only extra fasta seqs are added
time viral_usher build \
    -d angiehinrichs/viral_usher:development \
    -c $workdir/config.toml \
    -u \
    --no_genbank

# Check outputs exist
cd $workdir
source $testdir/check_output_files.sh

