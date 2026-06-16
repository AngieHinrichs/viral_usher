# EHDV segment 4 from scratch with extra_fasta but no extra_metadata
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
    -w $workdir \
    -c $workdir/config.toml

CONTAINER_IMAGE=${APPTAINER_CONTAINER:-angiehinrichs/viral_usher:development}
time viral_usher build \
    -d $CONTAINER_IMAGE \
    -c $workdir/config.toml

# Again, with --update
time viral_usher build \
    -d $CONTAINER_IMAGE \
    -c $workdir/config.toml \
    -u

# Check outputs exist
cd $workdir
source $testdir/check_output_files.sh

