# Quick test of simple build from scratch, with segment 4 of epizootic hemorrhagic disease virus
# which is small and doesn't have very many sequences
set -beEu -o pipefail

workdir=$1

testdir=$(realpath $(dirname "${BASH_SOURCE[0]}"))

mkdir -p $workdir
workdir=$(realpath $workdir)

time viral_usher init \
    -r NC_013399.1 \
    -s "epizootic hemorrhagic disease virus" \
    -t 3431252 \
    -w $workdir \
    -c $workdir/config.toml

CONTAINER_IMAGE=${APPTAINER_CONTAINER:-angiehinrichs/viral_usher:development}
time viral_usher build \
    -d $CONTAINER_IMAGE \
    -c $workdir/config.toml

# Now run build on the same config, but with --update
time viral_usher build \
    -d $CONTAINER_IMAGE \
    -c $workdir/config.toml \
    -u

# Check outputs exist
cd $workdir
source $testdir/check_output_files.sh

