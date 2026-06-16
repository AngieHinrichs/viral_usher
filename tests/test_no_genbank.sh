# Test a simple build with URL for ref_fasta, abs paths for ref_gbff and extra_fasta; run --no_genbank
set -beEu -o pipefail

workdir=$1

testdir=$(realpath $(dirname "${BASH_SOURCE[0]}"))

mkdir -p $workdir
workdir=$(realpath $workdir)

# Create test config with absolute paths
config=$(mktemp)
cat > $config <<EOF
ref_fasta = 'https://raw.githubusercontent.com/AngieHinrichs/viral_usher/refs/heads/main/tests/test_data/ref.fa'
ref_gbff = '$(pwd)/tests/test_data/ref.gbff'
extra_fasta = '$(pwd)/tests/test_data/sequences.fa'
workdir = '$workdir'
EOF

# Run build with test data
CONTAINER_IMAGE=${APPTAINER_CONTAINER:-angiehinrichs/viral_usher:development}
viral_usher build --config $config --docker_image $CONTAINER_IMAGE --no_genbank

# Check outputs exist
cd $workdir
source $testdir/check_output_files.sh

rm $config
