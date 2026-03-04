"""Tests for the original_metadata feature in finalize_metadata_from_original."""
import gzip
import json
import os
import tempfile
import pytest
from viral_usher.viral_usher_build import (
    finalize_metadata_from_original,
    get_extra_metadata,
    make_taxonium_config,
    open_maybe_decompress,
)


@pytest.fixture
def workdir():
    """Create a temporary working directory."""
    old_cwd = os.getcwd()
    with tempfile.TemporaryDirectory() as tmpdir:
        os.chdir(tmpdir)
        yield tmpdir
        os.chdir(old_cwd)


@pytest.fixture
def original_metadata_file(workdir):
    """Create a sample original metadata.tsv.gz with tree sequence metadata."""
    path = os.path.join(workdir, "original_metadata.tsv.gz")
    with gzip.open(path, 'wt') as f:
        f.write("strain\tnum_date\taccession\tcountry\tdate\n")
        f.write("USA/seq1|ACC001|2024-01-15\t2024.041096\tACC001\tUSA\t2024-01-15\n")
        f.write("GBR/seq2|ACC002|2024-03-20\t2024.218579\tACC002\tGBR\t2024-03-20\n")
        f.write("FRA/seq3|ACC003|2024-06-01\t2024.415301\tACC003\tFRA\t2024-06-01\n")
    return path


@pytest.fixture
def sample_names_file(workdir):
    """Create a sample_names file listing sequences in the final tree."""
    path = os.path.join(workdir, "tree_samples.txt")
    with open(path, 'w') as f:
        # Two existing sequences + one new sequence from extra_fasta
        f.write("USA/seq1|ACC001|2024-01-15\n")
        f.write("GBR/seq2|ACC002|2024-03-20\n")
        f.write("new_sample_1\n")
    return path


@pytest.fixture
def extra_metadata_file(workdir):
    """Create a simple extra_metadata file for new sequences."""
    path = os.path.join(workdir, "extra_metadata.tsv")
    with open(path, 'w') as f:
        f.write("name\tcountry\tdate\n")
        f.write("new_sample_1\tJPN\t2024-07-01\n")
    return path


def test_basic_original_metadata(workdir, original_metadata_file, sample_names_file):
    """Test that original_metadata produces valid output with correct columns and rows."""
    metadata_out, rename_out, date_min, date_max, extra_added, extra_mapped = \
        finalize_metadata_from_original(
            original_metadata_file, sample_names_file,
            extra_fasta_names=None, extra_metadata='', extra_metadata_date_column='')

    assert os.path.exists(metadata_out)
    assert os.path.exists(rename_out)

    # Read the output metadata
    with gzip.open(metadata_out, 'rt') as f:
        lines = f.readlines()

    # Header + 2 existing sequences (seq3 not in sample_names)
    assert len(lines) == 3
    header = lines[0].strip().split('\t')
    assert header == ['strain', 'num_date', 'accession', 'country', 'date']

    # Check first row
    row1 = lines[1].strip().split('\t')
    assert row1[0] == 'USA/seq1|ACC001|2024-01-15'
    assert row1[2] == 'ACC001'
    assert row1[3] == 'USA'

    # Check date range
    assert date_min is not None
    assert date_max is not None
    assert date_min < date_max

    # No extra cols since original_metadata defines the schema
    assert extra_added == []
    assert extra_mapped == {}


def test_original_metadata_with_extra_fasta(workdir, original_metadata_file, sample_names_file):
    """Test that new sequences from extra_fasta get rows with empty fields."""
    extra_fasta_names = {'new_sample_1', 'new_sample_not_in_tree'}

    metadata_out, rename_out, date_min, date_max, _, _ = \
        finalize_metadata_from_original(
            original_metadata_file, sample_names_file,
            extra_fasta_names=extra_fasta_names, extra_metadata='',
            extra_metadata_date_column='')

    with gzip.open(metadata_out, 'rt') as f:
        lines = f.readlines()

    # Header + 2 existing + 1 new (new_sample_not_in_tree not in tree)
    assert len(lines) == 4
    new_row = lines[3].strip().split('\t')
    assert new_row[0] == 'new_sample_1'
    # All other fields should be empty (no extra_metadata provided)
    assert all(v == '' for v in new_row[1:])


def test_original_metadata_with_extra_fasta_and_metadata(
        workdir, original_metadata_file, sample_names_file, extra_metadata_file):
    """Test that extra_metadata values are mapped into original_metadata columns for new sequences."""
    extra_fasta_names = {'new_sample_1'}

    metadata_out, rename_out, date_min, date_max, _, _ = \
        finalize_metadata_from_original(
            original_metadata_file, sample_names_file,
            extra_fasta_names=extra_fasta_names, extra_metadata=extra_metadata_file,
            extra_metadata_date_column='date')

    with gzip.open(metadata_out, 'rt') as f:
        lines = f.readlines()

    header = lines[0].strip().split('\t')
    assert len(lines) == 4
    new_row = lines[3].strip().split('\t')
    assert new_row[0] == 'new_sample_1'

    # 'country' column should be mapped from extra_metadata
    country_idx = header.index('country')
    assert new_row[country_idx] == 'JPN'

    # 'date' column should be mapped from extra_metadata
    date_idx = header.index('date')
    assert new_row[date_idx] == '2024-07-01'

    # 'num_date' should be computed from the date
    num_date_idx = header.index('num_date')
    assert float(new_row[num_date_idx]) > 2024.0


def test_original_metadata_filters_by_sample_names(workdir, original_metadata_file, sample_names_file):
    """Test that only sequences in sample_names are included in output."""
    metadata_out, _, _, _, _, _ = \
        finalize_metadata_from_original(
            original_metadata_file, sample_names_file,
            extra_fasta_names=None, extra_metadata='', extra_metadata_date_column='')

    with gzip.open(metadata_out, 'rt') as f:
        lines = f.readlines()

    # seq3 (FRA) is in original_metadata but NOT in sample_names -> excluded
    strains = [line.strip().split('\t')[0] for line in lines[1:]]
    assert 'FRA/seq3|ACC003|2024-06-01' not in strains
    assert 'USA/seq1|ACC001|2024-01-15' in strains
    assert 'GBR/seq2|ACC002|2024-03-20' in strains


def test_original_metadata_no_duplicate_columns(workdir, original_metadata_file, sample_names_file):
    """Test that the output header has no duplicate columns (the original bug)."""
    metadata_out, _, _, _, _, _ = \
        finalize_metadata_from_original(
            original_metadata_file, sample_names_file,
            extra_fasta_names=None, extra_metadata='', extra_metadata_date_column='')

    with gzip.open(metadata_out, 'rt') as f:
        header = f.readline().strip().split('\t')

    # No duplicates
    assert len(header) == len(set(header)), f"Duplicate columns found: {header}"


def test_make_taxonium_config_with_original_metadata_cols(workdir):
    """Test that make_taxonium_config uses original_metadata_cols for color_by_options."""
    original_cols = ['strain', 'num_date', 'accession', 'country', 'date']
    config_path = make_taxonium_config(
        date_min=2024.0, date_max=2024.5,
        nextclade_clade_columns='', got_extra_fasta=True,
        no_genbank=True, extra_added_cols=[], extra_mapped_cols={},
        original_metadata_cols=original_cols)

    with open(config_path, 'r') as f:
        config = json.load(f)

    color_by = config['colorBy']['colorByOptions']
    # Should include meta_ prefixed columns from original_metadata (excluding 'strain')
    assert 'meta_num_date' in color_by
    assert 'meta_country' in color_by
    assert 'meta_date' in color_by
    assert 'meta_accession' in color_by
    # Should NOT include strain (it's the index column)
    assert 'meta_strain' not in color_by


def test_make_taxonium_config_without_original_metadata(workdir):
    """Test that make_taxonium_config still works without original_metadata_cols."""
    config_path = make_taxonium_config(
        date_min=2024.0, date_max=2024.5,
        nextclade_clade_columns='', got_extra_fasta=True,
        no_genbank=True, extra_added_cols=[], extra_mapped_cols={})

    with open(config_path, 'r') as f:
        config = json.load(f)

    # Should still produce valid config
    assert 'colorBy' in config
    assert 'colorByOptions' in config['colorBy']
