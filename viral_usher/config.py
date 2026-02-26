import os
import shutil
import sys
import tempfile
import tomllib
import urllib.request


DEFAULT_MIN_LENGTH_PROPORTION = 0.8
DEFAULT_MAX_N_PROPORTION = 0.25
DEFAULT_MAX_PARSIMONY = 1000
DEFAULT_MAX_BRANCH_LENGTH = 10000
DEFAULT_DOCKER_IMAGE = "angiehinrichs/viral_usher"


def is_url(path):
    return path.startswith('http://') or path.startswith('https://')


def handle_path_or_url(config_value):
    """
    Resolve a config value that can be either a local file path or a URL.
    If it's a URL (starts with http:// or https://), download it to a temporary directory.
    Return the path to use (either the original path or the downloaded file path).
    """
    if not config_value:
        return config_value

    # Convert Path objects to strings
    config_value = str(config_value)

    if is_url(config_value):
        print(f"Downloading {config_value}...")
        filename = os.path.basename(config_value.split('?')[0])  # Remove query params if any
        if not filename:
            filename = "downloaded_file"

        temp_dir = tempfile.mkdtemp()
        temp_path = os.path.join(temp_dir, filename)

        try:
            with urllib.request.urlopen(config_value) as response:
                with open(temp_path, 'wb') as temp_file:
                    shutil.copyfileobj(response, temp_file)
            print(f"Downloaded to {temp_path}")
            return temp_path
        except Exception as e:
            shutil.rmtree(temp_dir)
            print(f"Error downloading {config_value}: {e}", file=sys.stderr)
            sys.exit(1)
    else:
        # It's a local file path, return as is
        return config_value


def read_config(config_path):
    resolved_path = handle_path_or_url(config_path)
    try:
        with open(resolved_path, 'rb') as f:
            config = tomllib.load(f)
    except FileNotFoundError:
        print(f"Config file '{config_path}' does not exist.", file=sys.stderr)
        raise
    return config


def check_refseq_vs_ref(config, check_paths=True):
    """Ensure that the config specifies either RefSeq info or local reference files, but not both."""
    refseq_acc = config.get('refseq_acc', '')
    refseq_assembly = config.get('refseq_assembly', '')
    ref_fasta = config.get('ref_fasta', '')
    ref_gbff = config.get('ref_gbff', '')
    if (refseq_acc or refseq_assembly) and (ref_fasta or ref_gbff):
        raise ValueError("Config must specify either refseq_acc and refseq_assembly, or ref_fasta and ref_gbff, but not both.")
    if (not refseq_acc or not refseq_assembly) and (not ref_fasta or not ref_gbff):
        raise ValueError("Config must specify either refseq_acc and refseq_assembly, or ref_fasta and ref_gbff.")
    if ref_fasta:
        if check_paths and not is_url(ref_fasta) and (not os.path.exists(ref_fasta) or not os.access(ref_fasta, os.R_OK)):
            raise ValueError(f"Reference FASTA file '{ref_fasta}' does not exist or is not readable.")
    if ref_gbff:
        if check_paths and not is_url(ref_gbff) and (not os.path.exists(ref_gbff) or not os.access(ref_gbff, os.R_OK)):
            raise ValueError(f"Reference GenBank file '{ref_gbff}' does not exist or is not readable.")


def check_refseq_acc(refseq_acc):
    if not refseq_acc or not isinstance(refseq_acc, str):
        raise ValueError("RefSeq accession must be a non-empty string.")
    if not (refseq_acc.startswith("NC_") or refseq_acc.startswith("NZ_") or refseq_acc.startswith("AC_")):
        raise ValueError("RefSeq accession must start with 'NC_', 'NZ_' or 'AC_'.")


def check_taxonomy_id(taxonomy_id):
    if not taxonomy_id or not isinstance(taxonomy_id, str):
        raise ValueError("Taxonomy ID must be a non-empty string.")
    if not taxonomy_id.isdigit():
        raise ValueError("Taxonomy ID must be a numeric string.")


def check_refseq_assembly(refseq_assembly):
    if not refseq_assembly or not isinstance(refseq_assembly, str):
        raise ValueError("RefSeq assembly must be a non-empty string.")
    if not refseq_assembly.startswith("GCA_") and not refseq_assembly.startswith("GCF_"):
        raise ValueError("RefSeq assembly must start with 'GCA_' or 'GCF_'.")


def check_config(config, check_paths=True, no_genbank=False):
    if not no_genbank:
        check_taxonomy_id(config['taxonomy_id'])
    elif config.get('taxonomy_id'):
        check_taxonomy_id(config['taxonomy_id'])
    check_refseq_vs_ref(config, check_paths)
    if config.get('refseq_acc', ''):
        check_refseq_acc(config['refseq_acc'])
        check_refseq_assembly(config['refseq_assembly'])
    extra_fasta = config.get('extra_fasta', '')
    if extra_fasta and check_paths and not (os.path.exists(extra_fasta) and os.access(extra_fasta, os.R_OK)):
        raise ValueError(f"Fasta file {extra_fasta} does not exist or is not readable.")


def parse_config(config_path, resolve_url_keys=True, no_genbank=False):
    """Read the config file and validate its contents.  Download URLs unless resolve_url_keys is False."""
    config = read_config(config_path)

    if resolve_url_keys:
        # Resolve paths that can be URLs
        url_resolvable_keys = ['ref_fasta', 'ref_gbff', 'extra_fasta', 'extra_metadata', 'update_tree_input', 'taxonium_overlay_html']
        for key in url_resolvable_keys:
            if key in config:
                config[key] = handle_path_or_url(config[key])

    check_config(config, check_paths=resolve_url_keys, no_genbank=no_genbank)
    return config


def write_config(config, config_path, check_paths=True, no_genbank=False):
    """Write a config TOML file to the specified path."""
    check_config(config, check_paths, no_genbank)
    with open(config_path, 'w') as f:
        if config.get('refseq_acc'):
            ref_desc = f"RefSeq {config['refseq_acc']}"
        else:
            ref_desc = f"reference genome {config['ref_fasta']}"
        taxid_desc = f", taxonomy ID {config['taxonomy_id']}" if config.get('taxonomy_id') else ""
        print(f"# viral_usher config for {ref_desc}{taxid_desc}\n", file=f)
        for key in config:
            print(f"{key} = '{config[key]}'", file=f)
