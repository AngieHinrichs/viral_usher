import argparse
import importlib.metadata
import sys
from .init import handle_init
from .build import handle_build
from . import config


def main():
    viral_usher_version = importlib.metadata.version('viral_usher')
    parser = argparse.ArgumentParser(
        prog="viral_usher",
        description=f"viral_usher {viral_usher_version} — Build viral phylogenies from NCBI data.",
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument("-v", "--version", action="store_true", help="display viral_usher version and exit")
    subparsers = parser.add_subparsers(dest="command")

    # Init subcommand
    init_parser = subparsers.add_parser("init", help="Initialize configuration for a new species/taxon and reference")
    init_parser.add_argument("-r", "--refseq", help="RefSeq accession to use as reference/root, if known")
    init_parser.add_argument("-t", "--taxonomy_id", help="NCBI Taxonomy ID of the viral species, if known")
    init_parser.add_argument("-s", "--species", help="Viral species name (if Taxonomy ID is not known)")
    init_parser.add_argument("-f", "--fasta", help="Additional sequences to include in tree")
    init_parser.add_argument("-m", "--metadata", help="Metadata file (TSV): first column must match names in --fasta file)")
    init_parser.add_argument("--metadata_date_column", help="Column in metadata file containing sequence collection date (default: look for 'date')")
    init_parser.add_argument("-x", "--nextclade_dataset", help="Nextclade dataset path (e.g. 'nextstrain/dengue/all') (default: auto-detect)")
    init_parser.add_argument("--annotate_allele_frequency", type=float, help="matUtils annotate --allele_frequency parameter (only used with --nextclade_dataset; default: no annotation)")
    init_parser.add_argument("--annotate_mask_frequency", type=float, help="matUtils annotate --mask_frequency parameter (only used with --nextclade_dataset; default: no annotation)")
    init_parser.add_argument("-l", "--min_length_proportion", type=float, help=f"Minimum proportion of RefSeq length to require for GenBank sequences (default: {config.DEFAULT_MIN_LENGTH_PROPORTION})")
    init_parser.add_argument("-n", "--max_N_proportion", type=float, help=f"Maximum proportion of N bases to allow in GenBank sequences (default: {config.DEFAULT_MAX_N_PROPORTION})")
    init_parser.add_argument("-a", "--max_parsimony", type=int, help=f"Remove sequences from the tree with parsimony score (private substitution count) greater than this (default: {config.DEFAULT_MAX_PARSIMONY})")
    init_parser.add_argument("-b", "--max_branch_length", type=int, help=f"Remove branches from the tree with branch length (substitution count) greater than this (default: {config.DEFAULT_MAX_BRANCH_LENGTH})")
    init_parser.add_argument(      "--title", type=str, help="Title for the tree display")
    init_parser.add_argument("-w", "--workdir", help="Directory in which tree files will be built")
    init_parser.add_argument("-c", "--config", type=str, help="Path to config file output")
    init_parser.add_argument("--ref_fasta", type=str, help="Path to local reference genome FASTA file (if not using RefSeq)")
    init_parser.add_argument("--ref_gbff", type=str, help="Path to local reference genome GenBank flat file (if not using RefSeq)")
    init_parser.add_argument("--use_viral_usher_trees", action="store_true", help="If the viral_usher_trees repository has a pre-built tree is available for the selected RefSeq, use it as a starting tree instead of building the tree from scratch")
    init_parser.add_argument("--taxonium_overlay_html", type=str, help="Path to HTML file to use as Taxonium overlay (instead of auto-generated file)")

    # Build subcommand
    build_parser = subparsers.add_parser("build", help="Run the pipeline to download sequences and build a tree")
    build_parser.add_argument("-c", "--config", type=str, required=True, help="Path to config file input")
    build_parser.add_argument("-d", "--docker_image", type=str, help=f"Use this docker image instead of docker.io/{config.DEFAULT_DOCKER_IMAGE}")
    build_parser.add_argument("-n", "--no_genbank", action="store_true", help="Do not download sequences from GenBank; build tree using only sequences in extra_fasta")
    build_parser.add_argument("-u", "--update", action="store_true", help="Update an existing tree with new sequences (requires optimized.pb.gz tree from previous build)")
    build_parser.add_argument("--threads", type=int, default=0, help="Number of threads to use for alignment and tree building steps (default: 0, use all available cores)")

    args = parser.parse_args()

    # Parse known args to check for -h or -v
    args, unknown = parser.parse_known_args()

    # Check for -v/--version
    if args.version:
        print(f"viral_usher {viral_usher_version}")
        sys.exit(0)

    # If -h/--help is present, let argparse print help and exit
    if "-h" in sys.argv or "--help" in sys.argv:
        print(f"unknown is {unknown}\n\n", file=sys.stderr)
        parser.parse_args()

    # A subcommand is required unless -h or -v were given
    if not args.command:
        parser.error("A subcommand {init,build} is required unless -h/--help or -v/--version is given.")

    # Parse all args for subcommand handlers
    args = parser.parse_args()

    if args.command == "init":
        handle_init(args)
    elif args.command == "build":
        handle_build(args)
    else:
        parser.print_help()
        exit(1)


if __name__ == "__main__":
    main()
