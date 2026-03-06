# 'init' subcommand: get user's preferences and write config file

import sys
import os
import importlib.metadata
import re
import logging
import requests
import tomllib
from . import ncbi_helper
from . import nextclade_helper
from . import config

viral_usher_trees_raw_base = "https://github.com/AngieHinrichs/viral_usher_trees/raw/refs/heads/main/"


def get_input(prompt):
    """Get input from the user.  Quit gracefully if it looks like the user wants out."""
    try:
        text = input(prompt)
    except (KeyboardInterrupt, EOFError):
        print("")
        sys.exit(1)
    if re.match(r"^[Qq](uit)?", text):
        print("Exiting...")
        sys.exit(0)
    return text


def prompt_int_choice(prompt, min_value=1, max_value=None):
    """Get an integer choice from the user, with optional min and max values.
    min_value will be presented as default in case of empty input."""
    while True:
        try:
            choice = get_input(f"{prompt} [{min_value}]: ")
            choice = int(choice) if choice else min_value
            if min_value <= choice and (max_value is None or choice <= max_value):
                return choice
            else:
                print(f"Please enter a number between {min_value} and {max_value}.")
        except ValueError:
            print("Invalid input. Please enter a number.")


def prompt_taxonomy_id(species_term, tax_entries):
    """Given a list of NCBI Datasets taxonomy match records, display to the user and ask them to choose one,
    also giving an option for None of the above / go back.  Return 1-based index of user's choice."""
    print(f"\nFound the following matches for '{species_term}':")
    for idx, entry in enumerate(tax_entries):
        print(f"{idx+1}. {entry['sci_name']}")
    print(f"{len(tax_entries)+1}. None of the above, go back")
    return prompt_int_choice("Enter a number to choose an option", 1, len(tax_entries)+1)


def prompt_refseq_id(taxid, refseq_entries):
    """Given a list of NCBI RefSeq records, display to the user and ask them to choose one,
    also giving an option for None of the above / go back.  Return 1-based index of user's choice."""
    print(f"\nFound the following RefSeq IDs for Taxonomy ID {taxid}:")
    for idx, entry in enumerate(refseq_entries):
        label = f"{entry['accession']}: {entry['title']}"
        if (entry["strain"] and entry["strain"] != "No strain" and entry["strain"] not in label):
            label += f" ({entry['strain']})"
        print(f"{idx+1}. {label}")
    print(f"{len(refseq_entries)+1}. None of the above, go back")
    return prompt_int_choice("Enter a number to choose an option", 1, len(refseq_entries)+1)


def get_species_taxonomy_refseq(ncbi):
    """Prompt user for a species search term.  Look up Taxonomy IDs for species,
    prompt the user to select a Taxonomy ID, look up RefSeqs for the selected Taxonomy ID,
    prompt the user to select a RefSeq.  Return Taxonomy ID, RefSeq accession and RefSeq Assembly ID."""
    species_term = ""
    while not species_term:
        species_term = get_input("\nWhat is your virus of interest? ")
    return lookup_taxonomy_refseq(ncbi, species_term)


def lookup_taxonomy_refseq(ncbi, species_term):
    """Look up Taxonomy IDs for species, prompt the user to select a Taxonomy ID,
    look up RefSeqs for the selected Taxonomy ID, prompt the user to select a RefSeq.
    Return Taxonomy ID, RefSeq accession and RefSeq Assembly ID."""
    print(f"Looking up NCBI Taxonomy entries for '{species_term}'...")
    tax_entries = ncbi.get_taxonomy_entries('"' + species_term + '"')
    if tax_entries:
        return get_taxonomy_refseq(ncbi, species_term, tax_entries)
    else:
        print(f"\nNo matches found for '{species_term}'.")
        return get_species_taxonomy_refseq(ncbi)


def get_taxonomy_refseq(ncbi, species_term, tax_entries):
    """Prompt the user to select a Taxonomy ID, look up RefSeqs for the selected Taxonomy ID,
    prompt the user to select a RefSeq.  Return Taxonomy ID, RefSeq accession and RefSeq Assembly ID."""
    choice = prompt_taxonomy_id(species_term, tax_entries)
    if choice > len(tax_entries):
        # User says none of the above, go back.
        return get_species_taxonomy_refseq(ncbi)
    else:
        taxid = tax_entries[choice-1]["tax_id"]
        return lookup_refseq(ncbi, species_term, tax_entries, taxid)


def lookup_refseq(ncbi, species_term, tax_entries, taxid):
    """Look up RefSeqs for the selected Taxonomy ID, prompt the user to select a RefSeq.
    Return Taxonomy ID and RefSeq ID."""
    print(f"Looking up RefSeqs associated with Taxonomy ID {taxid}...")
    refseq_entries = ncbi.get_refseqs_for_taxid(taxid)
    if refseq_entries:
        return get_refseq(ncbi, species_term, tax_entries, taxid, refseq_entries)
    else:
        print(f"No RefSeqs found for Taxonomy ID '{taxid}'.  Try a different match from species search.")
        return get_taxonomy_refseq(ncbi, species_term, tax_entries)


def get_refseq(ncbi, species_term, tax_entries, taxid, refseq_entries):
    """Prompt the user to select a RefSeq.  Return Taxonomy ID, RefSeq accession and RefSeq Assembly ID."""
    choice = prompt_refseq_id(taxid, refseq_entries)
    if choice > len(refseq_entries):
        # User says none of the above, go back.
        return get_taxonomy_refseq(ncbi, species_term, tax_entries)
    else:
        refseq_id = refseq_entries[choice-1]["accession"]
        print(f"Looking up NCBI Assembly accession for selected RefSeq '{refseq_id}'...")
        assembly_id = ncbi.get_assembly_acc_for_refseq_acc(refseq_id)
        if not assembly_id:
            print(f"Could not find assembly ID for RefSeq ID {refseq_id} -- can't download RefSeq.")
            print("Try choosing a different RefSeq.")
            return get_refseq(ncbi, species_term, tax_entries, taxid, refseq_entries)
        return species_term, taxid, refseq_id, assembly_id


def clade_columns_from_dataset(nextclade_dataset):
    """Return a comma-sep string of clade column names"""
    return ','.join(nextclade_dataset.get("clades", {}).keys())


def prompt_nextclade_path_columns(nextclade_datasets, species, is_interactive):
    """Prompt the user to select a Nextclade dataset from a list of paths."""
    if is_interactive:
        print(f"\nFound the following Nextclade datasets for '{species}':")
        for idx, dataset in enumerate(nextclade_datasets):
            print(f"{idx+1}. {dataset['path']}: {dataset['name']}")
        print(f"{len(nextclade_datasets)+1}. None of the above")
        choice = prompt_int_choice("Enter a number to choose an option", 1, len(nextclade_datasets)+1)
        if choice > len(nextclade_datasets):
            return "", ""
        else:
            dataset = nextclade_datasets[choice-1]
            return dataset["path"], clade_columns_from_dataset(dataset)
    else:
        print(f"Multiple Nextclade datasets found for '{species}':")
        for dataset in nextclade_datasets:
            print(f" - {dataset['path']}: {dataset['name']}")
        print("Please specify one of these using the -x/--nextclade_dataset option.")
        sys.exit(1)


def search_nextclade_datasets(datasets, term):
    """Search for Nextclade datasets that match the given term."""
    matches = []
    term_lower = term.lower()
    for dataset in datasets:
        if term_lower in dataset["name"].lower() or term_lower in dataset["path"].lower():
            matches.append(dataset)
    return matches


def get_nextclade_path_columns(args_nextclade, species, is_interactive):
    """Return a list of Nextclade dataset paths whose names match arg or species."""
    print(f"Looking up Nextclade datasets for '{species}'...")
    datasets = nextclade_helper.nextclade_get_index()
    if args_nextclade == "":
        return "", ""
    elif args_nextclade:
        for dataset in datasets:
            if dataset["path"] == args_nextclade:
                print(f"Using Nextclade dataset '{dataset['path']}'")
                return dataset['path'], clade_columns_from_dataset(dataset)
        print(f"Nextclade dataset '{args_nextclade}' not found. Available datasets:")
        for dataset in datasets:
            print(f"- {dataset['path']}: {dataset['name']}")
        print("Please try again with a different value for -x/--nextclade_dataset.")
        sys.exit(1)
    elif species:
        matches = search_nextclade_datasets(datasets, species)
        if len(matches) == 0 and ' ' in species:
            # Search each word in species separately; quit if a word has matches
            for word in species.lower().split(' '):
                if word in ["human", "virus", "fever", "genotype"] or len(word) < 3:
                    continue
                word_matches = search_nextclade_datasets(datasets, word)
                if word_matches:
                    matches = word_matches
                    break
        if len(matches) == 0:
            print(f"No Nextclade datasets found for '{species}'.")
            return "", ""
        elif len(matches) == 1:
            print(f"Found Nextclade dataset for '{species}': {matches[0]['path']}")
            return matches[0]["path"], clade_columns_from_dataset(matches[0])
        else:
            matches.reverse()
            return prompt_nextclade_path_columns(matches, species, is_interactive)
    else:
        return "", ""


def check_proportion(proportion):
    """Make sure the given value, while a string value, converts to a float between 0 and 1.  Empty/None is okay."""
    if proportion:
        try:
            value = float(proportion)
        except ValueError:
            return False, f"Sorry, need a value between 0 and 1, '{proportion}' doesn't parse"
        if value < 0:
            return False, "Sorry, need a value between 0 and 1, not a negative number"
        if value > 1.0:
            return False, f"Sorry, need a value between 0 and 1, {proportion} is too large"
    return True, None


def check_int(val, min=None, max=None):
    """Make sure the given value converts to an integer.  Empty/None is okay."""
    if val:
        try:
            value = int(val)
        except ValueError:
            return False, f"Sorry, need an integer, '{val}' doesn't parse"
        if min is not None:
            if value < min:
                return False, f"Sorry, need at least {min}, {val} is too low."
        if max is not None:
            if value > max:
                return False, f"Sorry, the maximum is {max}, {val} is too large."
    return True, None


def check_nonnegative_int(val):
    return check_int(val, min=0)


def check_optional_file_readable(filename):
    """If filename is empty, that's fine; if non-empty, make sure it's a readable file"""
    if filename:
        if not os.path.exists(filename):
            return False, f"Sorry, file '{filename}' does not exist"
        if not os.path.isfile(filename):
            return False, f"Sorry, file '{filename}' does not appear to be a regular file"
        else:
            try:
                f = open(filename, "r")
                f.close()
            except IOError:
                return False, f"Sorry, unable to read file '{filename}'."
    return True, None


def check_dir_exists_or_creatable(dirname):
    """If dirname exists, all good; otherwise try to create it.  Return False and error message if unable."""
    if not os.path.exists(dirname):
        print(f"Creating directory {dirname}...")
        try:
            os.makedirs(dirname)
        except OSError:
            return False, f"Sorry, can't create directory path '{dirname}', please enter a different one."
    elif not os.path.isdir(dirname):
        return False, f"Sorry, '{dirname}' exists but is not a directory, please enter a different name."
    return True, None


def prompt_with_checker(prompt, default, check_function):
    """Prompt the user for a value and check their answer.  If there's a problem then explain
    and prompt again until their answer meets criteria."""
    ok = False
    answer = None
    while not ok:
        if default is not None:
            answer = get_input(f"\n{prompt} [{default}]: ")
            if not answer:
                answer = default
        else:
            answer = get_input(f"\n{prompt}: ")
        ok, error_message = check_function(answer)
        if not ok:
            print(error_message)
    return answer


def check_write_config(config_contents, config_path):
    """Return True if able to write config to config_path, False otherwise."""
    try:
        config.write_config(config_contents, config_path)
        return True
    except (NotADirectoryError, OSError):
        return False


def search_viral_usher_trees(refseq_acc):
    """Search the viral_usher_trees repo for a tree with refseq_acc as its reference"""
    url = viral_usher_trees_raw_base + "tree_metadata.tsv"
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        lines = response.text.strip().split('\n')
        for line in lines[1:]:
            parts = line.split('\t')
            if parts[0] == refseq_acc:
                return parts[1]
    except requests.RequestException as e:
        print(f"Warning: Unable to access viral_usher_trees repository: {e}", file=sys.stderr)
    return None


def get_viral_usher_trees_name(refseq_acc, args_use_viral_usher_trees, is_interactive):
    """Check if there's a pre-built tree for this RefSeq in the viral_usher_trees repo.
    If so, ask the user if they want to use it.  Return the tree name if yes, '' otherwise."""
    if not refseq_acc:
        return ''
    tree_name = search_viral_usher_trees(refseq_acc)
    if tree_name:
        if args_use_viral_usher_trees:
            print(f"\nUsing pre-built tree for RefSeq {refseq_acc} from viral_usher_trees repository:\n{tree_name}")
            return tree_name
        elif is_interactive:
            print(f"\nA pre-built tree using RefSeq {refseq_acc} is available at:\n{tree_name}")
            choice = get_input("Would you like to use this tree instead of building a new one? (y/n) [y]: ") or "y"
            if choice.lower().startswith('y'):
                print(f"Okay, using {tree_name} as a starting tree.")
                return tree_name
            else:
                print("Okay, will build tree from scratch.")
                return ''
    return ''


def get_viral_usher_trees_refs(tree_name):
    """Download the config.toml file for the given tree and if it includes ref_fasta and ref_gbff entries then return those, otherwise return empty strings."""
    url = viral_usher_trees_raw_base + f"trees/{tree_name}/config.toml"
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        config_contents = tomllib.loads(response.text)
        ref_fasta = config_contents.get("ref_fasta", "").replace("./", viral_usher_trees_raw_base)
        ref_gbff = config_contents.get("ref_gbff", "").replace("./", viral_usher_trees_raw_base)
        return ref_fasta, ref_gbff
    except requests.RequestException as e:
        print(f"Warning: Unable to access viral_usher_trees repository to get config for {tree_name}: {e}", file=sys.stderr)
        sys.exit(1)


def get_min_length_proportion(args_min_length_proportion, is_interactive):
    min_length_proportion = config.DEFAULT_MIN_LENGTH_PROPORTION
    if args_min_length_proportion:
        min_length_proportion = args_min_length_proportion
        ok, error_message = check_proportion(min_length_proportion)
        if not ok:
            print(f"{error_message}\nPlease try again with a different value for --min_length_proportion", file=sys.stderr)
            sys.exit(1)
    elif is_interactive:
        min_length_proportion = prompt_with_checker("GenBank sequences will be filtered to retain only those whose length is at least this proportion of the RefSeq length",
                                                    min_length_proportion, check_proportion)
    return str(min_length_proportion)


def get_max_N_proportion(args_max_N_proportion, is_interactive):
    max_N_proportion = config.DEFAULT_MAX_N_PROPORTION
    if args_max_N_proportion:
        max_N_proportion = args_max_N_proportion
        ok, error_message = check_proportion(max_N_proportion)
        if not ok:
            print(f"{error_message}\nPlease try again with a different value for --max_N_proportion", file=sys.stderr)
            sys.exit(1)
    elif is_interactive:
        max_N_proportion = prompt_with_checker("GenBank sequences will be filtered to retain only those with at most this proportion of 'N' or ambiguous bases",
                                               max_N_proportion, check_proportion)
    return str(max_N_proportion)


def get_max_parsimony(args_max_parsimony, is_interactive):
    max_parsimony = config.DEFAULT_MAX_PARSIMONY
    if args_max_parsimony:
        max_parsimony = args_max_parsimony
        ok, error_message = check_nonnegative_int(max_parsimony)
        if not ok:
            print(f"{error_message}\nPlease try again with a different value for --max_parsimony", file=sys.stderr)
            sys.exit(1)
    elif is_interactive:
        max_parsimony = prompt_with_checker("After the tree is constructed, sequences with more than this number of private substitutions will be removed",
                                            max_parsimony, check_nonnegative_int)
    return str(max_parsimony)


def get_max_branch_length(args_max_branch_length, is_interactive):
    max_branch_length = config.DEFAULT_MAX_BRANCH_LENGTH
    if args_max_branch_length:
        max_branch_length = args_max_branch_length
        ok, error_message = check_nonnegative_int(max_branch_length)
        if not ok:
            print(f"{error_message}\nPlease try again with a different value for --max_branch_length", file=sys.stderr)
            sys.exit(1)
    elif is_interactive:
        max_branch_length = prompt_with_checker("After the tree is constructed, branches with more than this number of substitutions defining the branch will be removed",
                                                max_branch_length, check_nonnegative_int)
    return str(max_branch_length)


def get_user_fasta(args_fasta, is_interactive):
    fasta = ''
    if args_fasta is not None:
        fasta = os.path.abspath(args_fasta)
        ok, error_message = check_optional_file_readable(fasta)
        if not ok:
            print(f"{error_message}\nPlease try again with a different file for --fasta.", file=sys.stderr)
            sys.exit(1)
    elif is_interactive:
        fasta = prompt_with_checker("If you have your own fasta file, then enter its path", "", check_optional_file_readable)
    return fasta


def get_user_metadata(args_metadata, args_metadata_date_column, is_interactive):
    metadata_path = ''
    metadata_date_column = ''
    if args_metadata:
        metadata_path = os.path.abspath(args_metadata)
        ok, error_message = check_optional_file_readable(metadata_path)
        if not ok:
            print(f"{error_message}\nPlease try again with a different file for --metadata.", file=sys.stderr)
            sys.exit(1)
        if args_metadata_date_column:
            metadata_date_column = args_metadata_date_column
    elif args_metadata_date_column:
        print("Error: --metadata_date_column cannot be used without --metadata", file=sys.stderr)
        sys.exit(1)
    elif is_interactive:
        metadata_path = prompt_with_checker("If you have your own metadata file (TSV) whose first column matches names in the FASTA file, then enter its path", "", check_optional_file_readable)
        metadata_date_column = prompt_with_checker("If your metadata file has a column with sequence collection dates, enter its name (otherwise leave blank to look for a column named 'date')", "", lambda x: (True, None))
    return metadata_path, metadata_date_column


def get_workdir(args_workdir, is_interactive):
    if args_workdir:
        workdir = args_workdir
        ok, error_message = check_dir_exists_or_creatable(workdir)
        if not ok:
            print(f"{error_message}\nPlease try again with a different path for --workdir.", file=sys.stderr)
            sys.exit(1)
    else:
        workdir = prompt_with_checker("Enter directory where sequences should be downloaded and trees built", ".", check_dir_exists_or_creatable)
    return workdir


def make_config(config_contents, workdir, refseq_id, taxid, args_config, is_interactive):
    config_path_default = f"{workdir}/viral_usher_config_{refseq_id}_{taxid}.toml" if refseq_id else f"{workdir}/viral_usher_config_{taxid}.toml"
    if args_config:
        config_path = args_config
        if not check_write_config(config_contents, config_path):
            print(f"Unable to write to --config {config_path}.  Please try again with a different --config path.")
            sys.exit(1)
    elif is_interactive:
        config_ok = False
        while not config_ok:
            config_path = get_input(f"\nEnter path for config file [{config_path_default}]: ") or config_path_default
            config_ok = check_write_config(config_contents, config_path)
            if not config_ok:
                print(f"Unable to write config to {config_path}.")
            else:
                print(f"Wrote config file to {config_path}")
    else:
        config_path = config_path_default
        if not check_write_config(config_contents, config_path):
            print(f"Unable to write to default config path {config_path}.  Please try again with a different --config path.")
            sys.exit(1)
    return config_path


def handle_init(args):
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    ncbi = ncbi_helper.NcbiHelper()
    # First sort out the interdependent refseq and taxonomy id options.
    ref_fasta = ''
    ref_gbff = ''
    if args.refseq:
        if args.ref_fasta or args.ref_gbff:
            print("Error: --refseq cannot be used together with --ref_fasta or --ref_gbff", file=sys.stderr)
            sys.exit(1)
        refseq_id = args.refseq
        assembly_id = ncbi.get_assembly_acc_for_refseq_acc(refseq_id)
        if not assembly_id:
            print(f"Could not find assembly ID for RefSeq ID {refseq_id} -- can't download RefSeq.")
            sys.exit(1)
        species = args.species
        taxid = args.taxonomy_id
        if not species or not taxid:
            species, taxid = ncbi.get_species_taxid_for_refseq(refseq_id)
        if not species and args.species:
            species = args.species
        if not taxid and args.taxonomy_id:
            taxid = args.taxonomy_id
        if not species or not taxid:
            print(f"Could not find species or taxonomy ID for RefSeq ID {refseq_id}.", file=sys.stderr)
            sys.exit(1)
        taxid = ncbi.get_species_level_taxid(taxid)
        taxid = str(taxid)
        if args.taxonomy_id:
            if taxid != args.taxonomy_id:
                if taxid is None:
                    print(f"No taxid for refseq GI for {refseq_id}")
                # TODO: it's not a problem if the two taxids have an ancestor-descendant relationship -- look it up.
                print(f"\nNOTE: RefSeq ID {refseq_id} is associated with species-level Taxonomy ID {taxid}, not the provided Taxonomy ID {args.taxonomy_id}.\n")
                taxid = args.taxonomy_id
    elif args.ref_fasta or args.ref_gbff:
        if args.refseq:
            print("Error: --refseq cannot be used together with --ref_fasta or --ref_gbff", file=sys.stderr)
            sys.exit(1)
        if not args.ref_fasta or not args.ref_gbff:
            print("Error: both --ref_fasta and --ref_gbff must be provided when using one or the other", file=sys.stderr)
            sys.exit(1)
        ref_fasta = args.ref_fasta
        ref_gbff = args.ref_gbff
        refseq_id = ''
        assembly_id = ''
        ok, error_message = check_optional_file_readable(ref_fasta)
        if not ok:
            print(f"{error_message}\nPlease try again with a different file for --ref_fasta.", file=sys.stderr)
            sys.exit(1)
        ok, error_message = check_optional_file_readable(ref_gbff)
        if not ok:
            print(f"{error_message}\nPlease try again with a different file for --ref_gbff.", file=sys.stderr)
            sys.exit(1)
        taxid = args.taxonomy_id
        if not taxid:
            print("Error: when using --ref_fasta and --ref_gbff, -t/--taxonomy_id must be provided", file=sys.stderr)
            sys.exit(1)
        species = args.species
        if not species:
            species = ncbi.get_species_from_taxid(taxid)
            if not species:
                print(f"Could not find species for Taxonomy ID {taxid}.  Please provide a species name with --species.", file=sys.stderr)
                sys.exit(1)
    else:
        if args.taxonomy_id:
            species, taxid, refseq_id, assembly_id = lookup_refseq(ncbi, args.species, [], args.taxonomy_id)
        elif args.species:
            species, taxid, refseq_id, assembly_id = lookup_taxonomy_refseq(ncbi, args.species)
        else:
            species, taxid, refseq_id, assembly_id = get_species_taxonomy_refseq(ncbi)
        taxid = ncbi.get_species_level_taxid(taxid)

    # If --refseq/--ref_fasta and --workdir are given then accept defaults for other options
    is_interactive = not ((args.refseq or args.ref_fasta) and args.workdir)

    viral_usher_trees_name = get_viral_usher_trees_name(refseq_id, args.use_viral_usher_trees, is_interactive)
    if viral_usher_trees_name:
        update_tree_input = viral_usher_trees_raw_base + f"trees/{viral_usher_trees_name}/optimized.pb.gz"
        ref_fasta, ref_gbff = get_viral_usher_trees_refs(viral_usher_trees_name)
        if ref_fasta:
            if args.ref_fasta and args.ref_fasta != ref_fasta:
                print(f"Warning: ignoring --ref_fasta {args.ref_fasta} because it conflicts with the ref_fasta {ref_fasta} used for the viral_usher_trees tree {viral_usher_trees_name}.", file=sys.stderr)
            if args.ref_gbff and args.ref_gbff != ref_gbff:
                print(f"Warning: ignoring --ref_gbff {args.ref_gbff} because it conflicts with the ref_gbff {ref_gbff} used for the viral_usher_trees tree {viral_usher_trees_name}.", file=sys.stderr)
            print(f"Using reference fasta and GenBank file from viral_usher_trees {viral_usher_trees_name} config")
            refseq_id = ''
            assembly_id = ''
        update_metadata_input = viral_usher_trees_raw_base + f"trees/{viral_usher_trees_name}/metadata.tsv.gz"
    else:
        update_tree_input = ""
        update_metadata_input = ""
    nextclade_path, nextclade_columns = get_nextclade_path_columns(args.nextclade_dataset, species, is_interactive)
    min_length_proportion = get_min_length_proportion(args.min_length_proportion, is_interactive)
    max_N_proportion = get_max_N_proportion(args.max_N_proportion, is_interactive)
    max_parsimony = get_max_parsimony(args.max_parsimony, is_interactive)
    max_branch_length = get_max_branch_length(args.max_branch_length, is_interactive)
    fasta = get_user_fasta(args.fasta, is_interactive)
    metadata, metadata_date_column = get_user_metadata(args.metadata, args.metadata_date_column, is_interactive)
    workdir = get_workdir(args.workdir, is_interactive)

    viral_usher_version = importlib.metadata.version('viral_usher')
    config_contents = {
        "viral_usher_version": viral_usher_version,
        "refseq_acc": refseq_id,
        "refseq_assembly": assembly_id,
        "ref_fasta": ref_fasta,
        "ref_gbff": ref_gbff,
        "update_tree_input": update_tree_input,
        "update_metadata_input": update_metadata_input,
        "species": species,
        "taxonomy_id": taxid,
        "nextclade_dataset": nextclade_path,
        "nextclade_clade_columns": nextclade_columns,
        "min_length_proportion": min_length_proportion,
        "max_N_proportion": max_N_proportion,
        "max_parsimony": max_parsimony,
        "max_branch_length": max_branch_length,
        "extra_fasta": fasta,
        "extra_metadata": metadata,
        "extra_metadata_date_column": metadata_date_column,
        "taxonium_overlay_html": args.taxonium_overlay_html if args.taxonium_overlay_html else "",
        "workdir": os.path.abspath(workdir),
    }
    config_path = make_config(config_contents, workdir, refseq_id, taxid, args.config, is_interactive)

    print(f"\nReady to roll!  Next, try running 'viral_usher build --config {config_path}'\n")
