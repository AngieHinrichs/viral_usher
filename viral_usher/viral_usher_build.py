import argparse
import csv
import datetime
import gzip
import io
import json
import lzma
import os
import re
import shutil
import subprocess
import sys
import time
import zipfile
from Bio import SeqIO
from collections import namedtuple
from . import config
from . import ncbi_helper

update_tree_input = "optimized.pb.gz"
update_nextclade_input = "nextclade.clade.tsv"


def run_command(command, stdout_filename=None, stderr_filename=None, fail_ok=False):
    """Run a command, complain and raise if it fails (unless fail_ok, then just return False)"""
    success = False
    stdout = None
    stderr = None
    if stdout_filename:
        stdout = open(stdout_filename, 'w')
    if stderr_filename:
        stderr = open(stderr_filename, 'w')
    try:
        subprocess.run(command, check=True, stdout=stdout, stderr=stderr)
        success = True
    except subprocess.CalledProcessError as e:
        if fail_ok:
            pass
        else:
            if e.stderr:
                print(f"{command[0]} failed: {e.stderr}\nFailed command: {' '.join(command)}", file=sys.stderr)
            elif stderr_filename:
                print(f"{command[0]} failed, see {stderr_filename}\nFailed command: {' '.join(command)}", file=sys.stderr)
            else:
                print(f"{command[0]} failed\nFailed command: {' '.join(command)}", file=sys.stderr)
            sys.exit(1)
    if stdout:
        stdout.close()
    if stderr:
        stderr.close()
    return success


def run_command_get_stdout(command):
    """Run a command and return its stdout output."""
    try:
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return result.stdout.decode('utf-8')
    except subprocess.CalledProcessError as e:
        print(f"Command failed: {e.stderr}\nFailed command: {' '.join(command)}", file=sys.stderr)
        sys.exit(1)


def start_timing(message):
    print(message)
    return time.perf_counter()


def finish_timing(start_time):
    elapsed_time = time.perf_counter() - start_time
    print(f"... done in {elapsed_time:.1f}s")


def normalize_segment(segment):
    # Many viruses have segments S (small), M (medium), and L (large).  The vast majority of records use S, M or L
    # but there are a few oddballs that spell it out or use middle instead of medium.  Change those to S/M/L.
    if segment and len(segment) > 1:
        seg_lower = segment.lower()
        if seg_lower == 'small':
            segment = 'S'
        elif seg_lower == 'large':
            segment = 'L'
        elif seg_lower == 'medium' or seg_lower == 'middle':
            segment = 'M'
    return segment


def segment_from_gbff(gbff_path, acc):
    """Search for segment in gbff's record for acc"""
    segment = None
    with open(gbff_path, 'r') as gbff_in:
        for record in SeqIO.parse(gbff_in, 'genbank'):
            if record.id == acc:
                for feature in record.features:
                    if "segment" in feature.qualifiers:
                        segment = normalize_segment(feature.qualifiers["segment"][0])
                        print(f"Found segment {segment} for {acc} in {gbff_path}")
                        break
                break
    return segment


def unpack_refseq_zip(refseq_zip, refseq_acc):
    """Extract and rename the fasta and gbff files from the RefSeq zip."""
    start_time = start_timing(f"Unpacking {refseq_zip}...")
    with zipfile.ZipFile(refseq_zip, 'r') as zip_ref:
        fasta_path = 'refseq.fasta'
        gbff_path = 'refseq.gbff'
        length = 0
        fasta_found = False
        gbff_found = False
        for name in zip_ref.namelist():
            if name.endswith('.fna'):
                with io.TextIOWrapper(zip_ref.open(name, 'r'), encoding='utf-8') as fasta_in, open(fasta_path, 'w') as fasta_out:
                    for record in SeqIO.parse(fasta_in, 'fasta'):
                        if record.id == refseq_acc:
                            record.description = record.id
                            SeqIO.write(record, fasta_out, 'fasta')
                            length = len(record.seq)
                            fasta_found = True
            elif name.endswith('.gbff'):
                with zip_ref.open(name, 'r') as gbff_in, open(gbff_path, 'wb') as gbff_out:
                    shutil.copyfileobj(gbff_in, gbff_out)
                gbff_found = True
        if not fasta_found:
            raise ValueError(f"Failed to find .fna file with {refseq_acc} in {refseq_zip}.")
        if not gbff_found:
            raise ValueError(f"Failed to find .gbff file in {refseq_zip}.")
    # Search for segment in refseq.gbff's record for refseq_acc
    segment = segment_from_gbff(gbff_path, refseq_acc)
    finish_timing(start_time)
    return fasta_path, gbff_path, length, segment


def check_ref_fasta_gbff(ref_fasta, ref_gbff):
    """Make sure the given fasta and gbff files exist and contain exactly one sequence with the same accession.
    Return (ref_acc, ref_length, ref_segment)"""
    ref_acc = None
    ref_length = 0
    with open(ref_fasta, 'r') as fasta_in:
        records = list(SeqIO.parse(fasta_in, 'fasta'))
        if len(records) != 1:
            print(f"Error: expected exactly one sequence in {ref_fasta}, found {len(records)}.", file=sys.stderr)
            sys.exit(1)
        ref_acc = records[0].id
        ref_length = len(records[0].seq)
    if not ref_acc:
        print(f"Error: failed to find sequence in {ref_fasta}.", file=sys.stderr)
        sys.exit(1)
    # Check that the gbff has one record with the same accession
    with open(ref_gbff, 'r') as gbff_in:
        records = list(SeqIO.parse(gbff_in, 'genbank'))
    if len(records) != 1:
        print(f"Error: expected exactly one record in {ref_gbff}, found {len(records)}.", file=sys.stderr)
        sys.exit(1)
    if records[0].id != ref_acc:
        print(f"Error: failed to find record with accession {ref_acc} in {ref_gbff}.", file=sys.stderr)
        sys.exit(1)
    segment = segment_from_gbff(ref_gbff, ref_acc)
    return ref_acc, ref_length, segment


def get_reference(config_contents, ncbi):
    """Make sure config specifies either refseq_acc and refseq_assembly, or ref_fasta and ref_gbff, but not both.
    If refseq_acc is given, download the RefSeq genome and annotations.  Otherwise make sure the given files exist.
    Return (ref_acc, ref_fasta, ref_gbff, ref_length, ref_segment)"""
    refseq_acc = config_contents.get('refseq_acc')
    assembly_id = config_contents.get('refseq_assembly')
    ref_fasta = config_contents.get('ref_fasta')
    ref_gbff = config_contents.get('ref_gbff')
    if refseq_acc:
        # Download the RefSeq genome and annotations
        refseq_zip = f"{refseq_acc}.zip"
        start_time = start_timing(f"Downloading RefSeq {refseq_acc} (Assembly {assembly_id}) genome to {refseq_zip}...")
        ncbi.download_refseq(assembly_id, refseq_zip)
        finish_timing(start_time)
        refseq_fasta, refseq_gbff, refseq_length, refseq_segment = unpack_refseq_zip(refseq_zip, refseq_acc)
        return refseq_acc, refseq_fasta, refseq_gbff, refseq_length, refseq_segment
    else:
        # Use the user-provided reference genome and annotations; make sure the files are consistent
        if not os.path.exists(ref_fasta):
            print(f"Error: specified ref_fasta {ref_fasta} not found.", file=sys.stderr)
            sys.exit(1)
        if not os.path.exists(ref_gbff):
            print(f"Error: specified ref_gbff {ref_gbff} not found.", file=sys.stderr)
            sys.exit(1)
        ref_acc, ref_length, ref_segment = check_ref_fasta_gbff(ref_fasta, ref_gbff)
        return ref_acc, ref_fasta, ref_gbff, ref_length, ref_segment


def scan_ncbi_virus_metadata(ncbi_virus_metadata):
    """Get just a couple attributes that are useful for filtering sequences: length and segment."""
    acc_to_length_segment = {}
    with open(ncbi_virus_metadata, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            accession = row.get("accession")
            length = int(row.get("length", 0))
            segment = normalize_segment(row.get("segment", ""))
            acc_to_length_segment[accession] = (length, segment)
    return acc_to_length_segment


def get_genbank_metadata(ncbi, taxid, no_genbank):
    """Get metadata from NCBI Virus API and acc_to_length_segment mapping."""
    if no_genbank:
        return "/dev/null", {}
    start_time = start_timing(f"Querying NCBI Virus API for metadata for all GenBank sequences for taxid {taxid}...")
    ncbi_virus_metadata = "ncbi_virus_metadata.csv"
    ncbi.query_ncbi_virus_metadata(taxid, ncbi_virus_metadata)
    acc_to_length_segment = scan_ncbi_virus_metadata(ncbi_virus_metadata)
    finish_timing(start_time)
    return ncbi_virus_metadata, acc_to_length_segment


def get_genbank_fasta(update, no_genbank, ncbi, taxid, starting_tree_accessions,
                      acc_to_length_segment, min_length, max_N_proportion, ref_segment):
    """Download GenBank genomes for the given Taxonomy ID: only the new ones if --update, otherwise all of them"""
    if no_genbank:
        return "/dev/null", 0, 0
    got_genbank = False
    genbank_zip = f"genbank_{taxid}.zip"
    if (update):
        new_accessions = find_new_accessions(starting_tree_accessions, acc_to_length_segment, min_length, ref_segment)
        if len(new_accessions) > 0:
            start_time = start_timing(f"Downloading {len(new_accessions)} new GenBank genomes for taxid {taxid} to {genbank_zip}...")
            ncbi.download_genbank_accessions(new_accessions, genbank_zip)
            got_genbank = True
            finish_timing(start_time)
        else:
            print("No new GenBank accessions found, skipping GenBank download.")
    else:
        start_time = start_timing(f"Downloading all GenBank genomes for taxid {taxid} to {genbank_zip}...")
        ncbi.download_genbank(taxid, genbank_zip)
        got_genbank = True
        finish_timing(start_time)
    if got_genbank:
        genbank_fasta, gb_count, filtered_count = unpack_genbank_zip(genbank_zip, acc_to_length_segment, min_length, max_N_proportion, ref_segment)
    else:
        genbank_fasta = "/dev/null"
        gb_count = 0
        filtered_count = 0
    return genbank_fasta, gb_count, filtered_count


def passes_seq_filter(record, length, max_N_proportion):
    """Fasta record passes the filter if it has at most max_N_proportion of Ns"""
    return record.seq.count('N') / length <= max_N_proportion


def filter_genbank_sequences(fasta_in, fasta_out, min_length, max_N_proportion, acc_to_length_segment, ref_segment, fasta_out_path):
    # Filter fasta file and chop descriptions to get accession as name in nextclade
    start_time = start_timing(f"Filtering GenBank sequences by length >= {min_length} and proportion of Ns <= {max_N_proportion}...")
    record_count = 0
    passed_count = 0
    for record in SeqIO.parse(fasta_in, 'fasta'):
        record_count += 1
        (length, segment) = acc_to_length_segment.get(record.id, (0, None))
        if length == 0:
            length = len(record.seq)
        # Skip if segment is given but is different from ref_segment or sequence is shorter than min_length
        if length < min_length or (segment and segment != ref_segment):
            continue
        # Filter by proportion of ambiguous characters
        if not passes_seq_filter(record, length, max_N_proportion):
            continue
        passed_count += 1
        # Chop sequence name after first space to limit to accession only
        record.description = record.id
        SeqIO.write(record, fasta_out, 'fasta')
    finish_timing(start_time)
    print(f"Processed {record_count} sequences from GenBank and wrote {passed_count} sequences passing filters to {fasta_out_path}.")
    return record_count, passed_count


def unpack_genbank_zip(genbank_zip, acc_to_length_segment, min_length, max_N_proportion, ref_segment):
    """Process the fasta and data_report.jsonl files from the GenBank zip,
    filtering the fasta by length and proportion of ambiguous characters and
    keeping only the accession part of sequence names, and extracting basic
    metadata from data_report.jsonl into much more compact TSV.
    Compress fasta with lzma (xz) and TSV with gzip"""
    fasta_out_path = 'genbank.fasta.xz'
    with zipfile.ZipFile(genbank_zip, 'r') as zip_ref:
        namelist = zip_ref.namelist()
        genomic_fna_paths = [name for name in namelist if name.endswith('genomic.fna')]
        if len(genomic_fna_paths) != 1:
            raise ValueError(f"Expected exactly one genomic.fna file in {genbank_zip}, found: {genomic_fna_paths}")
        with io.TextIOWrapper(zip_ref.open(genomic_fna_paths[0], 'r'), encoding='utf-8') as fasta_in, lzma.open(fasta_out_path, 'wt') as fasta_out:
            genbank_count, passed_count = filter_genbank_sequences(fasta_in, fasta_out, min_length, max_N_proportion, acc_to_length_segment, ref_segment, fasta_out_path)
    return fasta_out_path, genbank_count, passed_count


def open_maybe_decompress(filename, mode='rt', encoding='utf-8'):
    """
    Open a file, automatically using gzip or lzma decompression based on the file extension.
    Supports .gz, .xz, and uncompressed files.
    """
    if filename.endswith('.gz'):
        return gzip.open(filename, mode, encoding=encoding)
    elif filename.endswith('.xz'):
        return lzma.open(filename, mode, encoding=encoding)
    else:
        return open(filename, mode, encoding=encoding)


def record_to_fasta_bytes(record):
    """Convert a SeqRecord to FASTA-encoded bytes (utf-8)"""
    buf = io.StringIO()
    SeqIO.write(record, buf, 'fasta')
    return buf.getvalue().encode('utf-8')


def sanitize_name(name):
    """Replace characters that could cause trouble, such as spaces or Newick special characters,
    with underscores."""
    for c in ['[', ']', '(', ')', ':', ';', ',', "'", ' ']:
        name = name.replace(c, '_')
    return name


def get_extra_fasta_names(extra_fasta):
    """Return None if extra_fasta is not given.  If it is, then read the extra fasta file and return a set of
    sanitized sequence names that may or may not appear in the tree."""
    if extra_fasta:
        extra_fasta_names = set()
        with open_maybe_decompress(extra_fasta) as ef:
            for record in SeqIO.parse(ef, 'fasta'):
                extra_fasta_names.add(sanitize_name(record.description))
        return extra_fasta_names
    else:
        return None


def align_sequences(ref_fasta, extra_fasta, genbank_fasta, ref_acc, min_length, max_N_proportion):
    """Run nextclade to align the filtered sequences to the reference, and pipe its output to faToVcf and gzip."""
    msa_vcf_gz = 'msa.vcf.gz'
    nextclade_err_txt = 'nextclade.align.err.log'
    start_time = start_timing(f"Aligning sequences with nextclade and converting to VCF, writing to {msa_vcf_gz}...")
    # Set up the pipeline: | nextclade | faToVcf | gzip > msa.vcf.gz
    nextclade_cmd = [
        'nextclade', 'run', '--input-ref', ref_fasta, '--include-reference', 'true', '--output-fasta', '/dev/stdout'
    ]
    fatovcf_cmd = ['faToVcf', '-includeNoAltN', '-ref=' + ref_acc, 'stdin', 'stdout']
    with gzip.open(msa_vcf_gz, 'wb') as vcf_out, open(nextclade_err_txt, 'wb') as nextclade_stderr:
        try:
            # Start nextclade
            nextclade_proc = subprocess.Popen(nextclade_cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=nextclade_stderr)
            # Start faToVcf, reading from nextclade's stdout
            fatovcf_proc = subprocess.Popen(fatovcf_cmd, stdin=nextclade_proc.stdout, stdout=subprocess.PIPE)
            nextclade_proc.stdout.close()  # Allow nextclade to receive SIGPIPE if faToVcf exits
            # Send fasta input(s) to nextclade_proc's stdin.  If extra_fasta duplicates sequences
            # also found in genbank_fasta, prefer the ones in extra_fasta.
            extra_ids = dict()
            if extra_fasta:
                with open_maybe_decompress(extra_fasta) as ef:
                    for record in SeqIO.parse(ef, 'fasta'):
                        record.description = sanitize_name(record.description)
                        record.id = record.description
                        if record.id in extra_ids:
                            print(f"Fasta ID '{record.id}' has already been used in {extra_fasta}, discarding duplicate.", file=sys.stderr)
                        else:
                            extra_ids[record.id] = True
                            length = len(record.seq)
                            if length < min_length or not passes_seq_filter(record, length, max_N_proportion):
                                continue
                            nextclade_proc.stdin.write(record_to_fasta_bytes(record))
            # genbank_fasta has already been filtered for length & Ns, but discard duplicates
            with open_maybe_decompress(genbank_fasta) as gf:
                for record in SeqIO.parse(gf, 'fasta'):
                    if record.id in extra_ids:
                        print(f"GenBank ID '{record.id}' was already used in {extra_fasta}, discarding duplicate.", file=sys.stderr)
                    else:
                        nextclade_proc.stdin.write(record_to_fasta_bytes(record))
            nextclade_proc.stdin.close()
            # Write faToVcf's output to gzip file
            for chunk in iter(lambda: fatovcf_proc.stdout.read(8192), b''):
                vcf_out.write(chunk)
            fatovcf_proc.stdout.close()
            fatovcf_proc.wait()
            nextclade_proc.wait()
        except subprocess.CalledProcessError as e:
            print(f"nextclade | faToVcf failed: {e.stderr}", file=sys.stderr)
            raise
    finish_timing(start_time)
    start_time = start_timing(f"Compressing {nextclade_err_txt}...")
    # gzip nextclade.align.err.log (note: using gzip.open above did not work, it was written uncompressed.)
    run_command(['gzip', '-f', nextclade_err_txt])
    # Count the number of samples in the #CHROM header line of msa.vcf.gz
    aligned_count = 0
    with gzip.open(msa_vcf_gz, 'rt') as vcf_in:
        for line in vcf_in:
            if line.startswith('#CHROM'):
                header = line.strip().split('\t')
                aligned_count = len(header) - 9  # Subtract fixed columns
                break
    finish_timing(start_time)
    return msa_vcf_gz, aligned_count


def make_empty_tree():
    """Create an empty tree file to be used as a starting point for usher-sampled."""
    empty_tree_path = 'empty_tree.nwk'
    with open(empty_tree_path, 'w') as f:
        f.write("()\n")
    return empty_tree_path


def run_usher_sampled(tree, vcf):
    """Build the tree.  Use docker --platform linux/amd64 so this will work even on Mac with ARM CPU. """
    pb_out = 'usher_sampled.pb.gz'
    start_time = start_timing(f"Running usher-sampled on {vcf}...")
    tree_flag = "-t" if tree.endswith('.nwk') else "-i"
    command = ['usher-sampled', '-A', '-e', '5', tree_flag, tree, '-v', vcf, '-o', pb_out,
               '--optimization_radius', '0', '--batch_size_per_process', '100']
    run_command(command, stdout_filename='usher-sampled.out.log', stderr_filename='usher-sampled.err.log')
    finish_timing(start_time)
    return pb_out


def run_matoptimize(pb_file, vcf_file, update):
    """Run matOptimize to clean up after usher-sampled"""
    pb_out = 'optimized.unfiltered.pb.gz'
    start_time = start_timing(f"Running matOptimize on {pb_file}...")
    command_no_vcf = ['matOptimize', '-m', '0.00000001', '-M', '1', '-i', pb_file, '-o', pb_out]
    # Unless we're doing an update, in which case the VCF covers only the new sequences,
    # first try with VCF, which is less tested but should give better results because it includes
    # info about which bases are ambiguous or N.  That info is lost when usher imputes values.
    if update:
        run_command(command_no_vcf, stdout_filename='matOptimize.out.log', stderr_filename='matOptimize.err.log')
    else:
        command_vcf = command_no_vcf.copy()
        command_vcf.append("-v")
        command_vcf.append(vcf_file)
        if not run_command(command_vcf, stdout_filename='matOptimize.out.log', stderr_filename='matOptimize.err.log', fail_ok=True):
            finish_timing(start_time)
            start_time = start_timing("matOptimize with VCF failed, trying again without VCF...")
            run_command(command_no_vcf, stdout_filename='matOptimize.out.log', stderr_filename='matOptimize.err.log')
    finish_timing(start_time)
    return pb_out


def run_matutils_filter(opt_unfiltered_tree, max_parsimony, max_branch_length):
    """Run matUtils extract to filter sequences/branches and collapse tree post-matOptimize"""
    pb_out = update_tree_input
    sample_names_out = 'tree_samples.txt'
    start_time = start_timing(f"Running matUtils to filter (--max-parsimony {max_parsimony} --max-branch-length {max_branch_length})...")
    command = ['matUtils', 'extract', '-i', opt_unfiltered_tree,
               '--max-parsimony', str(max_parsimony),
               '--max-branch-length', str(max_branch_length),
               '--collapse-tree',
               '--used-samples', sample_names_out,
               '-o', pb_out]
    run_command(command, stdout_filename='matUtils.filter.out.log', stderr_filename='matUtils.filter.err.log')
    finish_timing(start_time)
    # Run matUtils summary to get the final number of samples in the tree
    start_time = start_timing("Counting samples in the final tree...")
    command = ['matUtils', 'summary', '-i', pb_out]
    summary_output = run_command_get_stdout(command)
    for line in summary_output.split('\n'):
        if line.startswith('Total Samples in Tree:'):
            tree_tip_count = int(line.split(':')[1].strip())
            break
    finish_timing(start_time)
    return pb_out, sample_names_out, tree_tip_count


def get_existing_nextclade_assignments(update, starting_tree_accessions, nextclade_path):
    """Get existing Nextclade assignments, and the number of columns to expect, from the old TSV file about to be replaced."""
    if not update or not nextclade_path:
        return {}, None
    existing_assignments = {}
    nextclade_tsv = update_nextclade_input if os.path.exists(update_nextclade_input) else update_nextclade_input + ".gz"
    existing_column_count = None
    if os.path.exists(nextclade_tsv):
        with open_maybe_decompress(nextclade_tsv) as tsv_in:
            header_cols = tsv_in.readline().rstrip('\n').split('\t')
            # 'dynamic' might have added clade and/or phenotype, motif etc. columns; keep all of them
            col_count = len(header_cols) - 1
            if existing_column_count is None:
                existing_column_count = col_count
            elif existing_column_count != col_count:
                raise ValueError(f"Inconsistent Nextclade column count in {nextclade_tsv}: {existing_column_count} vs {col_count}")
            for line in tsv_in:
                fields = line.rstrip('\n').split('\t')
                seq_name = sanitize_name(fields[0])
                clades = fields[1:]
                existing_assignments[seq_name] = clades
    else:
        print(f"Nextclade TSV file {nextclade_tsv} not found. Only new items will have clade assignments.", file=sys.stderr)
    return existing_assignments, existing_column_count


def run_nextclade(nextclade_path, nextclade_clade_columns, existing_nextclade_assignments, existing_column_count, genbank_fasta, extra_fasta, ref_fasta):
    """Run nextclade to assign sequences to clades.  Return dict mapping accession to clade."""
    if not nextclade_path:
        return {}, ""
    nextclade_tsv = update_nextclade_input
    # Nextclade --output-columns-selection does not accept individual customClades names; instead use 'dynamic'
    # to get all of them.  See https://github.com/nextstrain/nextclade/issues/1667 .
    nextclade_col_list = nextclade_clade_columns.split(',')
    if len(nextclade_col_list) > 1:
        output_columns = 'seqName,clade,dynamic'
    else:
        output_columns = 'seqName,clade'
    start_time = start_timing(f"Running nextclade to assign clades (dataset: {nextclade_path})...")
    command = ['nextclade', 'run',
               '--dataset-name', nextclade_path,
               '--output-columns-selection', output_columns,
               '--output-tsv', nextclade_tsv,
               genbank_fasta, ref_fasta]
    if extra_fasta:
        command.append(extra_fasta)
    run_command(command, stdout_filename='nextclade.clade.out.log', stderr_filename='nextclade.clade.err.log')
    finish_timing(start_time)
    # Read the TSV file and return a dict mapping sequence names to clades
    start_time = start_timing(f"Reading nextclade assignments from {nextclade_tsv} to add to metadata...")
    new_assignments = {}
    with open(nextclade_tsv, 'r') as tsv_in:
        header_cols = tsv_in.readline().rstrip('\n').split('\t')
        # 'dynamic' might have added clade and/or phenotype, motif etc. columns; keep all of them
        nextclade_clade_columns = ','.join(header_cols[1:])
        nextclade_column_count = len(header_cols) - 1
        if existing_column_count is not None and nextclade_column_count != existing_column_count:
            print(f"Existing nextclade annotations had {existing_column_count} columns but new annotations have {nextclade_column_count} columns ({header_cols[1:]})." +
                  f"\nRemove {update_nextclade_input} and run again, to rebuild with the latest nextclade dataset '{nextclade_path}'", file=sys.stderr)
            sys.exit(1)
        for line in tsv_in:
            fields = line.rstrip('\n').split('\t')
            seq_name = sanitize_name(fields[0])
            clades = fields[1:]
            new_assignments[seq_name] = clades
    if existing_column_count is not None:
        nextclade_assignments = existing_nextclade_assignments
        nextclade_assignments.update(new_assignments)
    else:
        nextclade_assignments = new_assignments
    finish_timing(start_time)
    if existing_column_count is not None:
        # Append the existing assignments to the new little update nextclade.clade.tsv so they're not lost.
        with open(nextclade_tsv, 'a') as tsv_out:
            for acc, clade_columns in existing_nextclade_assignments.items():
                if acc not in new_assignments:
                    print(acc + "\t" + "\t".join(clade_columns), file=tsv_out)
    return nextclade_assignments, nextclade_clade_columns


def get_extra_metadata(extra_metadata):
    """If extra_metadata is given, parse it and return its header columns (skip the first which must be the name,
    matching the names in extra_fasta) and a dict mapping sanitized sequence name from the first column to its row
    (skipping the first column)."""
    if not extra_metadata:
        return [], {}
    extra_metadata_cols = []
    extra_metadata_dict = {}
    with open(extra_metadata, 'r') as em:
        header = em.readline().rstrip('\n').split('\t')
        extra_metadata_cols = header[1:]
        for line in em:
            fields = line.rstrip('\n').split('\t')
            name = sanitize_name(fields[0])
            extra_metadata_dict[name] = {k: v for k, v in zip(header[1:], fields[1:])}
    return extra_metadata_cols, extra_metadata_dict


def analyze_extra_columns(extra_metadata_cols, header, extra_metadata_date_column):
    """If extra metadata columns are given, try to match them to NCBI metadata columns; return a list of leftovers
    and a dict mapping matched column name to NCBI column name."""
    extra_added_cols = []
    extra_mapped_cols = {}
    if not extra_metadata_cols:
        return extra_added_cols, extra_mapped_cols
    header_set = set(header)
    if 'country_location' in header_set:
        header_set.add('country')
        header_set.add('location')
        header_set.remove('country_location')
    for col in extra_metadata_cols:
        if col == extra_metadata_date_column:
            extra_mapped_cols['date'] = col
        elif col in header_set:
            extra_mapped_cols[col] = col
        else:
            extra_added_cols.append(col)
    if extra_added_cols:
        print(f"Note: the following extra metadata columns were not found in NCBI metadata and will be added as new columns: {', '.join(extra_added_cols)}", file=sys.stderr)
    return extra_added_cols, extra_mapped_cols


def get_numerical_date(date):
    """Date should be YYYY-MM-DD or YYYY-MM or YYYY.  Return year.fraction"""
    if not date:
        return 0.0
    groups = re.match('^([0-9]{4})(-([0-9]{2})(-([0-9]{2}))?)?$', date)
    if not groups:
        return 0.0
    year = int(groups[1])
    if groups[3]:
        month = int(groups[3])
        day = int(groups[5]) if groups[5] else 15
    else:
        month = 7
        day = 1
    date_obj = datetime.date(year, month, day)
    day_of_year = date_obj.timetuple().tm_yday
    fraction = day_of_year / 365.25
    return year + fraction


def fudge_isolate(isolate, accession, date, country, year):
    """If isolate looks like a full /-separated descriptor then use it.
    Otherwise try to approximate it from available metadata bits and pieces."""
    full = None
    if '/' in isolate:
        full = isolate
    else:
        if country and isolate and year:
            full = '/'.join([country, isolate, year])
        elif country and year:
            full = '/'.join([country, year])
        elif isolate and year:
            full = '/'.join([isolate, year])
        elif isolate:
            full = isolate
    return full


def make_display_name(isolate, strain, date, accession, country):
    """Make a display name with whatever info is available in isolate or strain, the country and date in addition to accession."""
    if not isolate:
        isolate = strain
    elif strain and len(isolate) < len(strain):
        # Sometimes people use nice strain like "MVi/Marseille.FRA/21.19/20[D8]" but goofy isolate
        # like "9061923678/9051834309" for MZ031240.1.  If strain is longer, use it.
        isolate = strain
    groups = re.match('^([0-9]{4})(-[0-9]{2})?(-[0-9]{2})?$', date)
    year = None
    if groups:
        year = groups[1]
    if not date:
        date = '?'
    full = fudge_isolate(isolate, accession, date, country, year)
    if full:
        name = '|'.join([full, accession, date])
    else:
        name = '|'.join([accession, date])
    return name


def get_nextclade_column_list(nextclade_clade_columns):
    if not nextclade_clade_columns:
        return []
    nextclade_col_list = nextclade_clade_columns.split(',')
    nextclade_col_list = [col if col.startswith('nextclade') or col.startswith('Nextclade') else 'nextclade_' + col for col in nextclade_col_list]
    return nextclade_col_list


NcbiMetadataIndex = namedtuple('NcbiMetadataIndex',
                               ['accession', 'isolate', 'strain', 'date', 'country_location', 'segment'])


def make_metadata_header(ncbi_header, nextclade_clade_columns, remove_segment, extra_fasta_names, extra_added_cols, extra_mapped_cols):
    """Make the metadata header line.  Return an index into the NCBI metadata columns and the number of added nextclade columns."""
    if ncbi_header:
        accession_idx = ncbi_header.index('accession')
        isolate_idx = ncbi_header.index('isolate')
        strain_idx = ncbi_header.index('strain')
        date_idx = ncbi_header.index('date')
        country_loc_idx = ncbi_header.index('country_location')
        # In metadata header, make separate columns for country and location:
        header_mod = ncbi_header[:country_loc_idx] + ['country', 'location'] + ncbi_header[country_loc_idx+1:]
        # Rename 'strain' to 'gb_strain' because I add 'strain' at the beginning below, following the
        # nextstrain convention of 'strain' meaning the full descriptor
        header_mod[header_mod.index('strain')] = 'gb_strain'
        # Remove 'segment' if ref has no segment
        segment_idx = header_mod.index('segment')
        if remove_segment:
            header_mod.pop(segment_idx)
        metadata_header = 'strain\t' + '\t'.join(header_mod)
        midx = NcbiMetadataIndex(accession=accession_idx, isolate=isolate_idx, strain=strain_idx, date=date_idx,
                                 country_location=country_loc_idx, segment=segment_idx)
    else:
        metadata_header = 'strain'
        for col in extra_mapped_cols:
            metadata_header += '\t' + col
        midx = None
    # Add numerical date column
    metadata_header += '\tnum_date'
    nextclade_col_list = []
    if nextclade_clade_columns:
        nextclade_col_list = get_nextclade_column_list(nextclade_clade_columns)
        metadata_header += '\t' + '\t'.join(nextclade_col_list)
    if ncbi_header and extra_fasta_names is not None:
        # Add a column indicating the source of the sequence: GenBank or user-provided extra fasta
        metadata_header += '\tsource'
    if extra_added_cols:
        metadata_header += '\t' + '\t'.join(extra_added_cols)
    metadata_header += '\n'
    return metadata_header, midx, len(nextclade_col_list)


def make_metadata_and_rename(row, midx, remove_segment, nextclade_assignments, nextclade_col_count, extra_fasta_names, extra_added_col_count):
    """Make a metadata line from the given NCBI metadata row and nextclade clades, and a renaming line.  Also return numerical date."""
    accession = row[midx.accession].strip()
    isolate = sanitize_name(row[midx.isolate].strip())
    strain = sanitize_name(row[midx.strain].strip())
    date = row[midx.date].strip()
    country_loc = row[midx.country_location].strip()
    # Make separate columns for country and location
    country, location = country_loc.split(':', 1) if ':' in country_loc else (country_loc, '')
    country = sanitize_name(country.strip())
    location = location.strip()
    row_mod = row[:midx.country_location] + [country, location] + row[midx.country_location+1:]
    if remove_segment:
        row_mod.pop(midx.segment)
    # Add numerical date column (year.fraction)
    if date and re.match('^[0-9]{4}(-[0-9]{2}(-[0-9]{2})?)?$', date):
        num_date = get_numerical_date(date)
        row_mod.append(f"{num_date:.6f}")
    else:
        num_date = None
        row_mod.append("")
    name = make_display_name(isolate, strain, date, accession, country)
    metadata_line = name + '\t' + '\t'.join(row_mod)
    if nextclade_col_count > 0:
        clades = nextclade_assignments.get(accession, [''] * nextclade_col_count)
        metadata_line += '\t' + '\t'.join(clades)
    if extra_fasta_names is not None:
        # All of these are from GenBank
        metadata_line += '\tGenBank'
    if extra_added_col_count > 0:
        metadata_line += '\t' + '\t'.join([''] * extra_added_col_count)
    metadata_line += '\n'
    rename_line = "\t".join([accession, name]) + '\n'
    return metadata_line, rename_line, name, num_date


def add_extra_metadata_rows(m_out, midx, metadata_header, sample_names, nextclade_assignments, nextclade_col_count,
                            extra_fasta_names, extra_added_cols, extra_mapped_cols, extra_metadata_rows):
    """Add metadata lines for user-provided extra fasta sequences that are in the tree (i.e. in sample_names)."""
    if not extra_fasta_names:
        return
    sample_names_set = set()
    with open(sample_names, 'r') as sn:
        for line in sn:
            sample_names_set.add(line.strip())
    # For each name in extra_fasta_names that is in sample_names_set (i.e. in the tree), add a metadata line:
    # Name, followed by empty (or merged if user-provided) fields for all NCBI-derived metadata columns prior
    # to nextclade clades, then nextclade clades if any, then 'user-provided' in source column, then extra
    # user-provided metadata columns if any.
    metadata_header_cols = metadata_header.rstrip('\n').split('\t')
    # Subtract 1 each for strain, num_date, and source (if applicable)
    ncbi_col_count = len(metadata_header_cols) - 2 - nextclade_col_count - len(extra_added_cols)
    if midx is not None:
        ncbi_col_count -= 1
    date_min = None
    date_max = None
    for name in extra_fasta_names:
        if name not in sample_names_set:
            continue
        row = extra_metadata_rows.get(name, {})
        metadata = name
        if extra_mapped_cols:
            for col in metadata_header_cols[1:1+ncbi_col_count]:
                if col in extra_mapped_cols and extra_mapped_cols[col] in row:
                    metadata += '\t' + row[extra_mapped_cols[col]]
                else:
                    metadata += '\t'
        else:
            metadata += '\t' + '\t'.join([''] * ncbi_col_count)
        # Add numerical date column
        if extra_mapped_cols and 'date' in extra_mapped_cols and extra_mapped_cols['date'] in row:
            num_date = get_numerical_date(row[extra_mapped_cols['date']])
            metadata += f'\t{num_date:.6f}' if num_date else '\t'
        elif extra_added_cols and 'date' in extra_added_cols and 'date' in row:
            num_date = get_numerical_date(row['date'])
            metadata += f'\t{num_date:.6f}' if num_date else '\t'
        else:
            num_date = None
            metadata += '\t'
        if num_date is not None:
            if date_min is None or num_date < date_min:
                date_min = num_date
            if date_max is None or num_date > date_max:
                date_max = num_date
        if nextclade_col_count > 0:
            if name in nextclade_assignments:
                clades = nextclade_assignments[name]
            else:
                clades = [''] * nextclade_col_count
            metadata += '\t' + '\t'.join(clades)
        if midx is not None:
            metadata += '\tuser-provided'
        if extra_added_cols:
            for col in extra_added_cols:
                if col in row:
                    metadata += '\t' + row[col]
                else:
                    metadata += '\t'
        metadata += '\n'
        m_out.write(metadata)
    return date_min, date_max


def finalize_metadata(ncbi_virus_metadata, nextclade_assignments, nextclade_clade_columns, ref_segment, sample_names,
                      extra_fasta_names, extra_metadata, extra_metadata_date_column):
    """Make a file for renaming sequence names in tree to include country, isolate name and date when available.
    Prepare metadata for taxonium."""
    rename_out = "rename.tsv"
    metadata_out = "metadata.tsv.gz"
    start_time = start_timing(f"Finalizing sequence names for display and metadata, writing to {metadata_out}...")
    remove_segment = not ref_segment
    date_min = None
    date_max = None
    extra_metadata_cols, extra_metadata_rows = get_extra_metadata(extra_metadata)
    # Read header from ncbi_virus_metadata and add header columns if applicable
    with open(ncbi_virus_metadata, 'rt') as csv_in:
        line = csv_in.readline()
        if line:
            header = line.rstrip('\n').split(',')
        else:
            header = []
    extra_added_cols, extra_mapped_cols = analyze_extra_columns(extra_metadata_cols, header, extra_metadata_date_column)
    metadata_header, midx, nextclade_col_count = make_metadata_header(header, nextclade_clade_columns, remove_segment,
                                                                      extra_fasta_names, extra_added_cols, extra_mapped_cols)
    # Use a subprocess to get input from grep -Fwf {sample_names} {ncbi_virus_metadata}
    grep_cmd = ['grep', '-Fwf', sample_names, ncbi_virus_metadata]
    with subprocess.Popen(grep_cmd, stdout=subprocess.PIPE, text=True) as grep_proc, \
            open(rename_out, 'w') as r_out, gzip.open(metadata_out, 'wt') as m_out:
        # No header for rename_out; write header for metadata_out
        m_out.write(metadata_header)

        # Process filtered lines from grep of NCBI metadata
        reader = csv.reader(grep_proc.stdout, delimiter=',')
        for row in reader:
            metadata, rename, name, num_date = make_metadata_and_rename(row, midx, remove_segment, nextclade_assignments,
                                                                        nextclade_col_count, extra_fasta_names, len(extra_added_cols))
            if num_date is not None:
                if date_min is None or num_date < date_min:
                    date_min = num_date
                if date_max is None or num_date > date_max:
                    date_max = num_date
            # Write output to both renaming file and metadata file
            r_out.write(rename)
            m_out.write(metadata)
        # Add metadata lines for user-provided extra fasta sequences
        extra_date_min, extra_date_max = \
            add_extra_metadata_rows(m_out, midx, metadata_header, sample_names, nextclade_assignments, nextclade_col_count,
                                    extra_fasta_names, extra_added_cols, extra_mapped_cols, extra_metadata_rows)
        if extra_date_min is not None and (date_min is None or extra_date_min < date_min):
            date_min = extra_date_min
        if extra_date_max is not None and (date_max is None or extra_date_max > date_max):
            date_max = extra_date_max
    finish_timing(start_time)
    return metadata_out, rename_out, date_min, date_max, extra_added_cols, extra_mapped_cols


def rename_seqs(pb_in, rename_tsv):
    pb_out = "viz.pb.gz"
    start_time = start_timing(f"Renaming sequences in {pb_in} to make {pb_out}...")
    command = ['matUtils', 'mask', '-i', pb_in, '--rename-samples', rename_tsv, '-o', pb_out]
    run_command(command, stdout_filename='matUtils.rename.out.log', stderr_filename='matUtils.rename.err.log')
    finish_timing(start_time)
    return pb_out


def dump_newick(pb_in):
    """Dump the phylogenetic tree in Newick format."""
    newick_out = "viz.nwk"
    start_time = start_timing(f"Writing tree in Newick format to {newick_out}.gz...")
    command = ['matUtils', 'extract', '-i', pb_in, '--write-tree', newick_out]
    run_command(command, stdout_filename='/dev/null', stderr_filename='/dev/null')
    command = ['gzip', '-f', newick_out]
    run_command(command, stdout_filename='/dev/null', stderr_filename='/dev/null')
    finish_timing(start_time)


def get_header(tsv_in):
    with gzip.open(tsv_in, 'rt') as tsv:
        header = tsv.readline().split('\t')
        for idx, field in enumerate(header):
            header[idx] = field.strip()
    return header


def make_taxonium_config(date_min, date_max, nextclade_clade_columns, extra_fasta_names, no_genbank, extra_added_cols,
                         extra_mapped_cols):
    """Make a config file with a color gradient from date_min to date_max."""
    config_out = "taxonium_config.json"
    if no_genbank:
        color_by_options = []
        if extra_mapped_cols:
            color_by_options += ["meta_" + col for col in extra_mapped_cols.keys()]
        if 'date' in extra_mapped_cols:
            color_by_options.append("meta_num_date")
    else:
        color_by_options = ["meta_country", "meta_location", "meta_num_date", "meta_host", "meta_serotype"]
        if extra_fasta_names is not None:
            color_by_options.append("meta_source")
    if nextclade_clade_columns:
        color_by_options += ["meta_" + col for col in get_nextclade_column_list(nextclade_clade_columns)]
    if extra_added_cols:
        color_by_options += ["meta_" + col for col in extra_added_cols]
    color_by_options += ["genotype", "None"]
    config = {"colorRamps": {"meta_num_date": {"scale": [[date_min, "#0000F8"], [date_max, "#F80000"]]}},
              "customNames": {"meta_num_date": "Date (numeric)",
                              "meta_gb_strain": "Strain"},
              "colorBy": {"colorByOptions": color_by_options}}
    with open(config_out, 'w') as f:
        json.dump(config, f)
    return config_out


def usher_to_taxonium(pb_in, metadata_in, ref_gbff, tip_count, species, ref_acc, date_min, date_max,
                      nextclade_clade_columns, extra_fasta_names, no_genbank, extra_added_cols, extra_mapped_cols):
    jsonl_out = "tree.jsonl.gz"
    start_time = start_timing(f"Running usher_to_taxonium to make {jsonl_out}...")
    columns = ','.join(get_header(metadata_in))
    title = f"{tip_count} {species} sequences from GenBank aligned to {ref_acc}"
    config = make_taxonium_config(date_min, date_max, nextclade_clade_columns, extra_fasta_names, no_genbank,
                                  extra_added_cols, extra_mapped_cols)
    command = ['usher_to_taxonium', '--input', pb_in, '--metadata', metadata_in,
               '--columns', columns, '--title', title, '--config_json', config,
               '--genbank', ref_gbff, '--output', jsonl_out]
    if not run_command(command, stdout_filename='utt.out.log', stderr_filename='utt.err.log', fail_ok=True):
        finish_timing(start_time)
        start_time = start_timing("usher_to_taxonium failed with --genbank, trying again without --genbank...")
        command = ['usher_to_taxonium', '--input', pb_in, '--metadata', metadata_in,
                   '--columns', columns, '--title', title, '--config_json', config,
                   '--output', jsonl_out]
        run_command(command, stdout_filename='utt.out.log', stderr_filename='utt.err.log')
    finish_timing(start_time)
    return jsonl_out


def write_output_stats(ref_acc, ref_length, gb_count, filtered_count, aligned_count, tree_tip_count):
    output_stats = [["ref_acc", "ref_length", "gb_count", "filtered_count", "aligned_count", "tree_tip_count"],
                    [ref_acc, ref_length, gb_count, filtered_count, aligned_count, tree_tip_count]]
    print("Writing output stats to output_stats.tsv.")
    with open("output_stats.tsv", "w") as f:
        for row in output_stats:
            f.write("\t".join(map(str, row)) + "\n")


def get_tree_names(pb_in):
    """Return a set containing the names of the tips in pb_in."""
    tree_names = set()
    start_time = start_timing(f"Extracting tree names from {pb_in}...")
    command = ["matUtils", "extract", "-i", pb_in, "--used-samples", "tree_samples.txt"]
    run_command(command, stdout_filename="/dev/null", stderr_filename="/dev/null")
    with open("tree_samples.txt", "r") as f:
        for line in f:
            tree_names.add(line.strip())
    finish_timing(start_time)
    return tree_names


def find_new_accessions(tree_names, acc_to_length_segment, min_length, ref_segment):
    """Return a list of GenBank accessions in acc_to_length_segment that are not in the current tree and probably wouldn't be filtered later."""
    all_accessions = set(acc_to_length_segment.keys())
    # Exclude accessions whose metadata has a segment that is different from ref's
    exclude_accessions = {acc for acc, (length, seg) in acc_to_length_segment.items() if length < min_length or seg and seg != ref_segment}
    new_accessions = all_accessions - tree_names - exclude_accessions
    return list(new_accessions)


def get_new_extra_fasta(tree_names, extra_fasta):
    """Get a new extra_fasta file containing only sequences not already in the tree."""
    if not extra_fasta:
        return ""
    now = datetime.datetime.now()
    new_extra_fasta = "new_extra." + now.strftime("%Y-%m-%d_%H-%M-%S") + ".fasta"
    start_time = start_timing(f"Filtering {extra_fasta} to get new sequences not already in tree, writing to {new_extra_fasta}...")
    new_count = 0
    with open_maybe_decompress(extra_fasta) as ef, open(new_extra_fasta, 'w') as nef:
        for record in SeqIO.parse(ef, 'fasta'):
            record.description = sanitize_name(record.description)
            record.id = record.description
            if record.id not in tree_names:
                SeqIO.write(record, nef, 'fasta')
                new_count += 1
    finish_timing(start_time)
    if new_count == 0:
        print(f"No new sequences found in {extra_fasta}, skipping it.")
        return ""
    else:
        print(f"Wrote {new_count} new sequences to {new_extra_fasta}.")
        return new_extra_fasta


def run_nextclade_on_existing(ncbi, starting_tree_accessions, nextclade_path, nextclade_clade_columns,
                              extra_fasta, ref_fasta):
    """Download GenBank sequences that are already in the tree.  Run nextclade on those plus extra_fasta."""
    old_genbank_zip = "genbank_old.zip"
    start_time = start_timing("Downloading and unpacking GenBank genomes that were already in the tree, for nextclade...")
    ncbi.download_genbank_accessions(list(starting_tree_accessions), old_genbank_zip)
    old_genbank_fasta, _, _ = unpack_genbank_zip(old_genbank_zip, {}, 0, 1.0, "")
    finish_timing(start_time)
    print("Running nextclade on genomes that were already in the tree.")
    run_nextclade(nextclade_path, nextclade_clade_columns, {}, None, old_genbank_fasta, extra_fasta, ref_fasta)
    os.remove(old_genbank_zip)
    if not os.path.exists(update_nextclade_input):
        print(f"Expected {update_nextclade_input} to be created but it's not found.", sys.stderr)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(prog="viral_usher_build")
    parser.add_argument("--config", type=str, required=True, help="Path to config file input")
    parser.add_argument("--no_genbank", action="store_true", help="Skip downloading sequences from GenBank; use only extra_fasta")
    parser.add_argument("--update", action="store_true",
                        help="Add only new sequences to existing tree instead of building tree from scratch")

    args = parser.parse_args()

    print(f"Running pipeline with config file: {args.config}")
    config_contents = config.parse_config(args.config)
    taxid = config_contents['taxonomy_id']
    nextclade_path = config_contents.get('nextclade_dataset', '')
    nextclade_clade_columns = config_contents.get('nextclade_clade_columns', 'clade')
    min_length_proportion = float(config_contents.get('min_length_proportion', config.DEFAULT_MIN_LENGTH_PROPORTION))
    max_N_proportion = float(config_contents.get('max_N_proportion', config.DEFAULT_MAX_N_PROPORTION))
    max_parsimony = int(config_contents.get('max_parsimony', str(config.DEFAULT_MAX_PARSIMONY)))
    max_branch_length = int(config_contents.get('max_branch_length', str(config.DEFAULT_MAX_BRANCH_LENGTH)))
    extra_fasta = config_contents.get('extra_fasta', '')
    extra_metadata = config_contents.get('extra_metadata', '')
    extra_metadata_date_column = config_contents.get('extra_metadata_date_column', '')
    species = config_contents.get('species', None)
    ncbi = ncbi_helper.NcbiHelper()
    if not species:
        start_time = start_timing(f"Looking up species name for Taxonomy ID {taxid}...")
        species = ncbi.get_species_from_taxid(taxid)
        finish_timing(start_time)

    ref_acc, ref_fasta, ref_gbff, ref_length, ref_segment = get_reference(config_contents, ncbi)
    min_length = int(ref_length * min_length_proportion)

    extra_fasta_names = get_extra_fasta_names(extra_fasta)

    starting_tree_accessions = []
    if (args.update):
        # Make sure we have required inputs
        if not os.path.exists(update_tree_input):
            print(f"No existing {update_tree_input} found, cannot do update.  Run without --update to build tree from scratch.", file=sys.stderr)
            sys.exit(1)
        starting_tree_accessions = get_tree_names(update_tree_input)
        if nextclade_path and not os.path.exists(update_nextclade_input) and not os.path.exists(update_nextclade_input + ".gz"):
            print(f"No existing {update_nextclade_input} found; recreating it.")
            run_nextclade_on_existing(ncbi, starting_tree_accessions, nextclade_path, nextclade_clade_columns,
                                      extra_fasta, ref_fasta)

    # Get metadata from NCBI Virus API
    ncbi_virus_metadata, acc_to_length_segment = get_genbank_metadata(ncbi, taxid, args.no_genbank)

    genbank_fasta, gb_count, filtered_count = get_genbank_fasta(args.update, args.no_genbank, ncbi, taxid, starting_tree_accessions,
                                                                acc_to_length_segment, min_length, max_N_proportion, ref_segment)
    if (args.update):
        extra_fasta = get_new_extra_fasta(starting_tree_accessions, extra_fasta)

    # The core of the pipeline: align sequences, build tree, finalize metadata
    msa_vcf, aligned_count = align_sequences(ref_fasta, extra_fasta, genbank_fasta, ref_acc, min_length, max_N_proportion)
    if args.update:
        if aligned_count == 0:
            preopt_tree = "usher_sampled.pb.gz"
            shutil.copy(update_tree_input, preopt_tree)
        else:
            preopt_tree = run_usher_sampled(update_tree_input, msa_vcf)
    else:
        empty_tree = make_empty_tree()
        preopt_tree = run_usher_sampled(empty_tree, msa_vcf)
    opt_unfiltered_tree = run_matoptimize(preopt_tree, msa_vcf, args.update)
    opt_tree, sample_names, tree_tip_count = run_matutils_filter(opt_unfiltered_tree, max_parsimony, max_branch_length)
    existing_nextclade_assignments, existing_column_count = get_existing_nextclade_assignments(args.update, starting_tree_accessions, nextclade_path)
    nextclade_assignments, nextclade_clade_columns = run_nextclade(nextclade_path, nextclade_clade_columns,
                                                                   existing_nextclade_assignments, existing_column_count,
                                                                   genbank_fasta, extra_fasta, ref_fasta)
    metadata_tsv, rename_tsv, date_min, date_max, extra_added_cols, extra_mapped_cols = \
        finalize_metadata(ncbi_virus_metadata, nextclade_assignments, nextclade_clade_columns, ref_segment,
                          sample_names, extra_fasta_names, extra_metadata, extra_metadata_date_column)
    viz_tree = rename_seqs(opt_tree, rename_tsv)
    dump_newick(viz_tree)
    write_output_stats(ref_acc, ref_length, gb_count, filtered_count, aligned_count, tree_tip_count)
    usher_to_taxonium(viz_tree, metadata_tsv, ref_gbff, tree_tip_count, species, ref_acc, date_min, date_max,
                      nextclade_clade_columns, extra_fasta_names, args.no_genbank, extra_added_cols, extra_mapped_cols)


if __name__ == "__main__":
    main()
