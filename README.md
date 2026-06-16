# viral_usher

**viral_usher** is a command-line tool to set up and run a pipeline to build an [UShER](https://usher-wiki.readthedocs.io/en/latest/) tree for a new viral species (or type, subtype, etc.) using genomes downloaded from [NCBI](https://www.ncbi.nlm.nih.gov/).

---

## 🔧 Features

- Subcommands:
  - `init`: Generate a config file (interactive or via command line options)
  - `build`: Download sequences and build a tree, guided by the config file
- Uses a container (either [Apptainer](https://apptainer.org/) or [Docker](https://www.docker.com/) for portability to laptops, servers, or cloud platforms

---

## 📦 Installation

0. Install prerequisites (if not already installed)
- Either [Apptainer](https://apptainer.org/) or [Docker](https://www.docker.com/)
- [Python](https://www.python.org/) version 3.11 or later (we highly recommend using an environment manager such as [venv](https://docs.python.org/3/library/venv.html), [miniconda](https://www.anaconda.com/docs/getting-started/miniconda/main), [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) etc.)

1. Install with pip (again, we highly recommended using an environment manager):
    ```bash
    pip install viral_usher
    ```

---

## 🚀 Quickstart

### Create a config file with `viral_usher init`
If you want to start by just naming a virus, and let viral_usher interactively help you identify the right reference sequence, Taxonomy ID etc., then simply run
   ```bash
   viral_usher init
   ```
and reply to the prompts.

Alternatively, if you already know your parameters, then you can skip the interactive stuff by passing in command line options.  Run `viral_usher --help` to get a listing of options.  Here is an example that builds a tree for the Chikungunya virus using RefSeq [NC_004162.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_004162.2), all genomes available from GenBank for the Taxonomy ID associated with NC_004162.2 (Taxonomy ID [37124](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Tree&id=37124&lvl=3&lin=f&keep=1&srchmode=1&unlock)), plus additional sequences from example/hypothetical_chikungunya.fasta (in this repository):
   ```bash
   git clone https://github.com/AngieHinrichs/viral_usher.git
   cd viral_usher
   viral_usher init \
       --refseq NC_004162.2 \
       --workdir chikungunya \
       --fasta example/hypothetical_chikungunya.fasta \
       --config chikungunya/config.toml
   ```

### Build a tree using config file with `viral_usher build`:
Continuing the Chikungunya virus example:
   ```bash
   viral_usher build --config chikungunya/config.toml
   ```

That's all!  viral_usher will create the following files in workdir (`chikungunya` in our example):
- a tree in UShER protobuf format (optimized.pb.gz)
- a metadata file in TSV format (metadata.tsv.gz)
- a Taxonium tree file that you can view using https://taxonium.org/ (tree.jsonl.gz)

To view the example Chikungunya virus tree in Taxonium, [click here](https://taxonium.org/?protoUrl=https%3A%2F%2Fraw.githubusercontent.com%2FAngieHinrichs%2Fviral_usher%2Frefs%2Fheads%2Fmain%2Fexample%2Ftree.jsonl.gz&xType=x_dist).  Type or copy-paste "hypothetical" into Taxonium's Name search input to find the sequences from example/hypothetical_chikungunya.fasta.

---

## 🛠️ Troubleshooting

### MacOS Docker virtualization framework settings

If your `viral_usher build ...` command fails on MacOS with this error message:
```
usher-sampled failed, see usher-sampled.err.log
Failed command: usher-sampled -A -e 5 -t empty_tree.nwk -v msa.vcf.gz -o usher_sampled.pb.gz --optimization_radius 0 --batch_size_per_process 100
docker container angiehinrichs/viral_usher:v0.9.0 failed with exit code 1
```
and the file `usher-sampled.err.log` ends with this error:
```
assertion failed [rem_idx != -1]: Unable to find existing allocation for shared memory segment to unmap
```
then your Docker virtualization settings might need to be changed.  In the Docker app, click on the settings icon (gear), select the "General" category in the left column, and scroll down to the Virtual Machine Options section.  Under "Choose Virtual Machine Manager (VMM)", change the selection from "Apple Virtualization Framework" to "Docker VMM".  Then try running your `viral_usher build ...` command again.  For more on this topic see [issue #24](https://github.com/AngieHinrichs/viral_usher/issues/24) .

---

## 🧪 Development

```bash
# Clone the repo
git clone https://github.com/AngieHinrichs/viral_usher.git
cd viral_usher

# Install dev dependencies
pip install -e .[dev]

# Run tests
pytest

# Run build with a locally modified docker image
docker build --platform linux/amd64 -t viral_usher_test . \
    && viral_usher build --docker_image viral_usher_test --config path/to/config.toml

# Get a shell in the locally modified docker image to try out different commands
docker run --platform linux/amd64 --rm -it -u $(id -u):$(id -g) --net=host \
    -v $(pwd):/data \
    viral_usher_test
```
