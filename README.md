# BTEX-HMM: A database for the functional annotation of BTEX-degradation genes from isolate genomes and metagenomes

<a href="HMMs">
  <img src="img/toluene.png" alt="toluene button" width="60" />
</a>

This toolkit uses hmmscan from the HMMER suite to query BTEX-HMM profiles.
The following packages are needed for visualization and analysis:

- python 3.9 or newer
- pip
- Rscript
- hmmer 3.3 or newer
- prodigal
- kofamscan
- pandas 2.0 or newer
- numpy 1.24 or newer
- matplotlib 3.7 or newer
- biopython 1.81 or newer
- circos 0.69 or newer

If your system already satisfies these requirements, you can move directly to running the BTEX-HMM scripts. Otherwise, you may install everything through Conda as shown below.

## Install via Conda
Confirm that a working Conda installation is available. See [Conda installation](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) for more details.

### Development install
Create the BTEX-HMM environment and install all dependencies plus an editable install of the local repo via *btex_env.yml*:

```bash
cd /path/to/BTEX-HMMs
conda env create -n btex-hmm -f btex_env.yml
```

Activate the environment:

```bash
conda activate btex-hmm
```

Confirm that the KOfamScan executable is installed:

```bash
command -v exec_annotation
```

If this prints nothing, KOfamScan is not available on your `PATH`.

Database download:

The KOfam HMM database can be installed for users interested in the broad metabolic or degradation potential of their genomes

Database Download:

```bash
mkdir /path/to/databases/kofam
cd /path/to/databases/kofam

wget https://www.genome.jp/ftp/db/kofam/ko_list.gz
wget https://www.genome.jp/ftp/db/kofam/profiles.tar.gz

gunzip ko_list.gz
tar -xzf profiles.tar.gz
```

Setting the environment variable:

```bash
# shell scripts placed in activate.d folder to run automatically when env is activated
mkdir -p $CONDA_PREFIX/etc/conda/activate.d

cat > $CONDA_PREFIX/etc/conda/activate.d/kofam.sh <<'EOF'
export KOFAM_DB=/path/to/databases/kofam
EOF
```

Confirm the database path and executable:
```bash
conda deactivate
conda activate btex-hmm
echo $KOFAM_DB
command -v exec_annotation
```

`KOFAM_DB` should point to the database directory containing `profiles/` and `ko_list`.
`exec_annotation` is the KOfamScan program itself and must be installed separately from the database files. The Conda env files in this repo now include the `kofamscan` package so `exec_annotation` should be available after environment creation.

<!-- The `btex_env.yml` file includes `pip install -e .`, so it installs the local repository in editable mode and creates the `annotate-btex`, `vis-btex`, and `run-circos` command-line entry points in the environment's `bin/` directory. -->

<!-- ### End-user / release install
For a non-editable install that does not depend on the current working directory, build a wheel from the repo root:

```bash
cd /path/to/BTEX-HMMs
python -m build
```

Then create the environment from *btex_env_release.yml*:

```bash
conda env create -n btex-hmm -f btex_env_release.yml
```

This installs the packaged wheel from `./dist/` instead of the live source tree, which is more reproducible for end users. -->
## Usage
To run BTEX-HMMs, input should be a directory containing either genome DNA FASTA files or protein FASTA files.

## Example with protein files in *test_genomes*
```bash
annotate-btex -p btexhmm/test_genomes \
              -o path/to/output_dir \
              --evalue 1e-5 \
              --cpus 8
```

## Example with genome FASTA files
```bash
annotate-btex -p /path/to/genome_fastas \
              -o path/to/output_dir \
              -meta \
              --cpus 8
```
**Outputs**
- `btex_hmm_summary.csv` contains the output from running each input file against all the HMMs. 
- `prodigal_output` contains predicted protein FASTA files when genome DNA inputs are provided.
- `hmmscan_output` contains a sub-directory for each input file with the raw *.domtblout* output files produced before and after filtering by GA thresholds. 
- `kofam_abv_thres.tsv` is written for each sample when the KOfam step runs successfully.

**Notes:**
- The `annotate-btex` command accepts either genome DNA FASTA files or protein FASTA files, but all files in one run must be the same type.
- For genome DNA inputs, `annotate-btex` runs Prodigal first and requires either `-meta` or `-single`.
- BTEX-HMM reports the single best-scoring HMM hit per protein-coding region.

## Visualizations
BTEX-HMM supports visualization of detected genes on KEGG pathway modules or Circos plot.

**For visualization of hits on KEGG pathways:**

```bash
vis-btex --hmmscan /path/to/output_dir/btex_hmm_summary.csv \
         -s all \
         -o /path/to/vis-btex-outputs
```
 An example annotated pathway generated via the KEGG URL is available [here](https://www.kegg.jp/kegg-bin/show_pathway?map=map00642&multi_query=ko:K14748%20%23F08A8B,%23FF0000%0Ako:K14749%20%23F08A8B,%23FF0000%0Ako:K10700%20%23F08A8B,%23FF0000%20%23377EB8,%23FF0000%0Ako:K17048%20%23F08A8B,%23FF0000%20%23377EB8,%23FF0000%0Ako:K17049%20%23F08A8B,%23FF0000%20%23377EB8,%23FF0000%0Ako:K14579%20%23FFFFFF,%23FF0000%0Ako:K14580%20%23FFFFFF,%23FF0000
 ).

**Input:**
- `vis-btex` takes the btex_hmm_summary.csv file generated by `annotate-btex` as input. 
- Using `-s all` generates a pathway visualization that includes hits from all samples in the input file.
- You can also replace `all` with a specific sample name to visualize hits for a single sample. The sample name must exactly match the sample name in the input file.

**Output:**
- URLs for the visualization images are written to `KEGG_MAP_LINKS.txt` (By default, these links correspond to annotations on the KEGG pathway maps map00642, map00623, map00622, and map00362.). 
- The output directory also includes a legend file and a set of TSV files containing the KO equivalents of the BTEX-HMM hits.

**For visualization of hits on Circos:**

```bash
run-circos \
  --hmmscan /path/to/btex_hmm_summary.csv \
  --dna /test_genomes/Aromatoleum_bremense_PbN1T.fna \
  -o /path/to/output_dir \
  -s "Aromatoleum_bremense_PbN1T" 

  run-circos \
  --hmmscan /path/to/btex_hmm_summary.csv \
  --dna /path/to/genome.fna \
  --prodigal-gbk /path/to/sample_prodigal.gbk \
  --kofam-output /path/to/kofam_abv_thres.tsv \
  -o /path/to/outdir \
  -s sample_name
```

**Input:**
- `run-circos` takes btex_hmm_summary.csv together with the genome sequence file for a single sample. In the example above, the genome sequence file for Aromatoleum bremense PbN1T in the test_genomes folder is used.
- To display all hits identified by BTEX-HMMs, provide a sample name with -s that exactly matches the sample name in btex_hmm_summary.csv.

**Output:**
- `run-circos` generates .png and .pdf Circos plots, along with the Circos configuration files and supporting files used to produce the plots.

**Example Output:**

<!--## Visualization
For isolate genomes, BTEX-HMM hits can be rendered on a Circos plot together with a GenBank file describing the genomic regions containing the identified profiles.

Circos requires an etc configuration directory (for image settings, fonts, housekeeping). The `visualization_scripts/circos.py` helper auto-detects this using the `circos` executable on your PATH:

To obtain the correct etc directory (from a conda install):

```bash
ETCDIR=$(dirname "$(dirname "$(which circos)")")/etc
```

### Generate Circos Plot
The `visualization_scripts/circos.py` helper script builds a complete Circos project for a single genome, using BTEX-HMM hits and optional operon and pathway information.

#### Command-line help

```bash
run-circos --help
```

#### Example

To generate a Circos plot for the `aromatoleum_aromaticum_ebn1` genome:

```bash
run-circos \
  --hits /home/juneq/Toluene_test/main_output/btex_hmm_summary.csv \
  --outdir /home/juneq/Toluene_test/main_output \
  --dna /home/juneq/Toluene_test/fastas/aromatoleum_aromaticum_ebn1.fasta \
  --genome "aromatoleum_aromaticum_ebn1" 
```

**Notes:**

- `--hits` should point to the `btex_hmm_summary.csv` produced by hmmsearch.py.
- `--genome` must exactly match the `sample` value in `btex_hmm_summary.csv`.
- If `--contig-lengths` is omitted, lengths are generated from `--dna` into `<outdir>/<genome>_circos_plot/contig_length.tsv`. Provided contig-length files are copied into that directory for use.
- Circos outputs are written to `<outdir>/<genome>_circos_plot/`.
- The `--pathway-map` TSV supports an optional `color` column (`hmm<TAB>pathway<TAB>color`) to keep pathway colors stable across runs (colors can be `r,g,b` or `#RRGGBB`).
- If `--pathway-map` is omitted, `btexhmm/data/pathway_map.tsv` (relative to this repo) is used when available.
- Rendering runs automatically (`circos -conf circos.conf`) unless `--no-render` is provided. -->
