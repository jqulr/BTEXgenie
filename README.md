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
annotate-btex -g btexhmm/test_genomes \
              -o path/to/output_dir \
              --evalue 1e-5 \
              --cpus 8
```
> [!NOTE]
> For proper parsing of genomic coordinates, protein files produced from Prodigal are needed.

## Example with genome FASTA files
```bash
annotate-btex -g /path/to/genome_fastas \
              -o path/to/output_dir \
              --meta \
              --cpus 8
```

> [!NOTE]
> For FASTA sequence inputs, the program will run gene-calling with Prodigal with either the --meta or --single flag specified as input. 

**Main outputs**

1. `btex_hmm_summary.csv`  
   Reports individual BTEX HMM hits, including the matched HMM, the threshold used, the hit score, and the protein sequence header for the corresponding gene.

2. `btex_hmm_summary_counts.csv`  
   Summarizes BTEX HMM hits by HMM, reporting hit counts instead of individual protein sequence headers.

3. `prodigal_output/`  
   Generated when genome DNA sequences are used as input. This directory contains one subdirectory per genome and may include:  
   ` {genome}_prodigal.gbk `  
   ` {genome}.faa `  
   ` {genome}_kofam_abv_thres.tsv `, produced when the KOfam step runs successfully for that sample.

4. `hmmscan_output/`  
   Contains one subdirectory per input file with raw `.domtblout` results generated before and after filtering by GA thresholds.

5. `log_file_annotate-btex.txt`  
   Records run progress as well as detailed warning and error messages.

**Notes:**
- The `annotate-btex` command accepts either genome DNA FASTA files or protein FASTA files, but all files in one run must be the same type.
- BTEX-HMM reports the single best-scoring HMM hit per protein-coding region.

## Visualizations
BTEX-HMM supports visualization of detected genes on KEGG pathway modules or Circos plot.

**For visualization of BTEX-HMM hits on KEGG pathways:**

```bash
vis-btex --hmmscan /path/to/output_dir/btex_hmm_summary.csv \
         -o /path/to/vis-btex-outputs
```
 An example annotated pathway generated via the KEGG URL is available [here](https://www.kegg.jp/kegg-bin/show_pathway?map=map00642&multi_query=ko:K14748%20%23F08A8B,%23FF0000%0Ako:K14749%20%23F08A8B,%23FF0000%0Ako:K10700%20%23F08A8B,%23FF0000%20%23377EB8,%23FF0000%0Ako:K17048%20%23F08A8B,%23FF0000%20%23377EB8,%23FF0000%0Ako:K17049%20%23F08A8B,%23FF0000%20%23377EB8,%23FF0000%0Ako:K14579%20%23FFFFFF,%23FF0000%0Ako:K14580%20%23FFFFFF,%23FF0000
 ).

**For visualization of all KOfam hits on all KEGG pathways:**

```bash
vis-btex --hmmscan /path/to/output_dir/btex_hmm_summary.csv \
         -o /path/to/vis-btex-outputs \
         --pathway 00623
```

**Inputs:**
- `vis-btex` takes the `btex_hmm_summary.csv` file generated by `annotate-btex` as input. 
- Using `-s {genome}` generates the visualization for hits specific to {genome}. Note, {genome} must exactly match the sample name in `btex_hmm_summary.csv`. 

**Outputs:**
1. {output_dir}/`KEGG_MAP_LINKS.txt` 
   Contains the URLs for visualizing hits on KEGG pathways. 
   
   > [!Note]
   > By default, these links correspond to the KEGG pathways: map00642 (xylene degradation), map00623 (toluene degradation), map00622 (ethylbenzene degradation), and map00362 (benzoate degradation)
   > Users interested in other pathways should ensure the tool is ran successfully against the KOfam HMM database. Then other pathways such as map00626 (napthalene degradation) can be visualized properly.

2. {output_dir}/`sample_color_legend.tsv` 
   Contains legend file that provides the color code for each input genome. 

2. {output_dir}/`split_color...` 
   TODOTODO

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

