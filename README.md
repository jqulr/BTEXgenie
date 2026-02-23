# BTEX-HMM: A database for the functional annotation of BTEX-degradation genes from isolate genomes and metagenomes

<a href="HMMs">
  <img src="img/toluene.png" alt="toluene button" width="60" />
</a>

This toolkit uses hmmscan from the HMMER suite to query BTEX-HMM profiles.
The following packages are needed for visualization and analysis:

- python 3.11
- pip
- hmmer 3.3 or newer
- pandas 2.0 or newer
- numpy 1.24 or newer
- matplotlib 3.7 or newer
- biopython 1.81 or newer
- circos 0.69 or newer

If your system already satisfies these requirements, you can move directly to running the BTEX-HMM scripts. Otherwise, you may install everything through Conda as shown below.

## Install via Conda
Confirm that a working Conda installation is available.

Create the BTEX-HMM environment and install all dependencies via *btex_env.yml*:

```bash
conda env create -n btex-hmm -f btex_env.yml
```

Activate the environment:

```bash
conda activate btex-hmm
```
## Usage
To run BTEX-HMMs on genomes, input can be either a single protein coding file or a directory of protein files. 

## Example with protein files in *test_genomes*
```bash
annotate-btex -p btexhmm/test_genomes \
              -o path/to/output_dir \
              --evalue 1e-5 \
              --cpus 8
```
**Outputs**
- `btex-hmm-summary.csv` contains the output from running each input file against all the HMMs. 
- `hmmscan_output` contains a sub-directory for each input file with the raw *.domtblout* output files produced before and after filtering by GA thresholds. 

**Notes:**
- The `annotate-btex` command expects amino acid FASTA files as input. (If starting from genome assemblies, protein-coding sequences should be predicted first using a gene caller such as Prodigal)
- By default, only single best-scoring HMM hit for each protein-coding region is reported. To retain all HMM hits per protein-coding region, user should enable the `--all-hits-per-protein` option. With this flag, multiple matches to different HMMs may be reported for the same protein.

## Visualizations
BTEX-HMM supports visualization of detected genes on Circos plot or KEGG pathway modules. 

**For visualization of hits on KEGG pathways:**

```bash
vis-btex --hmmscan /path/to/output_dir/btex_hmm_summary.csv \
         -s all \
         -o /path/to/vis-btex-outputs
```
 An example annotated pathway generated via the KEGG URL is available [here](https://www.kegg.jp/kegg-bin/show_pathway?map=map00623&multi_query=ko:K03268%20%23FFFFFF,%23FF0000%0Ako:K03381%20%23FF7F00,%23FF0000%0Ako:K07540%20%23377EB8,%23FF0000%0Ako:K15757%20%23FFFFFF,%23FF0000%0Ako:K15758%20%23FFFFFF,%23FF0000%0Ako:K15760%20%23377EB8,%23FF0000%0Ako:K15761%20%23FFFFFF,%23FF0000%0Ako:K15762%20%23377EB8,%23FF0000%0Ako:K15763%20%23377EB8,%23FF0000%0Ako:K15764%20%23377EB8,%23FF0000%0Ako:K15765%20%23377EB8,%23FF0000%0Ako:K16242%20%23377EB8,%23FF0000%0Ako:K16243%20%23377EB8,%23FF0000%0Ako:K16244%20%23FFFFFF,%23FF0000%0Ako:K16245%20%23FFFFFF,%23FF0000%0Ako:K16246%20%23377EB8,%23FF0000%0Ako:K16268%20%23E41A1C,%23FF0000%0Ako:K16269%20%23FFFFFF,%23FF0000%0Ako:K18089%20%23FFFFFF,%23FF0000%0Ako:K18090%20%23FFFFFF,%23FF0000).

**For visualization of hits on Circos:**

```bash
run-circos \
  --hits /home/juneq/Toluene_test/main_output/btex_hmm_summary.csv \
  --outdir /home/juneq/Toluene_test/main_output \
  --dna /home/juneq/Toluene_test/fastas/aromatoleum_aromaticum_ebn1.fasta \
  --genome "aromatoleum_aromaticum_ebn1" 
```

<!--## Visualization
For isolate genomes, BTEX-HMM hits can be rendered on a Circos plot together with a GenBank file describing the genomic regions containing the identified profiles.

Circos requires an etc configuration directory (for image settings, fonts, housekeeping). The `circos.py` helper auto-detects this using the `circos` executable on your PATH:

To obtain the correct etc directory (from a conda install):

```bash
ETCDIR=$(dirname "$(dirname "$(which circos)")")/etc
```

### Generate Circos Plot
The `circos.py` helper script builds a complete Circos project for a single genome, using BTEX-HMM hits and optional operon and pathway information.

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
