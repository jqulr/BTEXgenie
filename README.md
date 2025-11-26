# BTEX-HMM: A database for functional annotations of isolate genomes or metagenomes

<a href="HMMs">
  <img src="img/toluene.png" alt="toluene button" width="60" />
</a>

This toolkit uses hmmsearch from the HMMER suite to query BTEX-HMM profiles. Ensure HMMER is installed and available on your PATH.
The following packages are needed for visualization and analysis:

- python 3.11
- pip
- hmmer 3.3 or newer
- pandas 2.0 or newer
- numpy 1.24 or newer
- matplotlib 3.7 or newer
- biopython 1.81 or newer
- circos 0.69 or newer

If your system already satisfies these requirements, you can move directly to running the BTEX-HMM scripts. Otherwise, you may install everything through Conda or pip as shown below.

## Install via Conda
Confirm that a working Conda installation is available.

Clone the repository, then install HMMER:

```bash
conda install -c bioconda hmmer
```

Create the BTEX-HMM environment:

```bash
conda env create -n btex-hmm -f btex_env.yml
```

Activate the environment:

```bash
conda activate btex-hmm
```

## Visualization
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
- Rendering runs automatically (`circos -conf circos.conf`) unless `--no-render` is provided.
