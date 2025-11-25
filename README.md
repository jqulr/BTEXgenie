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

Circos requires an etc configuration directory, supplied through the argument `--etcdir`, which provides the defaults for image settings, fonts and housekeeping parameters.

To obtain the correct etc directory (from a conda install):

```bash
ETCDIR=$(dirname "$(dirname "$(which circos)")")/etc
```

### Generate Circos Plot
The `circos.py` helper script builds a complete Circos project for a single genome, using BTEX-HMM hits and optional operon and pathway information.

#### Command-line help

```bash
python btexhmm/circos.py --help
```

#### Example

To generate a Circos plot for the `aromatoleum_aromaticum_ebn1` genome:

```bash
python /home/juneq/Toluene-HMM/btexhmm/circos.py \
  --outdir /home/juneq/Toluene_test/main_output \
  --etcdir /home/juneq/.conda/envs/btex-hmm/etc \
  --dna /home/juneq/Toluene_test/fastas/aromatoleum_aromaticum_ebn1.fasta \
  --project-subdir /home/juneq/Toluene_test/main_output/circos_plots \
  --genome "aromatoleum_aromaticum_ebn1" \
  --pathway-map /home/juneq/hmm/archetypes/hmm_cutoffs/pathway_map.tsv
```

**Notes:**

- `--genome` must exactly match the `sample` value in `btex_hmm_summary.csv`.
- If `--contig-lengths` is omitted, lengths are generated from `--dna` into `<outdir>/contig_length.tsv`.
