# BTEX-HMM: A database for functional annotations of isolate genomes or metagenomes

<a href="HMMs">
  <img src="img/toluene.png" alt="toluene button" width="60" />
</a>

This toolkit uses hmmsearch from the HMMER suite to query BTEX-HMM profiles. Ensure HMMER is installed and available on your PATH.
The following Python packages are needed for visualization and analysis:

python 3.11
pip
hmmer 3.3 or newer
pandas 2.0 or newer
numpy 1.24 or newer
matplotlib 3.7 or newer
biopython 1.81 or newer
circos 0.69 or newer

If your system already satisfies these requirements, you can move directly to running the BTEX-HMM scripts.
Otherwise, you may install everything through Conda or pip as shown below.

## Install via Conda
Confirm that a working Conda installation is available.

Clone the repository, then install HMMER:

```conda install -c bioconda hmmer```

Create the BTEX-HMM environment:

``` conda env create -n btex-hmm -f btex_env.yml ```

Activate the environment:

```conda activate btex-hmm```

## Visualization
For isolate genomes, BTEX-HMM hits can be rendered on a Circos plot together with a GenBank file describing the genomic regions containing the identified profiles.

Circos requires an etc configuration directory, supplied through the argument --etcdir, which provides the defaults for image settings, fonts and housekeeping parameters.

To obtain the correct etc directory:

After activating the Conda environment, set the etc path with:

```ETCDIR=$(dirname "$(dirname "$(which circos)")")/etc```

### Generate Circos Plot
The `circos.py` helper script builds a complete Circos project for a single genome, using BTEX-HMM hits and optional operon and pathway information.

#### Command-line help
Run:

```python circos.py --help```

to see the full summary of the commands.

usage: circos.py [-h] [--hits HITS] [--hmmsearch-outdir HMMSEARCH_OUTDIR] [--contig-lengths CONTIG_LENGTHS] --outdir OUTDIR --etcdir
                 ETCDIR --genome GENOME [--project-subdir PROJECT_SUBDIR] [--hmm-colors HMM_COLORS] [--only-hit-contigs] [--operon]
                 [--operon-defs OPERON_DEFS] [--dna DNA] [--pathway-map PATHWAY_MAP]

Generate Circos plot files for HMM hits on a single genome.

options:
  -h, --help            show this help message and exit
  --hits HITS           CSV: sample,hmm,hits,total_genes,hit_headers. If omitted, defaults to <hmmsearch-
                        outdir>/btex_hmm_summary.csv
  --hmmsearch-outdir HMMSEARCH_OUTDIR
                        Directory where hmmsearch.py wrote its outputs. Used to locate btex_hmm_summary.csv when --hits is not
                        given.
  --contig-lengths CONTIG_LENGTHS
                        TSV: sample<TAB>contig<TAB>length. If omitted, contig_length.tsv is generated under --outdir using --dna.
  --outdir OUTDIR       Output root dir
  --etcdir ETCDIR       Circos etc dir (for conf files)
  --genome GENOME       Sample/genome name to render
  --project-subdir PROJECT_SUBDIR
                        Subdir inside the genome folder
  --hmm-colors HMM_COLORS
                        Optional TSV: hmm<TAB>color (rgb 'r,g,b' or hex '#RRGGBB')
  --only-hit-contigs    Include only contigs that have >=1 hit
  --operon              Switch to operon-centric plotting mode with layered tracks for completeness.
  --operon-defs OPERON_DEFS
                        Required if --operon is used. TSV with header: operon_id,members,max_span_bp
  --dna DNA             Path to the genome nucleotide FASTA (.fna/.fasta). Used to generate contig lengths when --contig-lengths is
                        omitted and to populate GenBank sequences.
  --pathway-map PATHWAY_MAP
                        Optional TSV: hmm<TAB>pathway. Groups HMMs by pathway, each pathway has a different color. 

Example: 
To generate the circos plot for the **aromatoleum_aromaticum_ebn1** genome:
```python /home/juneq/Toluene-HMM/btexhmm/circos.py \
  --outdir /home/juneq/Toluene_test/main_output \
  --etcdir /home/juneq/.conda/envs/btex-hmm/etc \
  --dna /home/juneq/Toluene_test/fastas/aromatoleum_aromaticum_ebn1.fasta \
  --project-subdir /home/juneq/Toluene_test/main_output/circos_plots \
  --genome "aromatoleum_aromaticum_ebn1" \
  --pathway-map /home/juneq/hmm/archetypes/hmm_cutoffs/pathway_map.tsv
```
**Note:**The value passed to --genome must exactly match the entry in the Sample column of the btex_hmm_summary.csv file produced by the BTEX-HMM pipeline.

## things change
