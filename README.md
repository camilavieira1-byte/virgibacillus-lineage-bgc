# Phylogeny shapes BGC content and prophage dynamics in *Virgibacillus*

Analysis code for:

> de Souza Vieira C, Lemos L, Kumazawa Morais D, Rosado A, Pylro V.
> *Phylogeny shapes biosynthetic gene cluster content and prophage dynamics in the genus* Virgibacillus.
> [JOURNAL], [YEAR]. doi:[ARTICLE DOI]

Archived release: [![DOI](https://zenodo.org/badge/DOI/ZENODO.DOI.svg)](https://doi.org/ZENODO.DOI)

## What this does

We analysed the 138 publicly available isolate genomes deposited under the name *Virgibacillus*
to ask how biosynthetic gene clusters (BGCs) and prophages are distributed across the genus.
The pipeline reclassifies the genomes against GTDB, collapses near-identical strains into clone
groups, predicts and clusters BGCs, tests whether cluster content follows the phylogeny,
reconstructs ancestral states for individual BGC classes, builds the pangenome of GTDB
*Virgibacillus*, and detects prophages and relates their similarity to host divergence.

Every number, table and figure in the paper is produced by the scripts here.

## Layout

```
code/
  00_metadata/            genome retrieval and BioSample/BioProject metadata
  03_qc/                  assembly statistics and completeness (seqkit, QUAST, CheckM2)
  04_taxonomy/            GTDB-Tk classification and database download
  05_ani_clone_groups/    pairwise ANI and clone-group dereplication (skani)
  06_antismash_novelty/   BGC prediction and MIBiG comparison (antiSMASH)
  07_bigscape/            gene cluster families (BiG-SCAPE, Pfam download)
  08_tree/                marker-gene alignment and ML tree (GTDB-Tk align, IQ-TREE)
  09_gcf_lineage_stats/   GCF matrices, Mantel, PERMANOVA, per-family tree tests
  10_pangenome/           gene prediction and orthogroups (Prodigal, OrthoFinder)
  11_genomad_prophages/   prophage detection and BGC/prophage co-location tests
  12_phage_analyses/      vOTUs, AAI groups and prophage-host divergence (PhaBOX2)
  13_figures/             all main and supplementary figures, ancestral state reconstruction
  14_supplement/          builds the supplementary tables workbook
MANIFEST.tsv              every script with its length and first comment line
environment_qc.yml        the conda environment used for the analyses and figures
software_versions.md      external tools and the versions used
```

Folder numbers are the pipeline order and match the working directories named in the Methods.
`.slurm` files are the cluster jobs for the heavy steps; everything else is a Python or shell
script that runs on their outputs.

## Running it

The analyses ran on the KAUST Ibex cluster (SLURM) in a project folder called
`virgibacillus_v2`. That path was replaced throughout by the placeholder
`/PATH/TO/virgibacillus_v2`. Point it at your own folder first:

```bash
grep -rl "/PATH/TO/virgibacillus_v2" code | xargs sed -i 's#/PATH/TO/virgibacillus_v2#'"$PWD"'#g'
```

Then recreate the environment and work through the folders in order:

```bash
conda env create -f environment_qc.yml -n virgi
conda activate virgi
```

Two scripts in `00_metadata/` query NCBI E-utilities, which asks callers to identify
themselves. Set your own address before running them:

```bash
export NCBI_EMAIL="you@example.org"
```

The external tools in `software_versions.md` (GTDB-Tk, antiSMASH, BiG-SCAPE, OrthoFinder,
geNomad, PhaBOX2 and their databases) are not part of the conda environment and must be
installed separately. The `.slurm` headers carry the resources each step needed.

## Data

Genome sequences are public GenBank assemblies and are not redistributed here; the accessions
are in Supplementary Table S1. The derived tables behind every figure and statistic are in the
supplementary tables deposited with the article, and are rebuilt by
`code/14_supplement/build_supp_tables.py`.

That script writes the full working table set, in the order the analyses were run. The workbook
submitted with the article drops the sheets the manuscript does not cite and renumbers the rest
in order of first mention, so the sheet names differ:

| sheet written by the script | table in the article |
|---|---|
| `S1_genomes` | Table S1 |
| `S2_BGC_regions` | Table S2 |
| `S5_cluster_stability` | Table S3 |
| `S7_lineage_robustness` | Table S4 |
| `S6_GCF_tree_pattern` (+ `_995`, `_B`, `_B_995`) | Table S5 |
| `S8a/b/c` pangenome | Table S7a/b/c |
| `S3_GCFs_c0.5`, `S4_genome_x_GCF_c0.5`, `S9a/b` prophages | not submitted (uncited) |

Table S6 of the article, the ancestral state reconstruction, is not produced by that script: it
is `ancestral_summary_<tag>.tsv`, written by `code/13_figures/asr_all_classes.py`.

## License

MIT, see `LICENSE`. The code is released for reuse and inspection; if it helps your work,
please cite the article above.
