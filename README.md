
# Code Availability for Programmed Manipulation of RNA Targets using Human Argonaute 2

This repository contains the Snakemake pipeline, Python scripts, R scripts, and reference data used in the RNAseq analysis shown in the paper, as well as scripts and workflows for the MARIA de-immunization analysis.

---

## Repository Overview

| Directory                  | Contents                                                                                                                                                                                                                                                            |
| :------------------------- | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `1_rnaseq-final-analysis/` | Snakemake pipeline, scripts, reference files, and example data for the final RNAseq analysis. Includes raw variant calls and read count data for ALTER and CURE data which can be used to run the pipeline as an example or to reproduce data from the publication. |
| `2_MARIA_ Deimmunization/` | Scripts for de-immunization analysis using MARIA. Includes an example tutorial workflow.                                                                                                                                                                            |
| `3_initial-submission/`    | Code from the initial submission, prior to adapting to Snakemake. Preserved for transparency.                                                                                                                                                                       |

---

# RNAseq Analysis

The pipeline takes variant tables and Salmon quantification outputs from GATK/STAR/Salmon (derived from the raw read data as described in the Supplementary Information) and performs: off-target hit identification and filtering, transcript alignment and sequence context analysis, differential expression analysis (DESeq2, GSEA, GO enrichment), and cross-experiment miRNA seed alignment analysis. All code for this section is in `1_rnaseq-final-analysis/`.

## Dependencies

The pipeline is managed by Snakemake. Python and R dependencies are handled automatically through two Conda environment definition files in `1_rnaseq-final-analysis/env/`. These files can be reviewed for full information on libraries and versions used. Snakemake must be installed separately.

The initial repository with the reference files (see Reference Data) is ~2 GB. If the pipeline is run with the example data the final size will be ~5.5 GB after creation of the conda environments and analysis files. Unless otherwise indicated, code has been tested with the dependency versions listed below. Hardware: Apple M3 Pro, 32 GB unified memory, macOS Tahoe 26.5.1.

### Software

|   Name    | Version |
| :-------: | :-----: |
|   Conda   | 26.5.0  |
| Snakemake |  9.6.2  |
|  Python   | 3.12.13 |
|     R     |  4.4.3  |

## Reference Data

A custom combined genome reference (`GENCODEv48-ALTER-CURE`) was built from the GENCODE Human v48 primary assembly with added sequences and features for transfected plasmids. Due to GitHub's file size limits the full reference file directory `GENCODEv48-ALTER-CURE/` is distributed as a Release asset rather than included in the repository (~1.1 GB total). Data from utilized reference databases (miRBase, miRTarBase, GSEA Hallmark gene sets) are included in `2_reference-files/reference-dbs/`.

Before running the pipeline, download `GENCODEv48-ALTER-CURE-reference.tar.gz` from Releases and extract it into `2_reference-files/`:

```
tar -xzf GENCODEv48-ALTER-CURE-reference.tar.gz -C 1_rnaseq-final-analysis/2_reference-files/
```

The plasmid/reporter sequences used to build the combined reference are described in `GENCODEv48-ALTER-CURE.added.tsv` in the reference genome directory. The subdirectory `additions` includes the fasta and gtf files of the added elements. 

The base GENCODE files were:

| File Type |              GENCODE Content               | GENCODE Regions |                                                    Download Link                                                    |
| :-------: | :----------------------------------------: | :-------------: | :-----------------------------------------------------------------------------------------------------------------: |
|    gtf    |       Comprehensive gene annotation        |       PRI       | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.primary_assembly.annotation.gtf.gz |
|   fasta   | Genome sequence, primary assembly (GRCh38) |       PRI       |      https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/GRCh38.primary_assembly.genome.fa.gz      |
|   fasta   |            Transcript sequences            |       ALL       |         https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.transcripts.fa.gz          |

## Running the Pipeline on the Example Data

### Setup

1. Clone the repository and install Snakemake into a dedicated Conda environment.
2. Download `GENCODEv48-ALTER-CURE-reference.tar.gz` from Releases and extract into `2_reference-files/` (see Reference Data).
3. Open `1_rnaseq-final-analysis/config.yaml` and set `smk_dir` to the absolute path of the `1_rnaseq-final-analysis/` directory on your machine:

```yaml
smk_dir: '/path/to/ALTER-code/1_rnaseq-final-analysis'
```

All other config values (experiment names, reference genome name, analysis parameters) are pre-set to reproduce the paper results and do not need to be changed.

### Running

From within the `1_rnaseq-final-analysis/` directory, with the Snakemake Conda environment active:

```bash
snakemake --cores 1 --use-conda
```

Increase `--cores` to parallelize independent jobs. The `--use-conda` flag is required — without it Snakemake will ignore the environment definitions and use the active environment.

## Outputs

All outputs are written to `1_example-data/`. Paths below are relative to `1_rnaseq-final-analysis/`. `{experiment}` is `R098.2` or `CURE`; `{snp}` is `AG` or `CT`; `{condition}` refers to individual editor conditions as named in the sample map.

### Per-Experiment Variant Analysis

Located in `1_example-data/{1_R098.2 or 2_CURE}/4_variants/`.

| Output                                                    | Description                                                                                                                                                        |
| :-------------------------------------------------------- | :----------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `{experiment}.variants.annotated.parquet`                 | Full variant table annotated with coding strand, and per-sample allele frequencies. ref and alt sequences are changed to reflect the sequence in the coding strand |
| `{experiment}.variants.summary-stats.xlsx`                | Summary statistics per condition and SNP type                                                                                                                      |
| `{experiment}.variants.tgt-edit.tsv`                      | Variant table entries corresponding to the target edits of editors in the experiment                                                                               |
| `filter-counts-{snp}.tsv`                                 | Number of variants passing each filter at each step of off-target hit identification, per condition                                                                |
| `hit-id/{snp}/{condition}/`                               | Per-condition variant parquets at each filter stage (var, VOI, DP, GQ, non-wt)                                                                                     |
| `hit-id/{experiment}.overlap.{snp}.parquet`               | Final off-target hits passing all filters, with binary columns for each experimental condition indicating whether the entry was a hit for that condition           |
| `hit-id/{experiment}.overlap.{snp}.xlsx`                  | Tables showing how many off-target hits a given editor shares with other sets of editors.                                                                           |
| `hit-id/{experiment}.miRTar-freq.{snp}.xlsx`              | Tables of counts and frequencies for known miRNA target genes from miRTarBase appearing in off-target hits                                                         |
| `transcript-align/{experiment}.transcripts.{snp}.parquet` | Off-target hits aligned to transcript coordinates with sequence context                                                                                            |
| `transcript-align/{snp}/{condition}/`                     | Per-condition transcript-aligned hit parquets                                                                                                                      |
| `{experiment}.transcripts.{snp}.o10-consensus.xlsx`       | Sequence and structure context counts and frequencies in a ±10nt window centered on each off-target                                                                 |

### Per-Experiment Expression Analysis

Located in `1_example-data/{1_R098.2 or 2_CURE}/5_expression/deseq2/`.

| Output                                                                       | Description                                                                  |
| :--------------------------------------------------------------------------- | :--------------------------------------------------------------------------- |
| `2_txi-counts/{experiment}.raw-counts.gene-level.tsv.gz`                     | Gene-level raw counts aggregated by tximport from Salmon quantification      |
| `3_results/2_result-tables/1_norm-counts/`                                   | DESeq2 normalized count tables                                               |
| `3_results/2_result-tables/2_deg/{condition}.3_de-results.padj0.01-lfc1.tsv` | DEG results per condition (padj < 0.01, \|log2FC\| > 1)                      |
| `3_results/2_result-tables/3_gsea/{condition}.fgsea.tsv`                     | GSEA results per condition against Hallmark gene sets                        |
| `3_results/1_qc-plots/`                                                      | QC plots: geometric distance heatmap, Poisson distance heatmap, glmPCA (SVG) |
| `3_results/3_deg-plots/`                                                     | Volcano plots and DEG heatmaps per condition (SVG)                           |
| `3_results/4_go-plots/`                                                      | GO enrichment plots per condition (SVG)                                      |
| `3_results/{experiment}.degs.mirna-freq.xlsx`                                | Frequency of known miRNA targets in DEGs vs. expressed background            |
| `3_results/{experiment}.degs.off-tgt-freq.{snp}.xlsx`                        | Frequency of off-target genes in DEGs                                        |
| `deseq2.html`                                                                | Rendered DESeq2 analysis report                                              |

### Cross-Experiment Resources

Located in `1_example-data/3_analysis-resources/`. These files combine data from both experiments and are used in the aggregate results analyses.

| Output                                                                | Description                                                                                                  |
| :-------------------------------------------------------------------- | :----------------------------------------------------------------------------------------------------------- |
| `expr-transcripts/CURE.R098.2.transcripts.expressed.parquet`          | Expressed transcript table merged across R098.2 (ALTER) and CURE, used as the background for miRNA alignment |
| `transcript-overlap/CURE.R098.2.transcripts.{snp}.mirna.o100.parquet` | Off-target transcript table with miRNA seed alignment annotations, merged across both experiments            |

### Aggregate Results

Located in `1_example-data/4_aggregate_results/`.

| Output                                        | Description                                                                                                                                                               |
| :-------------------------------------------- | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `CURE.R098.2.mirna-align.{snp}.o100.xlsx`     | Total counts and frequencies of the number of miRNA target site alignments across the ±100 nucleotide context window for all off-target transcripts from both experiments |
| `CURE.R098.2.mirna-align.{snp}.o100.pos.xlsx` | Positional counts and frequencies of miRNA target site alignments across the ±100nt context window                                                                        |

## Initial Submission

Code from the initial submission is preserved in `3_initial-submission/` for transparency. This includes the original reference curation notebooks (`1_curating-refs/`), DEG analysis scripts (`2_deg-analysis/`), and off-target analysis notebooks (`3_off-tgt-analysis/`), along with rendered HTML example outputs (`4_example_notebooks/`).

---

# MARIA Analysis

## Dependencies

Unless otherwise indicated, code has only been tested with the dependency versions given in the following tables. Analysis code was run on a Microsoft Surface Pro 2 with Windows 10. Python dependencies were installed into dedicated Conda environments using Conda, including transitive dependencies. All installations could be completed in under half an hour. Running through the demo with the provided data may take 36 hours on a conventional laptop or as little as 4 hours on a dedicated computing cluster.

### Software

|Name|Version|
|:--:|:--:|
|R|4.3.2|
|Python|2.7.15|
|MARIA|2.20.2019|

### R Packages

|Name|Version|
|:--:|:--:|
|ggprism|1.0.5|
|ggplot2|3.5.1|
|progress|1.2.3|
|dplyr|1.1.4|
|data.table|1.16.2|

### Python Libraries

| Name | Version |
|:--:|:--:|
| certifi | 2019.11.28 |
| Cython | 0.29.15 |
| functools32 | 3.2.3.post2 |
| future | 0.18.2 |
| h5py | 2.9.0 |
| Keras | 2.0.3 |
| matplotlib | 2.2.3 |
| mkl-fft | 1.1.0 |
| mkl-service | 2.3.0 |
| pandas | 0.24.2 |
| pygpu | 0.7.6 |
| pyreadline | 2.1 |
| python-Levenshtein | 0.12.0 |
| tornado | 5.1.1 |
| unittest2 | 1.1.0 |
| wheel | 0.37.1 |
| wincertstore | 0.2 |

# Tutorial Workflow
A short description of the workflow can be found in `2_MARIA_ Deimmunization/1_examples/1_tutorial-workflow/README.md` 
Please see the accompanying tutorial document for step by step instructions in `Example De-Immunization Workflow.Rmd`
