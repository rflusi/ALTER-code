
# Code Availability for Programed Manipulation of RNA Targets using Human Argonaute 2

This repository contains jupyter notebooks, python scripts, and R scripts and rmd files used in the RNAseq and MARIA analyses shown in the paper.

# General Notes

Jupyter and Rmd files are uploaded with all outputs cleared, example outputs are given in the 4_example-notebooks directory as html files but with the directory variables cleared.

# RNAseq Analysis

## Dependencies

Unless otherwise indicated, code has only been tested with the dependency versions given in the following tables. For RNAseq analysis the code was run on an M4 MacBook Pro with the macOS Sequoia. Dependencies were installed into dedicated Conda environments using Conda, including transitive dependencies. All installations could be completed in under half an hour. Running through the demo with the provided data may take 1 hour.

### Software

|Name|Version|
|:--:|:--:|
|python|3.12|
|ipykernel|7.1.0|
|R|4.3|
|samtools|1.16.1|
|bcftools|1.16|

### Python Libraries

|Name|Version|
|:--:|:--:|
|tqdm|4.67.1|
|biopython|1.85|
|pysam|0.23.0|
|pandas|2.2.3|
|numpy|1.26.4|

### R Packages

|Name|Version|
|:--:|:--:|
|DESeq2|1.46.0|
|tximport|1.34.0|
|dplyr|1.1.4|
|ggplot2|3.5.2|
|pheatmap|1.0.13|
|RColorBrewer|1.1.3|
|PoiClaClu|1.0.2.1|
|glmpca|0.2.0|
|apeglm|1.28.0|
|ggrepel|0.9.6|
|dplyr|1.1.4|
|tidyr|1.3.1|
|purrr|1.1.0|
|svglite|2.2.1|
|fgsea|1.32.2|
|data.table|1.17.8|
|msigdbr|25.1.1|
|org.Hs.eg.db|3.20.0|
|clusterProfiler|4.14.0|
|AnnotationDbi|1.68.0|
|enrichplot|1.26.1|
|biomaRt|2.62.0|
|ggrepel|0.9.6|

## Reference Data

Genomic references were all based on GENCODE Human v48:

|File Type|GENCODE Content|GENCODE Regions|Download Link|
|:-------:|:-------------:|:-------------:|:-----------:|
|gtf|Comprehensive gene annotation|PRI|https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.primary_assembly.annotation.gtf.gz|
|fasta|Genome sequence, primary assembly (GRCh38)|PRI|https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/GRCh38.primary_assembly.genome.fa.gz|
|fasta|Transcript sequences|ALL|https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.transcripts.fa.gz|

For each sample transfected with plamid(s) the reference genome fasta, gtf, and transcript fasta for that sample were updated with the plasmid(s) sequences and information. Files with added data are located in ALTER-code/1_curating-refs/custom-ref-data

|Genome Reference|Added to Genome Fasta|Added to gtf|Added to Transcript Fasta|
|:-------:|:-------------:|:-------------:|:-------------:|
|ALTER-transfection-control|pRFL462.fa|pRFL462.gtf|pRFL462-transcript.fa|
|ALTER-4|pRFL382.fa|pRFL382.gtf|pRFL382-transcript.fa|
|ALTER-master-ref|pRFL382-462.fa|pRFL382-462.gtf|pRFL382-462-transcript.fa|
|A3A-overexpression and CURE "wtHEK"|CURE-reporter.fa|CURE-reporter.gtf|CURE-reporter-transcripts.fa|
|CURE-C2|CURE-C2.fa|CURE-C2.gtf|CURE-C2-transcript.fa|
|CURE-N|CURE-N.fa|CURE-N.gtf|CURE-N-transcript.fa|
|RESCUE-S|RESCUE-S.fa|RESCUE-S.gtf|RESCUE-S-transcript.fa|
|CURE-master-ref|CURE-all.fa|CURE-all.gtf|CURE-all-transcript.fa|

The following references were also used at different points in the analysis:

|Database/Resource|File Type|Release|Download Link|
|:-:|:-:|:-:|:-:|
|dbsnp (NCBI)|vcf|138 (obtained from GATK resource bundle)|https://storage.cloud.google.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf|
|ClinVar (NCBI)|vcf|2025-07-15|https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar_20250715.vcf.gz|
|Hallmark GSEA Gene Sets|gmt|v2025.1.Hs.symbols|https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2025.1.Hs/h.all.v2025.1.Hs.symbols.gmt
|miRtarBase|csv|Homo sapiens|https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2025/cache/download/10.0/hsa_MTI.csv

Unless otherwise indicated, all reference data files were compressed with:

```
bgzip -l 9
```

and indexed with one of:

```
samtools faidx
tabix -p gff
tabix -p vcf
```

## 1. Curating Reference Files

All code for this section is found in the 1_curating-refs directory

### Matching the Transcripts Fasta to the Primary Assembly GTF

The downloaded transcripts fasta contains transcript entries that are not in the comprehensive primary assembly gtf, which can cause compatability issues with salmon for quantifying read counts. The fasta was filtered to match the gtf using the jupyter notebook 1_filter-transcript-fastagz-to-gtf.ipynb.

tsv files which map transcript ids to gene ids and gene names can then be produced with 2_map-transcripts-to-gene.ipynb

### Changing ClinVar vcf Contigs

The vcf file for ClinVar uses different notation for chromosomes than GENCODE does, which leads to compatability issues. The contigs of the vcf file can be updated with 3_change-clinvar-contigs.sh which requires 3_clinvar-contig-map.txt, which maps the old contigs to the new.

## 2. DEG Analysis

### DESeq2 Analysis

Raw count data from salmon was padded with all transcripts in the master reference so that all samples would share a reference, 4_pad-salmon.py.

DESeq2 was used to aggregate raw counts to the gene level, normalize counts, and perform DEG analysis. A ranked list of assayed genes was used in GSEA analysis, 5_deseq2.Rmd.

## Analysis of Known miRNA Targets in DEGs

The final DEG hits were cross referenced with miRtarBase to assess the rate of known miRNA targets in DEGs, 6_miRNA-tgts-in-degs

## 3. C-to-U Off-Target Analysis

All code for this section is found in the 2_off-tgt-analysis directory. The starting point for the code provided here is the final variant table obtained after analysis of the raw read data. For details on read qc, alignment, and variant calling please see the Supplementary Information of the paper.

### Annotation of Variant Table from GATK Analysis

The final variant table from GATK is reformatted, and annotated with ClinVar data for variants found in ClinVar. The coding strand of each variant is inferred and sequences are altered to reflect the coding strand as required. Percentages of reads matching the reference nucleotide are calculated for all samples in all entries; percentages of reads matching the alternate nucleotide are calculated for samples matching the target SNP (C-to-U in this case). These analyses were performed with 7_var-table-annotation.ipynb

### Identification of Off-Target Hits

Off-target hits were defined as those variant calls meeting the following criteria, following the precedent of Huang, X. et al. Programmable C-to-U RNA editing using the human APOBEC3A deaminase. The EMBO Journal 39, e104741 (2020).

1. var    - sample was called as a variant by HapplotypeCaller
2. VOI    - sample variant matches the variant of interest
3. DP     - read depth >= 20
4. GQ     - depth by quality >=20 (corresponds to 99% confidence)
5. non-wt - pct_ref in the matched wt sample is >99%, meaning the SNP was introduced by editing

This analysis was performed with 8_off-tgt-hit-id.

## Sequence-Context Analysis of Off-Target Hits

For each of the final off-target hits, the transcript sequence context in a 21 base window centered on the off-target was pulled from genomic data. The lowest mfe secondary structure was calculated for each context window. Sequence and structure counts at each position in the window were tabulated. The context window was widened to 51 and the contexts were analyzed for guide sequence alignment. 9_off-tgt-seq-analysis.ipynb.

## Rate of Off-Targets in DEGs

The frequency of off-targets in DEGs was calculated with 10_DEG-edit-rate.py.

# Tutorial Workflow

See ALTER-code/5_tutorial-workflows/1_degs-off-tgts/TUTORIAL_READ_ME.md

# MARIA Analysis

## Dependencies

Unless otherwise indicated, code has only been tested with the dependency versions given in the following tables. Analysis code was run on a Microsoft Surface Pro 2 with Windows 10. Python dependencies were installed into dedicated Conda environments using Conda, including transitive dependencies. All installations could be completed in under half an hour. Running through the demo with the provided data may take 48 hours on a convential laptop or as little as 4 hours on a dedicated computing cluster.

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
A short description of the workflow can be found in 5_tutorial-workflows/2_MARIA_Pipeline/README.md
Please see the accompanying tutorial document for step by step instructions in 5_tutorial-workflows/2_MARIA_Pipeline/Example De-Immunization Workflow.Rmd
