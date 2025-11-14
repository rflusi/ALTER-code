
# Code Availability for Programed Manipulation of RNA Targets using Human Argonaute 2

This repository contains jupyter scripts, python scripts, and R scripts used in the RNAseq and MARIA analyses shown in the paper.

# General Notes

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

Jupyter and Rmd files are uploaded with all outputs cleared, example outputs are given in the examples directory as html files 

# RNAseq Analysis

## Reference Data

Genomic references were all based on GENCODE Human v48:

|File Type|GENCODE Content|GENCODE Regions|Download Link|
|:-------:|:-------------:|:-------------:|:-----------:|
|gtf|Comprehensive gene annotation|PRI|https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.primary_assembly.annotation.gtf.gz|
|fasta|Genome sequence, primary assembly (GRCh38)|PRI|https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/GRCh38.primary_assembly.genome.fa.gz|
|fasta|Transcript sequences|ALL|https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.transcripts.fa.gz|

The following references were also used at different points in the analysis:

|Database/Resource|File Type|Release|Download Link|
|:-:|:-:|:-:|:-:|
|dbsnp (NCBI)|vcf|138 (obtained from GATK resource bundle)|https://storage.cloud.google.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf|
|ClinVar (NCBI)|vcf|2025-07-15|https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar_20250715.vcf.gz|

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

## C-to-U Off-Target Analysis

### 1. Curating Reference Files

#### Matching the Transcripts Fasta to the Primary Assembly GTF

The downloaded transcripts fasta contains transcript entries that are not in the comprehensive primary assembly gtf, which can cause compatability issues with salmon for quantifying read counts. To fasta was filtered to match the gtf using the jupyter notebook 1_filter-transcript-fastagz-to-gtf.ipynb.

