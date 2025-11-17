# Provided Data

All provided data for this workflow is located in the ALTER-code/5_tutorial-workflows/1_degs-off-tgts directory, where this file should also be located:
|Data|Description|
|:--:|:----------|
|all-samples-variants.tsv.gz|final variant table for the ALTER-4 data from the HapplotypeCaller workflow described in the Supplementary Information|
|salmon-results|directory containing raw count data from the HapplotypeCaller workflow described in the Supplementary Information|
|sample-map.tsv|a tsv file with sample information 'sample', 'condition','replicate'| 

# Data not provided

Genomic reference data and the different reference databases are not provided. This workflow assumes that genomic references with the custom plasmid contigs are being used (see the main READ_ME for more info). It should be workable with standard GENCODE references but that has not been explicitly tested.

# Curate Reference Data

**NOTE:** Even if you use references without custom contigs it is still a good idea to standardize the transcripts fasta to the gtf. The transcript map tsv produced by step 2, and the ClinVar vcf with ENCODE compatable chromosome notation produced by step 3 will also be required in other parts of the workflow.

## 1. Filter the Transcript fasta

This should be done prior to building the custom references:

1. open ALTER-code/1_curating-refs/1_filter-transcript-fastagz-to-gtf.ipynb
2. complete cell 1 line 9 with the path to the reference transcripts fasta file.
3. complete cell 1 line 10 with the path to the reference gtf
4. run the notebook

Output: a new reference transcripts fasta file containing only transcripts found in the reference gtf

## 2. Generate Transcript Map

This should be done after all plasmid information has been added to the reference:

1. open ALTER-code/1_curating-refs/2_map-transcripts-to-gene.ipynb
2. complete cell 1 line 9 with the path to the reference gtf
3. complete cell 1 line 10 with the path to the reference transcripts fasta file.
4. run the notebook

Output: a tsv with one entry per transcript id in the reference data. Entries also include information on the corresponding gene. 

## 3. Change ClinVar Chromosome Notation

Make the chromosome notation in the clinvar vcf downloaded from NCBI match the notation used in the GENCODE references:

1. open ALTER-code/1_curating-refs/3_change-clinvar-contigs.sh
2. complete line 9 with the full path to the directory containing the vcf
3. complete line 12 with the name of the vcf file not including the file extension
4. complere line 15 with the full path to ALTER-code/1_curating-refs/3_clinvar-contig-map.txt
5. open a terminal
6. run the command

```
bash <path to ALTER-code/1_curating-refs/3_change-clinvar-contigs.sh>
```

Output: A new vcf with the corrected notation

# DEG Analysis

## 4. Pad Salmon Results

Salmon analysis was completed for each sample with a reference genome including only the plasmid relevant to that sample. To standardize the contigs of all the salmon outputs to include all plasmids:

1. open ALTER-code/2_deg-analysis/1_pad-salmon.py
2. complete line 157 with the full path to ALTER-code/5_tutorial-workflows/1_degs-off-tgts
3. complete line 159 with the full path to the transcript map tsv created in step 2

Outputs:

1. ALTER-code/5_tutorial-workflows/1_degs-off-tgts/combined-results which contains aggregated salmon results
2. ALTER-code/5_tutorial-workflows/1_degs-off-tgts/padded-results which contains one subdirectory for each sample with padded results for that sample
3. ALTER-code/5_tutorial-workflows/1_degs-off-tgts/combined-results/deseq2 which will be used as the output directory for the upcoming DESeq2 analysis

## 5. DEG Analysis with DESeq2

The DEG analysis can be completed as follows:

1. open ALTER-code/2_deg-analysis/5_deseq2.Rmd
2. complete line 48 with the full path to ALTER-code/5_tutorial-workflows/1_degs-off-tgts
3. complete line 50 with the full path to the transcript map
4. complete line 52 with the full path to the GSEA Hallmark Gene Sets (by symbol) gmt file
5. run the rmd file

Outputs:

1. raw count tables aggregated to the gene level in ALTER-code/5_tutorial-workflows/1_degs-off-tgts/salmon-results/deseq2/2_txi-counts
2. In deseq2/3_results/1_qc-plots various RNAseq qc plots
3. In deseq2/3_results/2_result-tables:
    1. 02_ALTER.4.PPIB.14L8_de-results.tsv - differential expression results
    2. 02_ALTER.4.PPIB.14L8_de-results-drop-na.tsv - differential expression results with genes that could not be assayed dropped
    3. 02_ALTER.4.PPIB.14L8_de-results-padj0.01-lfc1.tsv - differential expression results which surpassed the log2FC and padj cutoffs (DEGs)
    4. norm-gene-counts.tsv - normalized gene level transcript counts
    5. norm-gene-pseudocounts.tsv - normalized counts incremented by 1 for plotting on log scale
4. In 3_deg-plots volcano plots and heat map of differential expression in the ALTER-4 sample
5. In 4_go-plots GO term enrichment analysis dot plots

## 6. Rate of Known miRNA Targets in DEGs

1. open ALTER-code/2_deg-analysis/6_miRNA-tgts-in-degs.ipynb
2. complete line 6 with the full path to ALTER-code/5_tutorial-workflows/1_degs-off-tgts
3. complete line 8 with the full path to the miRtarBase csv
4. complete line 10 with the full path to the transcript map
5. run the notebook

Outputs:

1. ALTER-code/5_tutorial-workflows/1_degs-off-tgts/salmon-results/deseq2/4_mirna-tgt-freq containing data on the rate of miRNA targets in the DEG set

# C-to-U Off-Target Analysis

## 7. Annotate the variant calls tsv

Annotate the table of variant calls from the GATK results and prepare it for downstream analysis. Infer strand and adjust sequences, add per sample base data for reference and alt bases, convert GT calls to binary, cross reference variant entries to ClinVar, calculate %_ref for all samples in all entries, and calculate %_snp if the target SNP is detected in a sample in a given entry.

1. open ALTER-code/3_off-tgt-analysis/7_var-table-annotation.ipynb
2. complete line 8 with the full path to ALTER-code/5_tutorial-workflows/1_degs-off-tgts/all-samples-variants.tsv.gz
3. complete line 10 with the full path to the reference gtf
4. complete like 12 with the full path the the ClinVar vcf
5. run the notebook

Outputs: ALTER-code/5_tutorial-workflows/1_degs-off-tgts/init-processing

1. all-samples-variants-proc.tsv.gz - the annotated full variant table
2. all-samples-variants-summary-stats.xlsx - descriptive data on the variants and annotations

## 8. Call Off-Target Hits

Filter annotated variant table for final off-tgt hits. Off target hits are defined as passing the following filters, strategy similar to Huang, X. et al. Programmable C-to-U RNA editing using the human APOBEC3A deaminase. The EMBO Journal 39, e104741 (2020):

1. var    - sample was called as a variant by HapplotypeCaller
2. VOI    - sample variant matches the variant of interest
3. DP     - read depth >= 20
4. GQ     - depth by quality >=20 (corresponds to 99% confidence)
5. non-wt - pct_ref in the matched wt sample is >99%, meaning the SNP was introduced by editing 

To complete the analysis:

1. open ALTER-code/3_off-tgt-analysis/8_off-tgt-hit-id.ipynb
2. complete cell 2 line 8 with the full path to ALTER-code/5_tutorial-workflows/1_degs-off-tgts/all-samples-variants.tsv.gz
3. complete cell 2 line 10 with the full path to the reference gtf
4. run the notebook

Outputs: ALTER-code/5_tutorial-workflows/1_degs-off-tgts/init-processing/off-tgt-analysis/3-reps

1. a subdirectory for each condition containing tsvs of entries that passed each filter
2. a tsv containing a count matrix of number of entries that passed each filter
3. a tsv containing identity and %_snp data for each off-tgt hit

## 9. Analysis of Off-Target Sequence Context

Analyze the sequence and RNA secodary structural contexts of final off target hits. Also perform analysis of guide alignment proximal to the off-target hits.

1. open ALTER-code/4_example_notebooks/3_off-tgt-analysis/9_off-tgt-seq-analysis.html
2. complete cell 2 line 8 with the full path to ALTER-code/5_tutorial-workflows/1_degs-off-tgts/all-samples-variants.tsv.gz
3. complete cell 2 line 10 with the full path to the reference gtf
4. complete cell 2 line 12 with the full path to the reference genome fasta
5. complete cell 2 line 14 with the full path to the reference genome transcripts fasta
6. run the notebook

Outputs: ALTER-code/5_tutorial-workflows/1_degs-off-tgts/init-processing/off-tgt-analysis/3-reps/seq-analysis

1. off-tgt-seq-counts.xlsx - seequence and structure consensus counts and frequencies
2. alignment-counts.xlsx - guide alignment counts and frequencies

## 10. Rate of Off-Targets in DEGs

1. open ALTER-code/3_off-tgt-analysis/10_DEG-edit-rate.py
2. complete line 104 with the full path to ALTER-code/4_example_notebooks/3_off-tgt-analysis/9_off-tgt-seq-analysis.html
3. Run the script

Output: ALTER-code/5_tutorial-workflows/1_degs-off-tgts/init-processing/off-tgt-analysis/3-reps/deg-analysis

1. 02_ALTER.4.PPIB.14L8 - contains counts and rates of Off-Targets in DEGs
