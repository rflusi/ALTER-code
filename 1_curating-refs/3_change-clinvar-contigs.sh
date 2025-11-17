#!/bin/bash
#run as script

# --- load modules ---

# --- User-defined variables ---

# full path to directory where the vcfgz is located
PROJ_DIR=""

# name of the clinvar vcf, without file extension. It should be compressed with bgzip (vcf.gz) and indexed
INPUT_VCF_NAME=""

# full path to tab sepparated text file mapping old and new contigs
CONTIG_MAP=""

# --- Automatically defined variables ---
INPUT_VCFGZ="$PROJ_DIR/$INPUT_VCF_NAME.vcf.gz"
OUT_VCFGZ="$PROJ_DIR/$INPUT_VCF_NAME.new_contigs.vcf.gz"

# --- Step 1: change vcf sample ---
echo "[RUN] change $INPUT_VCF_NAME contigs"

bcftools annotate --rename-chrs "$CONTIG_MAP" -Oz -o "$OUT_VCFGZ" "$INPUT_VCFGZ"
tabix -p vcf "$OUT_VCFGZ"

echo "[DONE] change $INPUT_VCF_NAME contigs"