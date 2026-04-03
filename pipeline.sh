#!/bin/bash
# ============================================================
# Crohn's Disease Metagenomics Analysis Pipeline
# Carlos Simon Foundation Technical Assignment
# Author: Larisa Atanasiu
# Date: April 2026
# Data: HMP2 study (Lloyd-Price et al., Nature 2019)
# ============================================================

# ============================================================
# TOOLS AND VERSIONS
# ============================================================
# FastQC 0.12.1       - Quality control
# MultiQC 1.33        - Aggregate QC report
# fastp 1.3.1         - Adapter trimming
# Bowtie2 2.5.5       - Host decontamination
# MetaPhlAn 3.0.13    - Taxonomic profiling
# HUMAnN 3.9          - Functional profiling
# R 4.3.2             - Statistical analysis and visualization

# ============================================================
# SAMPLES
# ============================================================
# Healthy controls: SRR5983438, SRR5983451, SRR5983484
# Crohn's disease:  SRR5983306, SRR5983464, SRR5983287
# Source: NCBI SRA, BioProject PRJNA398089 (HMP2/IBDMDB)

# ============================================================
# PROJECT SETUP
# ============================================================

mkdir -p ~/crohns_mgx/raw_data
mkdir -p ~/crohns_mgx/qc
mkdir -p ~/crohns_mgx/trimmed
mkdir -p ~/crohns_mgx/host_removed
mkdir -p ~/crohns_mgx/taxonomy
mkdir -p ~/crohns_mgx/functional
mkdir -p ~/crohns_mgx/figures
mkdir -p ~/crohns_mgx/host_db
mkdir -p ~/crohns_mgx/humann_db

# ============================================================
# STEP 1: DOWNLOAD RAW DATA
# Tool: SRA Toolkit 3.4.1
# ============================================================

cd ~/crohns_mgx/raw_data

# Download SRA files
prefetch SRR5983438 SRR5983451 SRR5983484 \
         SRR5983306 SRR5983464 SRR5983287

# Convert to paired-end FASTQ
fastq-dump --split-files --gzip \
    SRR5983438 SRR5983451 SRR5983484 \
    SRR5983306 SRR5983464 SRR5983287

# ============================================================
# STEP 2: QUALITY CONTROL
# Tool: FastQC 0.12.1 + MultiQC 1.33
# Purpose: Assess raw read quality, detect adapter contamination,
#          identify problematic samples before analysis
# ============================================================

# Run FastQC on all 12 FASTQ files (paired-end = 2 files per sample)
fastqc ~/crohns_mgx/raw_data/*.fastq.gz \
    -o ~/crohns_mgx/qc/ \
    --threads 4

# Aggregate all FastQC reports into one MultiQC report
multiqc ~/crohns_mgx/qc/ \
    -o ~/crohns_mgx/qc/multiqc_report

# Key QC findings:
# - SRR5983484: only 0.5M reads vs 7-17M in other samples (LOW DEPTH)
# - SRR5983306: 42% duplication rate (elevated but acceptable)
# - All samples: adapter contamination detected (expected, addressed by fastp)
# - Per base sequence content failures in early bases (normal for metagenomics)

# ============================================================
# STEP 3: ADAPTER TRIMMING AND QUALITY FILTERING
# Tool: fastp 1.3.1
# Purpose: Remove adapter sequences, low quality bases,
#          and reads shorter than 50bp after trimming
# Parameters:
#   --detect_adapter_for_pe: auto-detect paired-end adapters
#   --qualified_quality_phred 20: remove bases with Q<20
#   --length_required 50: discard reads shorter than 50bp
# ============================================================

cd ~/crohns_mgx

for sample in SRR5983287 SRR5983306 SRR5983438 SRR5983451 SRR5983464 SRR5983484; do
    echo "Trimming ${sample}..."
    fastp \
        -i raw_data/${sample}_1.fastq.gz \
        -I raw_data/${sample}_2.fastq.gz \
        -o trimmed/${sample}_1.fastq.gz \
        -O trimmed/${sample}_2.fastq.gz \
        --json trimmed/${sample}_fastp.json \
        --html trimmed/${sample}_fastp.html \
        --thread 4 \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --length_required 50 \
        2> trimmed/${sample}_fastp.log
    echo "Done: ${sample}"
done

# ============================================================
# STEP 4: HOST DECONTAMINATION
# Tool: Bowtie2 2.5.5
# Reference: GRCh38 (hg38) human genome
# Purpose: Remove human-derived reads to retain only microbial reads
# Parameters:
#   --very-fast: fast preset (appropriate for stool with low host contamination)
#   --un-conc-gz: output unmapped (microbial) paired reads
# ============================================================

# Download GRCh38 Bowtie2 index (pre-built)
cd ~/crohns_mgx/host_db
curl -L "https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip" -o GRCh38.zip
unzip GRCh38.zip

# Run host decontamination for all samples
cd ~/crohns_mgx

for sample in SRR5983287 SRR5983306 SRR5983438 SRR5983451 SRR5983464 SRR5983484; do
    echo "Removing host reads from ${sample}..."
    bowtie2 \
        -x host_db/GRCh38_noalt_as/GRCh38_noalt_as \
        -1 trimmed/${sample}_1.fastq.gz \
        -2 trimmed/${sample}_2.fastq.gz \
        --un-conc-gz host_removed/${sample}_%.fastq.gz \
        -p 2 \
        --very-fast \
        2> host_removed/${sample}_bowtie2.log \
        > /dev/null
    echo "Done: ${sample}"
done

# Summary of human read proportions
echo "=== Host contamination summary ==="
for sample in SRR5983306 SRR5983438 SRR5983451 SRR5983484 SRR5983464 SRR5983287; do
    echo -n "${sample}: "
    grep "overall alignment rate" host_removed/${sample}_bowtie2.log
done

# Key findings:
# SRR5983306 (Crohn):   0.74% human reads - PASS
# SRR5983438 (Control): 1.12% human reads - PASS
# SRR5983451 (Control): 0.04% human reads - PASS
# SRR5983484 (Control): 0.10% human reads - PASS (already flagged: low depth)
# SRR5983464 (Crohn):  35.06% human reads - FLAGGED (high host contamination)
# SRR5983287 (Crohn):   0.07% human reads - PASS

# ============================================================
# STEP 5: TAXONOMIC PROFILING
# Tool: MetaPhlAn 3.0.13
# Database: mpa_v31_CHOCOPhlAn_201901
# Purpose: Identify microbial species and their relative abundances
#          using clade-specific marker genes
# ============================================================

conda activate mgx

for sample in SRR5983306 SRR5983438 SRR5983451 SRR5983484 SRR5983464 SRR5983287; do
    echo "Running MetaPhlAn on ${sample}..."
    metaphlan \
        host_removed/${sample}_1.fastq.gz,host_removed/${sample}_2.fastq.gz \
        --input_type fastq \
        --nproc 4 \
        -o taxonomy/${sample}_metaphlan.txt \
        --bowtie2out taxonomy/${sample}_bowtie2.bz2 \
        2> taxonomy/${sample}_metaphlan.log
    echo "Done: ${sample}"
done

# Merge all sample outputs into one abundance table
merge_metaphlan_tables.py taxonomy/*_metaphlan.txt > taxonomy/merged_abundance.txt

# Extract species level only (remove strain level)
grep "s__" taxonomy/merged_abundance.txt | grep -v "t__" > taxonomy/species_abundance.txt

# Extract genus level
grep "g__" taxonomy/merged_abundance.txt | grep -v "s__" > taxonomy/genus_abundance.txt

# Extract phylum level
grep "p__" taxonomy/merged_abundance.txt | grep -v "c__" > taxonomy/phylum_abundance.txt

echo "Taxonomic profiling complete"
echo "Species detected: $(wc -l < taxonomy/species_abundance.txt)"

# ============================================================
# STEP 6: FUNCTIONAL PROFILING
# Tool: HUMAnN 3.9
# Databases: ChocoPhlAn (nucleotide), UniRef90 (protein)
# Purpose: Quantify metabolic pathway abundances
#          to understand functional potential of the microbiome
# ============================================================

# Configure HUMAnN databases
humann_config --update database_folders nucleotide ~/crohns_mgx/humann_db/chocophlan
humann_config --update database_folders protein ~/crohns_mgx/humann_db/uniref

# Run HUMAnN on each sample
for sample in SRR5983306 SRR5983438 SRR5983451 SRR5983484 SRR5983464 SRR5983287; do
    echo "Running HUMAnN on ${sample}..."

    # Concatenate paired reads for HUMAnN input
    cat host_removed/${sample}_1.fastq.gz \
        host_removed/${sample}_2.fastq.gz \
        > functional/${sample}_combined.fastq.gz

    humann \
        --input functional/${sample}_combined.fastq.gz \
        --output functional/${sample}_humann \
        --metaphlan-options "--bowtie2db ~/crohns_mgx/taxonomy" \
        --threads 4 \
        2> functional/${sample}_humann.log

    echo "Done: ${sample}"
done

# Join pathway abundance tables
humann_join_tables \
    --input functional/ \
    --output functional/merged_pathabundance.tsv \
    --file_name pathabundance

# Normalize to relative abundance
humann_renorm_table \
    --input functional/merged_pathabundance.tsv \
    --output functional/merged_pathabundance_relab.tsv \
    --units relab

echo "Functional profiling complete"

# ============================================================
# DOWNSTREAM ANALYSIS IN R
# See analysis.R for:
# - Alpha diversity (Shannon, Richness)
# - Beta diversity (Bray-Curtis PCoA, PERMANOVA)
# - Taxonomic composition (stacked barplots)
# - Faecalibacterium prausnitzii analysis
# - Differential abundance (MaAsLin2)
# - Functional pathway analysis
# ============================================================

echo "Pipeline complete. Run analysis.R for statistical analysis."
