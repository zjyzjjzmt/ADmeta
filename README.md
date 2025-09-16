

# Anaerobic Digestion Metagenomics Analysis Pipeline

This repository contains a modular bioinformatics pipeline for analyzing metagenomic data from anaerobic digestion of livestock manure. The pipeline processes raw sequencing reads, performs quality control, assembles contigs, identifies viral communities, and analyzes metagenome-assembled genomes (MAGs) to investigate antimicrobial resistance genes (ARGs), microbial community structures, phage-host interactions, and antiviral defense systems. The pipeline is organized into four independent yet interconnected Bash scripts, each addressing a specific stage of analysis.

## Table of Contents
- [Overview](#overview)
- [Pipeline Modules](#pipeline-modules)
  - [01Reads_module.sh](#01reads_modulesh)
  - [02Contigs_module.sh](#02contigs_modulesh)
  - [03Virus_module.sh](#03virus_modulesh)
  - [04MAGs_module.sh](#04mags_modulesh)


## Overview
This pipeline is designed to process high-throughput metagenomic sequencing data from anaerobic digestion systems, providing insights into microbial and viral communities, ARGs, and phage-host interactions. It is tailored for researchers studying microbial ecology and antibiotic resistance in livestock manure anaerobic digestion environments. The pipeline is divided into four modules, each handling a distinct stage of analysis, from read preprocessing to MAG-based functional analyses.

## Pipeline Modules

### 01Reads_module.sh
**Purpose**: Performs quality control (QC) on raw metagenomic reads, detects ARGs at the read level, and profiles microbial community structures.

**Key Functions**:
- **Quality Control**: Uses MetaWrap to trim adapters and remove low-quality reads, ensuring high-quality input for downstream analyses.
- **ARG Detection**: Identifies ARGs using SARG and CARD databases via Diamond BLASTX, summarizing resistance gene families across samples.
- **Microbial Community Profiling**: Uses Kraken2 (16S rRNA) and MetaPhlAn4 to analyze microbial community composition, followed by Graphlan for visualization of taxonomic profiles and diversity metrics.

**Use Case**: Provides an initial assessment of raw data quality, ARG prevalence, and microbial diversity, setting the stage for contig assembly and deeper analyses.

### 02Contigs_module.sh
**Purpose**: Assembles reads into contigs, performs binning to recover MAGs, and conducts contig-level ARG analysis.

**Key Functions**:
- **Assembly**: Uses MetaWrap (MEGAHIT-based) to assemble clean reads into contigs.
- **Binning and Refinement**: Applies MetaBat2, CONCOCT, and MaxBin2 for binning, followed by MetaWrap bin refinement to generate high-quality MAGs.
- **ARG Detection**: Identifies ARGs in contigs using SARG, CARD, ICEs, and MGEs databases, with CD-HIT for deduplication.
- **Taxonomic Classification**: Uses CAT and PlasFlow to classify contigs and identify plasmid sequences.
- **Abundance Estimation**: Maps reads back to ARGs-containing contigs using CoverM to quantify abundance (RPKM, counts, and covered bases).

**Use Case**: Generates assembled contigs and preliminary MAGs, enabling ARG detection at the contig level and providing input for viral and MAG analyses.

### 03Virus_module.sh
**Purpose**: Identifies and characterizes viral communities from assembled contigs, including ARG detection and taxonomic classification.

**Key Functions**:
- **Viral Detection**: Uses VirFinder, VirSorter2, and geNomad to identify viral sequences (>5000 bp) from contigs.
- **Quality Control**: Validates viral sequences with CheckV to filter high-quality viral contigs.
- **Clustering**: Performs viral sequence clustering using BLAST and ANI-based clustering (anicalc.py, aniclust.py) to identify viral clusters (VCs).
- **ARG Detection**: Identifies ARGs in viral sequences using SARG and CARD databases.
- **Taxonomic and Functional Analysis**: Uses CAT, vConTACT2, RefSeq, and IMG/VR for taxonomic classification and novelty assessment.
- **Abundance Estimation**: Maps reads to viral sequences using CoverM to quantify viral abundance.
- **Functional Annotation**: Uses VIBRANT to annotate viral functions.

**Use Case**: Mines viral communities, assesses their role in ARG dissemination, and provides viral sequences for phage-host interaction analyses.

### 04MAGs_module.sh
**Purpose**: Analyzes MAGs for ARGs, phage-host interactions, and antiviral defense systems.

**Key Functions**:
- **MAG Dereplication**: Uses dRep to deduplicate MAGs with high completeness (>70%) and low contamination (<10%).
- **Taxonomic Classification**: Uses GTDB-Tk to classify MAGs taxonomically.
- **ARG Detection**: Identifies ARGs in MAGs using SARG, CARD, BacMet, ICEs, and Victors databases, with MGEs analysis via BLAST.
- **Phage-Host Linkage**: Detects CRISPR spacers (using CRT) and tRNA sequences (using Aragorn) to infer phage-host interactions by mapping to viral sequences from `03Virus_module.sh`.
- **Antiviral Defense Systems**: Identifies defense systems in MAGs using PADLOC and DefenseFinder.
- **Homology Analysis**: Performs BLAST-based homology analysis between MAGs and viral sequences to identify shared genomic regions.

**Use Case**: Provides detailed insights into microbial genomes, their ARGs, interactions with phages, and defense mechanisms against viral infections.

## Dependencies
The pipeline requires the following tools and databases, installable via Conda or manually:

- **Quality Control and Reads Analysis**:
  - MetaWrap (v1.3.2+), FastQC, Trimmomatic
  - Kraken2 (v2.1.3+), MetaPhlAn4 (v4.1+), Graphlan
- **Assembly and Binning**:
  - MetaWrap (MEGAHIT, MetaBat2, CONCOCT, MaxBin2)
- **Viral Analysis**:
  - VirFinder, VirSorter2 (v2.2.4+), geNomad (v1.7+), CheckV (v1.0.3+), VIBRANT, vConTACT2
- **MAG Analysis**:
  - dRep (v3.4+), GTDB-Tk (v2.3.2+), Prokka (v1.14+), PADLOC (v2.0+), DefenseFinder (v1.3+)
- **ARG Detection**:
  - Diamond (v2.1.8+), BLAST+ (v2.15+)
  - Databases: SARG (v3.0+), CARD (v3.2+), BacMet (v2.0+), ICEberg (v3.0+), Victors, MGEs
- **Phage-Host Analysis**:
  - Aragorn (v1.2.36), CRT (v1.2; consider CRISPRCasFinder v2.1+)
- **General Utilities**:
  - seqtk, cd-hit-est, taxonkit, csvtk, datamash, pigz

## Contact
For questions or support, please contact:
- **Maintainer**: Junya Zhang
- **Email**: jyzhang@rcees.ac.cn
