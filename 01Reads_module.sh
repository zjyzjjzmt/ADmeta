#!/bin/bash
#SBATCH --job-name=AD_reads_analysis
#SBATCH --output=AD_reads_analysis-%j.log
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=50
#SBATCH --mem=200G

# Set strict error handling
set -euo pipefail

# ---------------------- Path Settings ----------------------
WORKDIR="/media/kunyu/data2/03AD"
RAWREADS="${WORKDIR}/01RAW_READS"
QC_DIR="${WORKDIR}/01READ_QC"
CLEANREADS="${WORKDIR}/01Reads_fq"
OUTDIR="${WORKDIR}/04reads_based_analysis"
SARG_DIR="${OUTDIR}/01SARG"
CARD_DIR="${OUTDIR}/02card"
KRAKEN_DIR="${OUTDIR}/07kraken2"
METAPHLAN_DIR="${WORKDIR}/10metaphlan4"
TMP_DIR="${WORKDIR}/tmp"

# ---------------------- Utility Functions ----------------------
# Activate Conda environment
activate_conda() {
    local env=$1
    conda deactivate >/dev/null 2>&1 || true
    conda activate "$env" || { echo "Error: Failed to activate Conda env: $env"; exit 1; }
}

# Check if file exists
check_file() {
    local file=$1
    if [[ ! -f "$file" ]]; then
        echo "Error: File $file not found"
        exit 1
    fi
}

# Ensure directories exist
ensure_dir() {
    for dir in "$@"; do
        mkdir -p "$dir" || { echo "Error: Failed to create directory: $dir"; exit 1; }
    done
}

# ---------------------- Initialize ----------------------
ensure_dir "$QC_DIR" "$CLEANREADS" "$SARG_DIR" "$CARD_DIR" "$KRAKEN_DIR" "$METAPHLAN_DIR" "$TMP_DIR"
activate_conda base

# 1. Read Quality Control
echo "Running quality control on raw reads..."
for F in "${RAWREADS}"/*_1.fastq; do
    [[ -f "$F" ]] || continue
    R="${F%_*}_2.fastq"
    check_file "$F"
    check_file "$R"
    SAMPLE=$(basename "${F%_*}")

    echo "Processing sample: $SAMPLE"
    metawrap read_qc -1 "$F" -2 "$R" -t 20 --skip-bmtagger -o "${QC_DIR}/${SAMPLE}"

    # Compress raw reads
    pigz -p 10 "$F" "$R" || echo "Warning: Failed to compress $F or $R"
done

# Move clean reads
echo "Moving cleaned reads to $CLEANREADS..."
for dir in "${QC_DIR}"/*; do
    [[ -d "$dir" ]] || continue
    b=$(basename "$dir")
    mv "${dir}/final_pure_reads_1.fastq" "${CLEANREADS}/${b}_1.fastq" || echo "Warning: Failed to move reads_1 for $b"
    mv "${dir}/final_pure_reads_2.fastq" "${CLEANREADS}/${b}_2.fastq" || echo "Warning: Failed to move reads_2 for $b"
done

# Count clean reads
echo "Counting clean reads..."
activate_conda tools
: > "${CLEANREADS}/AD_Meta_reads_number.txt"
for file in "${CLEANREADS}"/*_1.fastq; do
    [[ -f "$file" ]] || continue
    count=$(grep -c "^@" "$file")
    echo "$(basename "$file" _1.fastq): $count reads" >> "${CLEANREADS}/AD_Meta_reads_number.txt"
done

# 2. SARG Analysis
echo "Running SARG analysis..."
check_file "${SARG_DIR}/SARG.2.2.fasta"
diamond makedb --in "${SARG_DIR}/SARG.2.2.fasta" -d "${SARG_DIR}/SARG.2.2_nr" --threads 10

for fq in "${CLEANREADS}"/*_1.fastq; do
    [[ -f "$fq" ]] || continue
    out="${SARG_DIR}/$(basename "$fq" _1.fastq)-SARG.txt"
    echo "Running Diamond BLASTX for SARG on $(basename "$fq")"
    diamond blastx -d "${SARG_DIR}/SARG.2.2_nr" -q "$fq" -o "$out" \
        --evalue 1e-5 --query-cover 75 --id 90 -k 1 --threads 10
done

# Summarize SARG results
echo "Summarizing SARG results..."
echo -e "SARG_family\t$(ls "${SARG_DIR}"/*-SARG.txt | sed 's|.*/||;s|-SARG.txt||;s|^|count |' | tr '\n' '\t')" \
    > "${SARG_DIR}/AD_SARG.xls"

while IFS= read -r i; do
    [[ -n "$i" ]] || continue
    echo -ne "$i\t" >> "${SARG_DIR}/AD_SARG.xls"
    for file in "${SARG_DIR}"/*-SARG.txt; do
        count=$(grep -c -F "$i" "$file")
        echo -ne "$count\t" >> "${SARG_DIR}/AD_SARG.xls"
    done
    echo >> "${SARG_DIR}/AD_SARG.xls"
done < <(cut -f 1 "${SARG_DIR}/SARG.2.2.txt" | sort -u)

awk -F'\t' 'BEGIN {OFS="\t"} NR==1 {print} NR>1 {sum=0; for (i=2; i<=NF; i++) sum+=$i; if (sum!=0) print}' \
    "${SARG_DIR}/AD_SARG.xls" > "${SARG_DIR}/AD_SARG_final.xls"

csvtk join -t --left-join -f 1 "${SARG_DIR}/AD_SARG_final.xls" "${SARG_DIR}/SARG_mapping-20220906.txt" \
    > "${SARG_DIR}/AD_SARG_mapping.xls"

# 3. CARD Analysis
echo "Running CARD analysis..."
check_file "${CARD_DIR}/CARD3.1.4.fasta"
diamond makedb --in "${CARD_DIR}/CARD3.1.4.fasta" -d "${CARD_DIR}/card3.1.4_nr" --threads 10

for fq in "${CLEANREADS}"/*_1.fastq; do
    [[ -f "$fq" ]] || continue
    out="${CARD_DIR}/$(basename "$fq" _1.fastq)-CARD.txt"
    echo "Running Diamond BLASTX for CARD on $(basename "$fq")"
    diamond blastx -d "${CARD_DIR}/card3.1.4_nr" -q "$fq" -o "$out" \
        --evalue 1e-5 --query-cover 75 --id 90 -k 1 --threads 10
done

# Summarize CARD results
echo "Summarizing CARD results..."
echo -e "CARD_family\t$(ls "${CARD_DIR}"/*-CARD.txt | sed 's|.*/||;s|-CARD.txt||;s|^|count |' | tr '\n' '\t')" \
    > "${CARD_DIR}/AD_CARD.xls"

while IFS= read -r i; do
    [[ -n "$i" ]] || continue
    echo -ne "$i\t" >> "${CARD_DIR}/AD_CARD.xls"
    for file in "${CARD_DIR}"/*-CARD.txt; do
        count=$(grep -c -F "$i" "$file")
        echo -ne "$count\t" >> "${CARD_DIR}/AD_CARD.xls"
    done
    echo >> "${CARD_DIR}/AD_CARD.xls"
done < <(cut -f 3 "${CARD_DIR}/CARD3.1.4.txt" | sort -u)

awk -F'\t' 'BEGIN {OFS="\t"} NR==1 {print} NR>1 {sum=0; for (i=2; i<=NF; i++) sum+=$i; if (sum!=0) print}' \
    "${CARD_DIR}/AD_CARD.xls" > "${CARD_DIR}/AD_CARD_final.xls"

csvtk join -t --left-join -f 1 "${CARD_DIR}/AD_CARD_final.xls" "${CARD_DIR}/CARD_mapping-20220905.txt" \
    > "${CARD_DIR}/AD_CARD_mapping.xls"

# 4. Kraken2 16S Analysis
echo "Running Kraken2 16S analysis..."
for fq in "${CLEANREADS}"/*_1.fastq; do
    [[ -f "$fq" ]] || continue
    base=$(basename "$fq" _1.fastq)
    echo "Running Kraken2 for $base"
    kraken2 --db "/media/kunyu/data2/db/metawrap/kraken2/16S_Greengenes_k2db" \
        "$fq" \
        --output "${KRAKEN_DIR}/${base}.txt" \
        --report "${KRAKEN_DIR}/${base}.report.txt" \
        --classified-out "${KRAKEN_DIR}/${base}.16s.fasta" \
        --threads 20
done

# Count 16S reads
echo "Counting 16S reads..."
: > "${KRAKEN_DIR}/16s_reads_number.txt"
for fasta in "${KRAKEN_DIR}"/*.16s.fasta; do
    [[ -f "$fasta" ]] || continue
    count=$(grep -c "^>" "$fasta")
    echo "$(basename "$fasta" .16s.fasta): $count reads" >> "${KRAKEN_DIR}/16s_reads_number.txt"
done

# 5. MetaPhlAn4 and Graphlan
echo "Running MetaPhlAn4 analysis..."
activate_conda metaphlan4
for F in "${CLEANREADS}"/*_1.fastq; do
    [[ -f "$F" ]] || continue
    R="${F%_*}_2.fastq"
    check_file "$R"
    SAMPLE=$(basename "${F%_*}")

    echo "Running MetaPhlAn4 for $SAMPLE"
    metaphlan "$F,$R" \
        --input_type fastq \
        --nproc 50 \
        --bowtie2out "${METAPHLAN_DIR}/${SAMPLE}.bowtie2out" \
        --output_file "${METAPHLAN_DIR}/${SAMPLE}.profile.txt" \
        --bowtie2db "/media/kunyu/data2/db/metaphlan4" \
        --tmp_dir "$TMP_DIR"
done

# Merge MetaPhlAn4 tables
echo "Merging MetaPhlAn4 abundance tables..."
merge_metaphlan_tables.py "${METAPHLAN_DIR}"/*.profile.txt > "${METAPHLAN_DIR}/merged_abundance_table.tsv"

# Calculate diversity metrics
echo "Calculating diversity metrics..."
for metric in "beta bray-curtis" "alpha richness" "alpha shannon" "alpha simpson" "alpha gini"; do
    set -- $metric
    Rscript /home/kunyu/miniconda3/envs/metaphlan4/lib/python3.7/site-packages/metaphlan/utils/calculate_diversity.R \
        -f "${METAPHLAN_DIR}/merged_abundance_table.tsv" -d "$1" -m "$2"
done

# Run Graphlan
echo "Running Graphlan visualization..."
activate_conda graphlan
export2graphlan.py --skip_rows 1,2 \
    -i "${METAPHLAN_DIR}/merged_abundance_table.tsv" \
    --tree "${METAPHLAN_DIR}/merged_abundance.tree.txt" \
    --annotation "${METAPHLAN_DIR}/merged_abundance.annot.txt" \
    --most_abundant 100 --abundance_threshold 1 \
    --least_biomarkers 10 --annotations 5,6 \
    --external_annotations 7 --min_clade_size 1

graphlan_annotate.py --annot "${METAPHLAN_DIR}/merged_abundance.annot.txt" \
    "${METAPHLAN_DIR}/merged_abundance.tree.txt" \
    "${METAPHLAN_DIR}/merged_abundance.xml"

graphlan.py --dpi 300 \
    "${METAPHLAN_DIR}/merged_abundance.xml" \
    "${METAPHLAN_DIR}/merged_abundance.pdf" --external_legends

echo "Reads module completed successfully!"