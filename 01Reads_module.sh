#!/bin/bash
#SBATCH --job-name=AD_pipeline
#SBATCH --output=AD_pipeline-%j.log
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=50
#SBATCH --mem=200G

set -euo pipefail

# ---------------------- 0. Project Path Settings ---------------------- #
WORKDIR="/media/kunyu/data2/03AD"

RAWREADS="$WORKDIR/01RAW_READS"
QC_DIR="$WORKDIR/01READ_QC"
CLEANREADS="$WORKDIR/01Reads_fq"
OUTDIR="$WORKDIR/04reads_based_analysis"

SARG_DIR="$OUTDIR/01SARG"
CARD_DIR="$OUTDIR/02card"
KRAKEN_DIR="$OUTDIR/07kraken2"
METAPHLAN_DIR="$WORKDIR/10metaphlan4"
TMP_DIR="$WORKDIR/tmp"

mkdir -p "$QC_DIR" "$CLEANREADS" "$SARG_DIR" "$CARD_DIR" "$KRAKEN_DIR" "$METAPHLAN_DIR" "$TMP_DIR"

# ---------------------- 1. Read QC ---------------------- #
conda deactivate || true
conda activate base

for F in "$RAWREADS"/*_1.fastq; do 
    R=${F%_*}_2.fastq
    SAMPLE=$(basename "${F%_*}")

    echo "Running metawrap read_qc on sample: $SAMPLE"
    metawrap read_qc -1 "$F" -2 "$R" -t 20 --skip-bmtagger -o "$QC_DIR/$SAMPLE"

    pigz "$F" "$R"
done

#  clean reads
for i in "$QC_DIR"/*; do 
    b=$(basename "$i")
    mv "$i/final_pure_reads_1.fastq" "$CLEANREADS/${b}_1.fastq"
    mv "$i/final_pure_reads_2.fastq" "$CLEANREADS/${b}_2.fastq"
done

#  reads number
conda deactivate
for file in "$CLEANREADS"/*.fastq; do
    echo "$(basename "$file" .fastq): $(grep -c "^@" "$file") reads" \
        >> "$CLEANREADS/AD_Meta_reads_number.txt"
done

# ---------------------- 2. SARG Analysis ---------------------- #
conda activate tools

diamond makedb \
    --in "$SARG_DIR/SARG.2.2.fasta" \
    -d "$SARG_DIR/SARG.2.2_nr"

for fq in "$CLEANREADS"/*.fastq; do
    out="$SARG_DIR/$(basename "$fq")-SARG.txt"
    diamond blastx -d "$SARG_DIR/SARG.2.2_nr" \
        -q "$fq" -o "$out" --evalue 1e-5 --query-cover 75 --id 90 -k 1
done

echo -e "SARG_family\t$(ls "$SARG_DIR"/*-SARG.txt | sed 's/.*\///;s/-SARG.txt//;s/^/count /' | tr '\n' '\t')" \
    > "$SARG_DIR/AD_SARG.xls"

for i in $(cut -f 1 "$SARG_DIR/SARG.2.2.txt" | sort -u); do
    echo -ne "$i\t" >> "$SARG_DIR/AD_SARG.xls"
    for file in "$SARG_DIR"/*-SARG.txt; do
        count=$(grep -c -F "$i" "$file")
        echo -ne "$count\t" >> "$SARG_DIR/AD_SARG.xls"
    done
    echo >> "$SARG_DIR/AD_SARG.xls"
done

awk -F'\t' 'BEGIN {OFS="\t"} NR==1 {print} NR>1 {sum=0; for (i=2; i<=NF; i++) {sum+=$i} if (sum!=0) {print}}' \
    "$SARG_DIR/AD_SARG.xls" > "$SARG_DIR/AD_SARG_final.xls"

csvtk join -t --left-join -f 1 "$SARG_DIR/AD_SARG_final.xls" "$SARG_DIR/SARG_mapping-20220906.txt" \
    > "$SARG_DIR/AD_SARG_mapping.xls"

# ---------------------- 3. CARD Analysis ---------------------- #
diamond makedb \
    --in "$CARD_DIR/CARD3.1.4.fasta" \
    -d "$CARD_DIR/card3.1.4_nr"

for fq in "$CLEANREADS"/*.fastq; do
    out="$CARD_DIR/$(basename "$fq")-CARD.txt"
    diamond blastx -d "$CARD_DIR/card3.1.4_nr" \
        -q "$fq" -o "$out" --evalue 1e-5 --query-cover 75 --id 90 -k 1
done

echo -e "CARD_family\t$(ls "$CARD_DIR"/*-CARD.txt | sed 's/.*\///;s/-CARD.txt//;s/^/count /' | tr '\n' '\t')" \
    > "$CARD_DIR/AD_CARD.xls"

for i in $(cut -f 3 "$CARD_DIR/CARD3.1.4.txt" | sort -u); do
    echo -ne "$i\t" >> "$CARD_DIR/AD_CARD.xls"
    for file in "$CARD_DIR"/*-CARD.txt; do
        count=$(grep -c -F "$i" "$file")
        echo -ne "$count\t" >> "$CARD_DIR/AD_CARD.xls"
    done
    echo >> "$CARD_DIR/AD_CARD.xls"
done

awk -F'\t' 'BEGIN {OFS="\t"} NR==1 {print} NR>1 {sum=0; for (i=2; i<=NF; i++) {sum+=$i} if (sum!=0) {print}}' \
    "$CARD_DIR/AD_CARD.xls" > "$CARD_DIR/AD_CARD_final.xls"

csvtk join -t --left-join -f 1 "$CARD_DIR/AD_CARD_final.xls" "$CARD_DIR/CARD_mapping-20220905.txt" \
    > "$CARD_DIR/AD_CARD_mapping.xls"

# ---------------------- 4. Kraken2 16S ---------------------- #
for fq in "$CLEANREADS"/*.fastq; do
    base=$(basename "$fq")
    kraken2 \
        --db /media/kunyu/data2/db/metawrap/kraken2/16S_Greengenes_k2db \
        "$fq" \
        --output "$KRAKEN_DIR/$base.txt" \
        --report "$KRAKEN_DIR/$base.report.txt" \
        --classified-out "$KRAKEN_DIR/$base.16s.fasta" \
        --threads 20
done

for fasta in "$KRAKEN_DIR"/*.fasta; do
    echo "$(basename "$fasta" .fasta): $(grep -c "^@" "$fasta") reads" >> "$KRAKEN_DIR/16s_reads_number.txt"
done

# ---------------------- 5. MetaPhlAn4 + Graphlan ---------------------- #
conda deactivate
conda activate metaphlan4

for F in "$CLEANREADS"/*_1.fastq; do
    R=${F%_*}_2.fastq
    SAMPLE=$(basename "${F%_*}")
    
    metaphlan "$F","$R" \
        --input_type fastq \
        --nproc 50 \
        --bowtie2out "$METAPHLAN_DIR/$SAMPLE.bowtie2out" \
        --output_file "$METAPHLAN_DIR/$SAMPLE.profile.txt" \
        --bowtie2db /media/kunyu/data2/db/metaphlan4 \
        --tmp_dir "$TMP_DIR"
done

merge_metaphlan_tables.py "$METAPHLAN_DIR"/*.txt > "$METAPHLAN_DIR/merged_abundance_table.tsv"

for metric in "beta bray-curtis" "alpha richness" "alpha shannon" "alpha simpson" "alpha gini"; do
    set -- $metric
    Rscript /home/kunyu/miniconda3/envs/metaphlan4/lib/python3.7/site-packages/metaphlan/utils/calculate_diversity.R \
        -f "$METAPHLAN_DIR/merged_abundance_table.tsv" -d "$1" -m "$2"
done

conda deactivate
conda activate graphlan

export2graphlan.py --skip_rows 1,2 \
    -i "$METAPHLAN_DIR/merged_abundance_table.tsv" \
    --tree "$METAPHLAN_DIR/merged_abundance.tree.txt" \
    --annotation "$METAPHLAN_DIR/merged_abundance.annot.txt" \
    --most_abundant 100 --abundance_threshold 1 \
    --least_biomarkers 10 --annotations 5,6 \
    --external_annotations 7 --min_clade_size 1

graphlan_annotate.py --annot "$METAPHLAN_DIR/merged_abundance.annot.txt" \
    "$METAPHLAN_DIR/merged_abundance.tree.txt" \
    "$METAPHLAN_DIR/merged_abundance.xml"

graphlan.py --dpi 300 \
    "$METAPHLAN_DIR/merged_abundance.xml" \
    "$METAPHLAN_DIR/merged_abundance.pdf" --external_legends
