#!/bin/bash
#SBATCH --job-name=AD_contigs_analysis
#SBATCH --output=AD_contigs_analysis-%j.log
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=50
#SBATCH --mem=200G

# Set strict error handling
set -euo pipefail

# ---------------------- Path Settings ----------------------
WORKDIR="/media/kunyu/data2/03AD"
READS_DIR="${WORKDIR}/01Reads_fq"
ASSEMBLY_DIR="${WORKDIR}/05ASSEMBLY"
BINS_DIR="${WORKDIR}/06MAGs"
BINS_STAGE1="${BINS_DIR}/01bins"
BINS_REFINED="${BINS_DIR}/02refine_bins"
ARCS_DIR="${WORKDIR}/07ARCs_analysis"
ORFS_GENES="${ARCS_DIR}/ORFs_genes"
ORFS_PROT="${ARCS_DIR}/ORFs_protein"
ORFS_NUCL="${ARCS_DIR}/ORFs_nucl"
SARG_DIR="${ORFS_NUCL}/01SARG"
CARD_DIR="${ORFS_NUCL}/02card"
ICE_DIR="${ORFS_NUCL}/04ICEs"
MGES_DIR="${ORFS_NUCL}/05MGEs"
CAT_DIR="${ARCS_DIR}/07CAT"
PLASFLOW_DIR="${ARCS_DIR}/08plasflow"
METAWRAP_OUT="${WORKDIR}/11metawrap-classify"
TMP_DIR="${WORKDIR}/tmp"

# ---------------------- Utility Functions ----------------------
activate_conda() {
    local env=$1
    conda deactivate >/dev/null 2>&1 || true
    conda activate "$env" || { echo "Error: Failed to activate Conda env: $env"; exit 1; }
}

check_file() {
    local file=$1
    if [[ ! -f "$file" ]]; then
        echo "Error: File $file not found"
        exit 1
    fi
}

ensure_dir() {
    for dir in "$@"; do
        mkdir -p "$dir" || { echo "Error: Failed to create directory: $dir"; exit 1; }
    done
}

# ---------------------- Initialize ----------------------
ensure_dir "$READS_DIR" "$ASSEMBLY_DIR" "$BINS_STAGE1" "$BINS_REFINED" "$ARCS_DIR" \
           "$ORFS_GENES" "$ORFS_PROT" "$ORFS_NUCL" "$SARG_DIR" "$CARD_DIR" \
           "$ICE_DIR" "$MGES_DIR" "$CAT_DIR" "$PLASFLOW_DIR" "$METAWRAP_OUT" "$TMP_DIR"
activate_conda tools

# 1. Assembly and Binning
echo "Running assembly and binning..."
for F in "${READS_DIR}"/*_1.fastq; do
    [[ -f "$F" ]] || continue
    R="${F%_*}_2.fastq"
    check_file "$F"
    check_file "$R"
    SAMPLE=$(basename "${F%_*}")

    echo "Assembling sample: $SAMPLE"
    metawrap assembly -1 "$F" -2 "$R" -m 600 -t 30 -o "${ASSEMBLY_DIR}/${SAMPLE}"

    echo "Binning sample: $SAMPLE"
    metawrap binning -o "${BINS_STAGE1}/${SAMPLE}" -t 30 -m 900 \
        -a "${ASSEMBLY_DIR}/${SAMPLE}/final_assembly.fasta" \
        --metabat2 --concoct --maxbin2 "${READS_DIR}/${SAMPLE}"*.fastq --universal

    echo "Refining bins for sample: $SAMPLE"
    metawrap bin_refinement -o "${BINS_REFINED}/${SAMPLE}" -t 30 \
        -A "${BINS_STAGE1}/${SAMPLE}/maxbin2_bins/" \
        -B "${BINS_STAGE1}/${SAMPLE}/metabat2_bins/" \
        -C "${BINS_STAGE1}/${SAMPLE}/concoct_bins/" -c 50 -x 20
done

# 2. Prodigal Gene Prediction
echo "Running Prodigal gene prediction..."
for F in "${READS_DIR}"/*_1.fastq; do
    [[ -f "$F" ]] || continue
    SAMPLE=$(basename "${F%_*}")
    ASSEMBLY="${ASSEMBLY_DIR}/${SAMPLE}/final_assembly.fasta"

    if [[ ! -f "$ASSEMBLY" ]]; then
        echo "Warning: Assembly not found for $SAMPLE at $ASSEMBLY. Skipping."
        continue
    fi

    echo "Predicting genes for $SAMPLE"
    prodigal -i "$ASSEMBLY" \
        -o "${ORFS_GENES}/${SAMPLE}.genes" \
        -a "${ORFS_PROT}/${SAMPLE}.pro.fa" \
        -d "${ORFS_NUCL}/${SAMPLE}.fa" -p meta
done

# 3. SARG Analysis
echo "Running SARG analysis on ORFs..."
check_file "${SARG_DIR}/SARG.2.2.fasta"
diamond makedb --in "${SARG_DIR}/SARG.2.2.fasta" -d "${SARG_DIR}/SARG.2.2_nr" --threads 10

for F in "${ORFS_NUCL}"/*.fa; do
    [[ -f "$F" ]] || continue
    SAMPLE=$(basename "${F%.*}")
    echo "Running Diamond BLASTX for SARG on $SAMPLE"
    diamond blastx -d "${SARG_DIR}/SARG.2.2_nr" -q "$F" -o "${SARG_DIR}/${SAMPLE}.txt" \
        --evalue 1e-10 --query-cover 70 --id 80 -k 1 --threads 10

    awk -F'\t' '{print $1 "\t" $2}' "${SARG_DIR}/${SAMPLE}.txt" | \
        awk 'BEGIN{FS=OFS="\t"}{gsub("_","\t",$1)}1' | \
        cut -f 1-6,8 | \
        awk -F " " '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6" "$7}' | \
        sort -k 1n > "${SARG_DIR}/${SAMPLE}_ARGs.txt"

    sort "${SARG_DIR}/${SAMPLE}_ARGs.txt" > "${SARG_DIR}/${SAMPLE}_ARGs_contigs.txt"
    cut -d' ' -f1 "${SARG_DIR}/${SAMPLE}_ARGs_contigs.txt" | sort -u > "${SARG_DIR}/${SAMPLE}_ARGs_list.txt"

    if [[ -s "${SARG_DIR}/${SAMPLE}_ARGs_list.txt" ]]; then
        seqtk subseq "${ASSEMBLY_DIR}/${SAMPLE}/final_assembly.fasta" "${SARG_DIR}/${SAMPLE}_ARGs_list.txt" | \
            sed "s/k141/${SAMPLE}/g" > "${SARG_DIR}/${SAMPLE}_ARGs.fa"
    else
        echo "No ARG contigs found for $SAMPLE (SARG)"
    fi
done

if compgen -G "${SARG_DIR}/*_ARGs.fa" >/dev/null; then
    cat "${SARG_DIR}"/*_ARGs.fa > "${SARG_DIR}/AD_SARG_ARCs.fa"
fi

# 4. CARD Analysis
echo "Running CARD analysis on ORFs..."
check_file "${CARD_DIR}/CARD3.1.4.fasta"
diamond makedb --in "${CARD_DIR}/CARD3.1.4.fasta" -d "${CARD_DIR}/card3.1.4_nr" --threads 10

for F in "${ORFS_NUCL}"/*.fa; do
    [[ -f "$F" ]] || continue
    SAMPLE=$(basename "${F%.*}")
    echo "Running Diamond BLASTX for CARD on $SAMPLE"
    diamond blastx -d "${CARD_DIR}/card3.1.4_nr" -q "$F" -o "${CARD_DIR}/${SAMPLE}.txt" \
        --evalue 1e-5 --query-cover 70 --id 80 -k 1 --threads 10

    awk -F'\t' '{print $1 "\t" $2}' "${CARD_DIR}/${SAMPLE}.txt" | \
        awk 'BEGIN{FS=OFS="\t"}{gsub("_","\t",$1)}1' | \
        cut -f 1-6,8 | \
        awk -F " " '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6" "$7}' | \
        sort -k 1n > "${CARD_DIR}/${SAMPLE}_ARGs.txt"

    sort "${CARD_DIR}/${SAMPLE}_ARGs.txt" > "${CARD_DIR}/${SAMPLE}_ARGs_contigs.txt"
    cut -d' ' -f1 "${CARD_DIR}/${SAMPLE}_ARGs_contigs.txt" | sort -u > "${CARD_DIR}/${SAMPLE}_ARGs_list.txt"

    if [[ -s "${CARD_DIR}/${SAMPLE}_ARGs_list.txt" ]]; then
        seqtk subseq "${ASSEMBLY_DIR}/${SAMPLE}/final_assembly.fasta" "${CARD_DIR}/${SAMPLE}_ARGs_list.txt" | \
            sed "s/k141/${SAMPLE}/g" > "${CARD_DIR}/${SAMPLE}_ARGs.fa"
    else
        echo "No ARG contigs found for $SAMPLE (CARD)"
    fi
done

if compgen -G "${CARD_DIR}/*_ARGs.fa" >/dev/null; then
    cat "${CARD_DIR}"/*_ARGs.fa > "${CARD_DIR}/AD_card_ARCs.fa"
fi

# Combine ARCs
echo "Combining SARG and CARD ARCs..."
if [[ -f "${SARG_DIR}/AD_SARG_ARCs.fa" && -f "${CARD_DIR}/AD_card_ARCs.fa" ]]; then
    cat "${SARG_DIR}/AD_SARG_ARCs.fa" "${CARD_DIR}/AD_card_ARCs.fa" > "${ARCS_DIR}/AD_ARCs.fasta"
elif [[ -f "${SARG_DIR}/AD_SARG_ARCs.fa" ]]; then
    cp "${SARG_DIR}/AD_SARG_ARCs.fa" "${ARCS_DIR}/AD_ARCs.fasta"
elif [[ -f "${CARD_DIR}/AD_card_ARCs.fa" ]]; then
    cp "${CARD_DIR}/AD_card_ARCs.fa" "${ARCS_DIR}/AD_ARCs.fasta"
else
    echo "Warning: No ARCs generated from SARG or CARD"
fi

# 5. CD-HIT and Prodigal on ARCs
if [[ -f "${ARCS_DIR}/AD_ARCs.fasta" ]]; then
    echo "Running CD-HIT on combined ARCs..."
    cd-hit-est -i "${ARCS_DIR}/AD_ARCs.fasta" -o "${ARCS_DIR}/AD_ARCs_nr.fasta" \
        -c 1.0 -n 10 -M 0 -T 8

    echo "Predicting proteins on non-redundant ARCs..."
    prodigal -i "${ARCS_DIR}/AD_ARCs_nr.fasta" \
        -o "${ARCS_DIR}/AD_ARCs_nr.genes" \
        -a "${ARCS_DIR}/AD_ARCs_nr.pro.fa" \
        -d "${ARCS_DIR}/AD_ARCs_nr.nul.fa" -p meta
fi

# 6. ARCs Annotation (SARG, CARD, ICEs, MGEs)
declare -A databases=(
    ["01SARG/SARG.2.2.fasta"]="SARG.2.2_nr"
    ["02card/CARD3.1.4.fasta"]="card3.1.4_nr"
    ["04ICEs/ICE3.0.fasta"]="ICE3.0_nr"
    ["05MGEs/MGEs.fasta"]="MGEs"
)

for db_path in "${!databases[@]}"; do
    [[ -f "${WORKDIR}/${db_path}" ]] || { echo "Error: Database ${WORKDIR}/${db_path} not found"; continue; }
    db_name=${databases[$db_path]}
    out_dir="${WORKDIR}/${db_path%/*}"
    if [[ "$db_name" == "MGEs" ]]; then
        makeblastdb -dbtype nucl -in "${WORKDIR}/${db_path}" -input_type fasta \
            -title MGEs -out "${out_dir}/${db_name}"
    else
        diamond makedb --in "${WORKDIR}/${db_path}" -d "${out_dir}/${db_name}" --threads 10
    fi
done

if [[ -f "${ARCS_DIR}/AD_ARCs_nr.nul.fa" ]]; then
    echo "Annotating ARCs against SARG, CARD, ICEs, and MGEs..."
    for db_type in "${!databases[@]}"; do
        db_name=${databases[$db_type]}
        out_dir="${WORKDIR}/${db_type%/*}"
        out_file="${out_dir}/AD_ARCs-${db_name%.nr}.txt"
        if [[ "$db_name" == "MGEs" ]]; then
            blastn -query "${ARCS_DIR}/AD_ARCs_nr.nul.fa" -db "${out_dir}/${db_name}" \
                -out "$out_file" -evalue 1e-10 -perc_identity 80 -num_threads 10 -outfmt 6 -max_target_seqs 1
        else
            diamond blastx -d "${out_dir}/${db_name}" -q "${ARCS_DIR}/AD_ARCs_nr.nul.fa" \
                -o "$out_file" --evalue 1e-10 --query-cover 70 --id 80 -k 1 --threads 10
        fi

        awk -F'\t' '{print $1 "\t" $2}' "$out_file" | \
            awk 'BEGIN{FS=OFS="\t"}{gsub("_","\t",$1)}1' | \
            cut -f 1-6,8 | \
            awk -F " " '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6" "$7}' | \
            sort -k 1n > "${out_file%.txt}_list.txt"

        if command -v datamash >/dev/null 2>&1; then
            datamash -sW -g1 collapse 2 < "${out_file%.txt}_list.txt" > "${out_file%.txt}_list2.txt"
        fi
    done
fi

# 7. CAT Classification
echo "Running CAT classification..."
CAT_DB="/media/kunyu/data2/db/CAT_prepare_20210107/2021-01-07_CAT_database"
CAT_TAX="/media/kunyu/data2/db/CAT_prepare_20210107/2021-01-07_taxonomy"

if [[ -f "${ARCS_DIR}/AD_ARCs_nr.fasta" ]]; then
    CAT contigs -c "${ARCS_DIR}/AD_ARCs_nr.fasta" -d "$CAT_DB" -t "$CAT_TAX" \
        -p "${ARCS_DIR}/AD_ARCs_nr.pro.fa" --cpus 10
    CAT add_names -i out.CAT.contig2classification.txt -o "${CAT_DIR}/CAT_contigs_classification.txt" -t "$CAT_TAX"
    CAT add_names -i out.CAT.ORF2LCA.txt -o "${CAT_DIR}/CAT_ORFs_classification.txt" -t "$CAT_TAX"

    mv out.CAT.* "${CAT_DIR}/" || true
    rm -rf out.CAT.alignment.diamond || true
fi

# 8. PlasFlow Analysis
echo "Running PlasFlow analysis..."
activate_conda plasflow
if [[ -f "${ARCS_DIR}/AD_ARCs_nr.fasta" ]]; then
    PlasFlow.py --input "${ARCS_DIR}/AD_ARCs_nr.fasta" --output "${PLASFLOW_DIR}/AD_ARCs_nr.plasflow.txt" --threshold 0.7
fi

# 9. MetaWrap Classification
echo "Running MetaWrap classification..."
activate_conda metawrap-env
if [[ -f "${ARCS_DIR}/AD_ARCs_nr.fasta" ]]; then
    cp "${ARCS_DIR}/AD_ARCs_nr.fasta" "${METAWRAP_OUT}/" || true
    metawrap classify_bins -b "${METAWRAP_OUT}" -o "${METAWRAP_OUT}" -t 10
fi

# 10. CoverM Analysis
echo "Running CoverM analysis..."
activate_conda tools
for F in "${READS_DIR}"/*_1.fastq; do
    [[ -f "$F" ]] || continue
    R="${F%_*}_2.fastq"
    check_file "$R"
    SAMPLE=$(basename "${F%_*}")

    echo "Mapping reads for sample: $SAMPLE"
    TMPDIR="$TMP_DIR" coverm make -r "${ARCS_DIR}/AD_ARCs_nr.fasta" -1 "$F" -2 "$R" -o "${COVERM_DIR}" -t 30

    BAM_IN=$(ls "${COVERM_DIR}"/*.bam 2>/dev/null | grep "${SAMPLE}_1.fastq" || true)
    if [[ -f "$BAM_IN" ]]; then
        coverm filter -b "$BAM_IN" -o "${COVERM_DIR}/${SAMPLE}.bam" \
            --min-read-percent-identity 95 --min-read-aligned-percent 75 --threads 30
        coverm contig --bam-files "${COVERM_DIR}/${SAMPLE}.bam" \
            -o "${COVERM_DIR}/${SAMPLE}.rpkm.txt" --trim-min 10 --trim-max 90 \
            --min-read-percent-identity 95 --min-read-aligned-percent 75 -m rpkm -t 30
        coverm contig --bam-files "${COVERM_DIR}/${SAMPLE}.bam" \
            -o "${COVERM_DIR}/${SAMPLE}.count.txt" --min-covered-fraction 0 -m count -t 30
        coverm contig --bam-files "${COVERM_DIR}/${SAMPLE}.bam" \
            -o "${COVERM_DIR}/${SAMPLE}.bases.txt" --min-covered-fraction 0 -m covered_bases -t 30
    else
        echo "Warning: BAM not found for sample $SAMPLE"
    fi
done

if compgen -G "${COVERM_DIR}/*.rpkm.txt" >/dev/null; then
    paste "${COVERM_DIR}"/*.rpkm.txt > "${COVERM_DIR}/AD_SSU.rpkm.txt"
    paste "${COVERM_DIR}"/*.count.txt > "${COVERM_DIR}/AD_SSU.count.txt"
    paste "${COVERM_DIR}"/*.bases.txt > "${COVERM_DIR}/AD_SSU.bases.txt"
fi

echo "Contigs module completed successfully!"