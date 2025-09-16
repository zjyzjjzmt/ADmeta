#!/bin/bash
#SBATCH --job-name=myjob_name
#SBATCH --chdir=/work/<USERNAME>
#SBATCH --output=/work/%u/%x-%j.log
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=50
#SBATCH --mem=200G

set -euo pipefail
IFS=$'\n\t'

# ------------------- Project path settings (modify this single line) -------------------
WORKDIR="/media/kunyu/data2/03AD"

# ------------------- Derived directories -------------------
RAWREADS="$WORKDIR/01RAW_READS"        # raw reads
READS_DIR="$WORKDIR/01Reads_fq"        # clean reads
QC_DIR="$WORKDIR/01READ_QC"
ASSEMBLY_DIR="$WORKDIR/05ASSEMBLY"
BINS_DIR="$WORKDIR/06MAGs"
BINS_STAGE1="$BINS_DIR/01bins"
BINS_REFINED="$BINS_DIR/02refine_bins"
OUTDIR="$WORKDIR/04reads_based_analysis"
ORFS_GENES="$ASSEMBLY_DIR/ORFs_genes"
ORFS_PROT="$ASSEMBLY_DIR/ORFs_protein"
ORFS_NUCL="$ASSEMBLY_DIR/ORFs_nucl"
SARG_DIR="$ORFS_NUCL/01SARG"
CARD_DIR="$ORFS_NUCL/02card"
ICE_DIR="$ORFS_NUCL/04ICEs"
MGES_DIR="$ORFS_NUCL/05MGEs"
ARCS_DIR="$WORKDIR/07ARCs_analysis"
CAT_DIR="$ARCS_DIR/07CAT"
PLASFLOW_DIR="$ARCS_DIR/08plasflow"
METAWRAP_CLASSIFY_DIR="$WORKDIR/10contigs"  # copy AD_ARCs_nr.fasta here
METAWRAP_OUT="$WORKDIR/11metawrap-classify"
COVERM_DIR="$ARCS_DIR/12coverm"
METAPHLAN_DIR="$WORKDIR/10metaphlan4"
TMP_DIR="$WORKDIR/tmp"

# ------------------- Make sure directories exist -------------------
mkdir -p "$READS_DIR" "$QC_DIR" "$ASSEMBLY_DIR" "$BINS_STAGE1" "$BINS_REFINED" \
         "$OUTDIR" "$ORFS_GENES" "$ORFS_PROT" "$ORFS_NUCL" "$SARG_DIR" "$CARD_DIR" \
         "$ICE_DIR" "$MGES_DIR" "$ARCS_DIR" "$CAT_DIR" "$PLASFLOW_DIR" "$METAWRAP_CLASSIFY_DIR" \
         "$METAWRAP_OUT" "$COVERM_DIR" "$METAPHLAN_DIR" "$TMP_DIR"

# ------------------- Conda / environment utilities -------------------
source deactivate
conda activate base
conda activate tools

# ------------------- Assembly & binning (integrated before Prodigal) -------------------
for F in "$READS_DIR"/*_1.fastq; do
    R=${F%_*}_2.fastq
    BASE=$(basename "$F")
    SAMPLE=${BASE%_*}

    echo "=== Assembly and binning for sample: $SAMPLE ==="

    # metawrap assembly
    metawrap assembly -1 "$F" -2 "$R" -m 600 -t 30 -o "$ASSEMBLY_DIR/$SAMPLE"

    # metawrap binning
    metawrap binning -o "$BINS_STAGE1/$SAMPLE" -t 30 -m 900 \
        -a "$ASSEMBLY_DIR/$SAMPLE/final_assembly.fasta" \
        --metabat2 --concoct --maxbin2 "$READS_DIR/$SAMPLE"*.fastq --universal

    # bin refinement
    metawrap bin_refinement -o "$BINS_REFINED/$SAMPLE" -t 30 \
        -A "$BINS_STAGE1/$SAMPLE/maxbin2_bins/" \
        -B "$BINS_STAGE1/$SAMPLE/metabat2_bins/" \
        -C "$BINS_STAGE1/$SAMPLE/concoct_bins/" -c 50 -x 20
done

# ------------------- Contigs Prodigal -------------------
for F in "$READS_DIR"/*_1.fastq; do
    BASE=$(basename "$F")
    SAMPLE=${BASE%_*}

    echo "=== Prodigal on assembly of sample: $SAMPLE ==="
    if [[ ! -f "$ASSEMBLY_DIR/$SAMPLE/final_assembly.fasta" ]]; then
        echo "Warning: assembly not found for $SAMPLE at $ASSEMBLY_DIR/$SAMPLE/final_assembly.fasta. Skipping."
        continue
    fi

    prodigal -i "$ASSEMBLY_DIR/$SAMPLE/final_assembly.fasta" \
        -o "$ASSEMBLY_DIR/$SAMPLE/$SAMPLE.genes" \
        -a "$ASSEMBLY_DIR/$SAMPLE/$SAMPLE.pro.fa" \
        -d "$ASSEMBLY_DIR/$SAMPLE/$SAMPLE.fa" -p meta

    mv "$ASSEMBLY_DIR/$SAMPLE/$SAMPLE.genes" "$ORFS_GENES/" || true
    mv "$ASSEMBLY_DIR/$SAMPLE/$SAMPLE.pro.fa" "$ORFS_PROT/" || true
    mv "$ASSEMBLY_DIR/$SAMPLE/$SAMPLE.fa" "$ORFS_NUCL/" || true
done

# ------------------- AD-SARG (on ORFs_nucl) -------------------
echo "=== AD-SARG analysis on predicted ORFs (nucleotide) ==="
diamond makedb --in "$SARG_DIR/SARG.2.2.fasta" -d "$SARG_DIR/SARG.2.2_nr"

for F in "$ORFS_NUCL"/*.fa; do
    [[ -f "$F" ]] || continue
    BASE=$(basename "$F")
    SAMPLE=${BASE%.*}
    echo "Running diamond blastx for SARG on $BASE"
    diamond blastx -d "$SARG_DIR/SARG.2.2_nr" -q "$F" -o "$SARG_DIR/$SAMPLE.fa.txt" \
        --evalue 1e-10 --query-cover 70 --id 80 -k 1

    cut -f 1,2 "$SARG_DIR/$SAMPLE.fa.txt" \
      | awk 'BEGIN{FS=OFS="\t"}{gsub("_","\t",$1)}1' \
      | cut -f 1,2,3,4,5,6,8 \
      | awk -F " " '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6" "$7}' \
      | sort -k 1n > "$SARG_DIR/${SAMPLE}_ARGs.txt"

    sort "$SARG_DIR/${SAMPLE}_ARGs.txt" > "$SARG_DIR/${SAMPLE}_ARGs_contigs.txt"
    cut -d' ' -f1 "$SARG_DIR/${SAMPLE}_ARGs_contigs.txt" | sort -u > "$SARG_DIR/${SAMPLE}_ARGs_list.txt"

    if [[ -s "$SARG_DIR/${SAMPLE}_ARGs_list.txt" ]]; then
        seqtk subseq "$ASSEMBLY_DIR/$SAMPLE/final_assembly.fasta" "$SARG_DIR/${SAMPLE}_ARGs_list.txt" \
            | sed "s/k141/$SAMPLE/g" > "$SARG_DIR/${SAMPLE}_ARGs.fa"
    else
        echo "No ARG contigs found for $SAMPLE (SARG)."
    fi
done

if compgen -G "$SARG_DIR/*_ARGs.fa" >/dev/null; then
    cat "$SARG_DIR"/*_ARGs.fa > "$SARG_DIR/AD_SARG_ARCs.fa"
else
    echo "No SARG ARCs fasta files to concatenate."
fi

# ------------------- AD-CARD (on ORFs_nucl) -------------------
echo "=== AD-CARD analysis on predicted ORFs (nucleotide) ==="
diamond makedb --in "$CARD_DIR/CARD3.1.4.fasta" -d "$CARD_DIR/card3.1.4_nr"

for F in "$ORFS_NUCL"/*.fa; do
    [[ -f "$F" ]] || continue
    BASE=$(basename "$F")
    SAMPLE=${BASE%.*}
    echo "Running diamond blastx for CARD on $BASE"
    diamond blastx -d "$CARD_DIR/card3.1.4_nr" -q "$F" -o "$CARD_DIR/$SAMPLE.fa.txt" \
        --evalue 1e-5 --query-cover 70 --id 80 -k 1

    cut -f 1,2 "$CARD_DIR/$SAMPLE.fa.txt" \
      | awk 'BEGIN{FS=OFS="\t"}{gsub("_","\t",$1)}1' \
      | cut -f 1,2,3,4,5,6,8 \
      | awk -F " " '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6" "$7}' \
      | sort -k 1n > "$CARD_DIR/${SAMPLE}_ARGs.txt"

    sort "$CARD_DIR/${SAMPLE}_ARGs.txt" > "$CARD_DIR/${SAMPLE}_ARGs_contigs.txt"
    cut -d' ' -f1 "$CARD_DIR/${SAMPLE}_ARGs_contigs.txt" | sort -u > "$CARD_DIR/${SAMPLE}_ARGs_list.txt"

    if [[ -s "$CARD_DIR/${SAMPLE}_ARGs_list.txt" ]]; then
        seqtk subseq "$ASSEMBLY_DIR/$SAMPLE/final_assembly.fasta" "$CARD_DIR/${SAMPLE}_ARGs_list.txt" \
            | sed "s/k141/$SAMPLE/g" > "$CARD_DIR/${SAMPLE}_ARGs.fa"
    else
        echo "No ARG contigs found for $SAMPLE (CARD)."
    fi
done

if compgen -G "$CARD_DIR/*_ARGs.fa" >/dev/null; then
    cat "$CARD_DIR"/*_ARGs.fa > "$CARD_DIR/AD_card_ARCs.fa"
else
    echo "No CARD ARCs fasta files to concatenate."
fi

if [[ -f "$SARG_DIR/AD_SARG_ARCs.fa" ]] && [[ -f "$CARD_DIR/AD_card_ARCs.fa" ]]; then
    cat "$SARG_DIR/AD_SARG_ARCs.fa" "$CARD_DIR/AD_card_ARCs.fa" > "$ARCS_DIR/AD_ARCs.fasta"
elif [[ -f "$SARG_DIR/AD_SARG_ARCs.fa" ]]; then
    cp "$SARG_DIR/AD_SARG_ARCs.fa" "$ARCS_DIR/AD_ARCs.fasta"
elif [[ -f "$CARD_DIR/AD_card_ARCs.fa" ]]; then
    cp "$CARD_DIR/AD_card_ARCs.fa" "$ARCS_DIR/AD_ARCs.fasta"
else
    echo "No ARCs generated from SARG or CARD. Exiting ARCs downstream steps."
fi

# ------------------- CD-HIT & prodigal on combined ARCs -------------------
if [[ -f "$ARCS_DIR/AD_ARCs.fasta" ]]; then
    echo "Running cd-hit-est on combined ARCs"
    cd-hit-est -i "$ARCS_DIR/AD_ARCs.fasta" -o "$ARCS_DIR/AD_ARCs_nr.fasta" -c 1.0 -n 10 -M 0 -T 8

    echo "Predict proteins on non-redundant ARCs"
    prodigal -i "$ARCS_DIR/AD_ARCs_nr.fasta" \
        -o "$ARCS_DIR/AD_ARCs_nr.genes" \
        -a "$ARCS_DIR/AD_ARCs_nr.pro.fa" \
        -d "$ARCS_DIR/AD_ARCs_nr.nul.fa" -p meta
else
    echo "Skipping CD-HIT and ARCs prodigal because $ARCS_DIR/AD_ARCs.fasta not found."
fi

# ------------------- AD-01SARG (annotate ARCs_nr against SARG) -------------------
if [[ -f "$ARCS_DIR/AD_ARCs_nr.nul.fa" ]]; then
    diamond makedb --in "$WORKDIR/01SARG/SARG.2.2.fasta" -d "$WORKDIR/01SARG/SARG.2.2_nr"
    diamond blastx -d "$WORKDIR/01SARG/SARG.2.2_nr" -q "$ARCS_DIR/AD_ARCs_nr.nul.fa" \
        -o "$WORKDIR/01SARG/AD_ARCs-SARG.txt" --evalue 1e-10 --query-cover 70 --id 80 -k 1

    cut -f 1,2 "$WORKDIR/01SARG/AD_ARCs-SARG.txt" \
      | awk 'BEGIN{FS=OFS="\t"}{gsub("_","\t",$1)}1' \
      | cut -f 1,2,3,4,5,6,8 \
      | awk -F " " '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6" "$7}' \
      | sort -k 1n > "$WORKDIR/01SARG/AD_ARCs-SARG_list.txt"

    # collapse using datamash 
    if command -v datamash >/dev/null 2>&1; then
        cat "$WORKDIR/01SARG/AD_ARCs-SARG_list.txt" | datamash -sW -g1 collapse 2 > "$WORKDIR/01SARG/AD_ARCs-SARG_list2.txt"
    fi
else
    echo "Skipping AD-01SARG: ARCs_nr.nul.fa not found."
fi

# ------------------- AD-02CARD (annotate ARCs_nr against CARD) -------------------
if [[ -f "$ARCS_DIR/AD_ARCs_nr.nul.fa" ]]; then
    diamond makedb --in "$WORKDIR/02card/CARD3.1.4.fasta" -d "$WORKDIR/02card/card3.1.4_nr"
    diamond blastx -d "$WORKDIR/02card/card3.1.4_nr" -q "$ARCS_DIR/AD_ARCs_nr.nul.fa" \
        -o "$WORKDIR/02card/AD_ARCs-card.txt" --evalue 1e-10 --query-cover 70 --id 80 -k 1

    cut -f 1,2 "$WORKDIR/02card/AD_ARCs-card.txt" \
      | awk 'BEGIN{FS=OFS="\t"}{gsub("_","\t",$1)}1' \
      | cut -f 1,2,3,4,5,6,8 \
      | awk -F " " '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6" "$7}' \
      | sort -k 1n > "$WORKDIR/02card/AD_ARCs-card_list.txt"

    awk '{split ($2, T, "|"); $2 = T[4]}1' OFS="\t" "$WORKDIR/02card/AD_ARCs-card_list.txt" \
        | ( command -v datamash >/dev/null 2>&1 && datamash -sW -g1 collapse 2 > "$WORKDIR/02card/AD_ARCs-card_list2.txt" || cat > "$WORKDIR/02card/AD_ARCs-card_list2.txt" )
else
    echo "Skipping AD-02CARD: ARCs_nr.nul.fa not found."
fi

# ------------------- AD-04ICEs (annotate ARCs against ICEs) -------------------
if [[ -f "$ARCS_DIR/AD_ARCs_nr.nul.fa" ]]; then
    diamond makedb --in "$WORKDIR/04ICEs/ICE3.0.fasta" -d "$WORKDIR/04ICEs/ICE3.0_nr"
    diamond blastx -d "$WORKDIR/04ICEs/ICE3.0_nr" -q "$ARCS_DIR/AD_ARCs_nr.nul.fa" \
        -o "$WORKDIR/04ICEs/AD_ARCs-ICE.txt" --evalue 1e-10 --query-cover 70 --id 80 -k 1

    cut -f 1,2 "$WORKDIR/04ICEs/AD_ARCs-ICE.txt" \
      | awk 'BEGIN{FS=OFS="\t"}{gsub("_","\t",$1)}1' \
      | cut -f 1,2,3,4,5,6,8 \
      | awk -F " " '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6" "$7}' \
      | sort -k 1n > "$WORKDIR/04ICEs/AD_ARCs-ICE_list.txt"

    if command -v datamash >/dev/null 2>&1; then
        cat "$WORKDIR/04ICEs/AD_ARCs-ICE_list.txt" | datamash -sW -g1 collapse 2 > "$WORKDIR/04ICEs/AD_ARCs-ICE_list2.txt"
    fi
else
    echo "Skipping AD-04ICEs: ARCs_nr.nul.fa not found."
fi

# ------------------- AD-05MGEs (blastn) -------------------
if [[ -f "$ARCS_DIR/AD_ARCs_nr.nul.fa" ]]; then
    if [[ ! -f "$MGES_DIR/MGEs.nhr" ]]; then
        makeblastdb -dbtype nucl -in "$MGES_DIR/MGEs.fasta" -input_type fasta -title MGEs -out "$MGES_DIR/MGEs"
    fi
    blastn -query "$ARCS_DIR/AD_ARCs_nr.nul.fa" -db "$MGES_DIR/MGEs" \
        -out "$MGES_DIR/AD_ARCs-MGEs.txt" -evalue 1e-10 -perc_identity 80 -num_threads 10 -outfmt 6 -max_target_seqs 1

    cut -f 1,2 "$MGES_DIR/AD_ARCs-MGEs.txt" \
      | awk 'BEGIN{FS=OFS="\t"}{gsub("_","\t",$1)}1' \
      | cut -f 1,2,3,4,5,6,8 \
      | awk -F " " '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6" "$7}' \
      | sort -k 1n > "$MGES_DIR/AD_ARCs-MGEs_list.txt"

    if command -v datamash >/dev/null 2>&1; then
        cat "$MGES_DIR/AD_ARCs-MGEs_list.txt" | datamash -sW -g1 collapse 2 > "$MGES_DIR/AD_ARCs-MGEs_list2.txt"
    fi
fi

# ------------------- AD-09CAT (classification with CAT) -------------------
# Ensure CAT database paths are correct
CAT_DB="/media/kunyu/data2/db/CAT_prepare_20210107/2021-01-07_CAT_database"
CAT_TAX="/media/kunyu/data2/db/CAT_prepare_20210107/2021-01-07_taxonomy"

if [[ -f "$ARCS_DIR/AD_ARCs_nr.fasta" ]]; then
    CAT contigs -c "$ARCS_DIR/AD_ARCs_nr.fasta" -d "$CAT_DB" -t "$CAT_TAX" -p "$ARCS_DIR/AD_ARCs_nr.pro.fa"
    CAT add_names -i out.CAT.contig2classification.txt -o "$CAT_DIR/CAT_contigs_classification.txt" -t "$CAT_TAX"
    CAT add_names -i out.CAT.ORF2LCA.txt -o "$CAT_DIR/CAT_ORFs_classification.txt" -t "$CAT_TAX"

    # move outputs to CAT_DIR
    mv out.CAT.contig2classification.txt "$CAT_DIR/" || true
    mv out.CAT.predicted_proteins.faa "$CAT_DIR/" || true
    mv out.CAT.predicted_proteins.gff "$CAT_DIR/" || true
    mv out.CAT.ORF2LCA.txt "$CAT_DIR/" || true
    mv out.CAT.log "$CAT_DIR/" || true
    rm -rf out.CAT.alignment.diamond || true
else
    echo "Skipping CAT: $ARCS_DIR/AD_ARCs_nr.fasta not found."
fi

# ------------------- AD-08plasflow -------------------
conda deactivate || true
conda activate plasflow
if [[ -f "$ARCS_DIR/AD_ARCs_nr.fasta" ]]; then
    PlasFlow.py --input "$ARCS_DIR/AD_ARCs_nr.fasta" --output "$PLASFLOW_DIR/AD_ARCs_nr.plasflow.txt" --threshold 0.7
else
    echo "Skipping PlasFlow: $ARCS_DIR/AD_ARCs_nr.fasta not found."
fi

# ------------------- AD-08metawrap classify_bins -------------------
conda deactivate || true
conda activate metawrap-env
if [[ -f "$ARCS_DIR/AD_ARCs_nr.fasta" ]]; then
    cp "$ARCS_DIR/AD_ARCs_nr.fasta" "$METAWRAP_CLASSIFY_DIR/" || true
    metawrap classify_bins -b "$METAWRAP_CLASSIFY_DIR" -o "$METAWRAP_OUT" -t 10
else
    echo "Skipping metawrap classify_bins: $ARCS_DIR/AD_ARCs_nr.fasta not found."
fi

# ------------------- AD-ARCs-coverm (mapping reads back to ARCs) -------------------
conda deactivate || true
conda activate tools

# run coverm per sample: make, filter, contig metrics
for F in "$READS_DIR"/*_1.fastq; do
    [[ -f "$F" ]] || continue
    R=${F%_*}_2.fastq
    BASE=$(basename "$F")
    SAMPLE=${BASE%_*}
    echo "=== COVERM mapping for sample: $SAMPLE ==="

    # make BAMs and index
    TMPDIR="$TMP_DIR" coverm make -r "$ARCS_DIR/AD_ARCs_nr.fasta" -1 "$F" -2 "$R" -o "$COVERM_DIR" -t 30

    # the original script referenced a specific bam name pattern:
    BAM_IN="$COVERM_DIR/AD_ARCs_nr.fasta.$BASE.bam"
    if [[ ! -f "$BAM_IN" ]]; then
        # try to find any bam produced for this sample
        BAM_IN=$(ls "$COVERM_DIR"/*.bam 2>/dev/null | grep "$BASE" || true)
    fi

    if [[ -f "$BAM_IN" ]]; then
        coverm filter -b "$BAM_IN" -o "$COVERM_DIR/$SAMPLE.bam" --min-read-percent-identity 95 --min-read-aligned-percent 75 --threads 30

        coverm contig --bam-files "$COVERM_DIR/$SAMPLE.bam" -o "$COVERM_DIR/$SAMPLE.rpkm.txt" --trim-min 10 --trim-max 90 --min-read-percent-identity 95 --min-read-aligned-percent 75 -m rpkm -t 30

        coverm contig --bam-files "$COVERM_DIR/$SAMPLE.bam" -o "$COVERM_DIR/$SAMPLE.count.txt" --min-covered-fraction 0 -m count -t 30

        coverm contig --bam-files "$COVERM_DIR/$SAMPLE.bam" -o "$COVERM_DIR/$SAMPLE.bases.txt" --min-covered-fraction 0 -m covered_bases -t 30
    else
        echo "BAM not found for sample $SAMPLE (expected $BAM_IN). Skipping coverm contig steps."
    fi
done

if compgen -G "$COVERM_DIR/*.rpkm.txt" >/dev/null; then
    paste "$COVERM_DIR"/*.rpkm.txt > "$COVERM_DIR/AD_SSU.rpkm.txt"
fi
if compgen -G "$COVERM_DIR/*.count.txt" >/dev/null; then
    paste "$COVERM_DIR"/*.count.txt > "$COVERM_DIR/AD_SSU.count.txt"
fi
if compgen -G "$COVERM_DIR/*.bases.txt" >/dev/null; then
    paste "$COVERM_DIR"/*.bases.txt > "$COVERM_DIR/AD_SSU.bases.txt"
fi

echo "=== Pipeline completed ==="
