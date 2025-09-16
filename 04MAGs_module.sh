#!/bin/bash
#SBATCH --job-name=AD_mags_analysis
#SBATCH --output=AD_mags_analysis-%j.log
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=50
#SBATCH --mem=64G

# Set strict error handling
set -euo pipefail

# ---------------------- Path Settings ----------------------
WORKDIR="/media/kunyu/data2/03AD"
READS_DIR="${WORKDIR}/01Reads_fq"
MAGS_DIR="${WORKDIR}/06MAGs"
VIRUS_DIR="${WORKDIR}/06VCs"
DB_DIR="/media/kunyu/data2/db"
TMP_DIR="${MAGS_DIR}/tmp"

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
    }
}

ensure_dir() {
    for dir in "$@"; do
        mkdir -p "$dir" || { echo "Error: Failed to create directory: $dir"; exit 1; }
    done
}

# ---------------------- Initialize ----------------------
ensure_dir "${MAGS_DIR}/03beforedRep" "${MAGS_DIR}/04dRep" "${MAGS_DIR}/06gtdbtk_classify" \
           "${MAGS_DIR}/ORFs_genes" "${MAGS_DIR}/ORFs_pro" "${MAGS_DIR}/ORFs_nucl" \
           "${MAGS_DIR}/01MAGs_crisprcas" "${VIRUS_DIR}/10host-link/01crispr" "${VIRUS_DIR}/10host-link/02tRNA" \
           "${VIRUS_DIR}/10host-link/03homology" "${MAGS_DIR}/04dRep/03MAGs_prokka" "${MAGS_DIR}/04dRep/04Defense_system"
activate_conda base
cd "$WORKDIR"

# 1. MAGs Dereplication
echo "Running dRep for MAGs dereplication..."
activate_conda dRep
: > "${MAGS_DIR}/file_list.txt"
ls "${READS_DIR}"/*_1.fastq > "${MAGS_DIR}/file_list.txt"

while read -r F; do
    [[ -f "$F" ]] || continue
    R="${F%_*}_2.fastq"
    check_file "$F"
    check_file "$R"
    SAMPLE=$(basename "${F%_*}")

    BIN_DIR="${MAGS_DIR}/02refine_bins/${SAMPLE}/metawrap_50_20_bins"
    if [[ -d "$BIN_DIR" ]]; then
        cd "$BIN_DIR"
        for bin_file in *.fa; do
            [[ -f "$bin_file" ]] && mv "$bin_file" "${SAMPLE}-${bin_file}"
        done
        cd "$WORKDIR"
        cp "${BIN_DIR}"/*.fa "${MAGS_DIR}/03beforedRep/"
    else
        echo "Warning: Directory $BIN_DIR not found for sample $SAMPLE"
    fi
done < "${MAGS_DIR}/file_list.txt"

dRep dereplicate "${MAGS_DIR}/04dRep" -g "${MAGS_DIR}/03beforedRep"/*.fa -comp 70 -con 10 -p 50 --S_algorithm gANI

# 2. GTDB-Tk Classification
echo "Running GTDB-Tk classification..."
activate_conda gtdbtk
gtdbtk classify_wf --skip_ani_screen --genome_dir "${MAGS_DIR}/04dRep/dereplicated_genomes" \
    -x fa --out_dir "${MAGS_DIR}/06gtdbtk_classify" --cpus 50 --mash_db "${MAGS_DIR}/06gtdbtk_classify/gtdbtk.msh"

# 3. ARGs Analysis
echo "Running ARGs analysis on MAGs..."
activate_conda tools
ls "${MAGS_DIR}/04dRep/dereplicated_genomes"/*.fa > "${MAGS_DIR}/file_list1.txt"

while read -r genome; do
    [[ -f "$genome" ]] || continue
    SAMPLE=$(basename "${genome%.*}")
    prodigal -i "$genome" -o "${MAGS_DIR}/ORFs_genes/${SAMPLE}.genes" \
        -a "${MAGS_DIR}/ORFs_pro/${SAMPLE}.pro.fa" -d "${MAGS_DIR}/ORFs_nucl/${SAMPLE}.fa" -p meta
done < "${MAGS_DIR}/file_list1.txt"

# Create Diamond/BLAST databases
declare -A databases=(
    ["01SARG/SARG.2.2.fasta"]="SARG.2.2_nr"
    ["02card/CARD3.1.4.fasta"]="card3.1.4_nr"
    ["03bacmet/BacMet2.fasta"]="bacmet2_nr"
    ["04ICEs/ICE3.0.fasta"]="ICE3.0_nr"
    ["06victors/victors_pro.fasta"]="victors_nr"
    ["05MGEs/MGEs.fasta"]="MGEs"
)

for db_path in "${!databases[@]}"; do
    ensure_dir "${MAGS_DIR}/ORFs_nucl/${db_path%/*}"
    check_file "${MAGS_DIR}/ORFs_nucl/${db_path}"
    if [[ "${databases[$db_path]}" == "MGEs" ]]; then
        makeblastdb -dbtype nucl -in "${MAGS_DIR}/ORFs_nucl/${db_path}" -input_type fasta \
            -title MGEs -out "${MAGS_DIR}/ORFs_nucl/${db_path%/*}/${databases[$db_path]}"
    else
        diamond makedb --in "${MAGS_DIR}/ORFs_nucl/${db_path}" -d "${MAGS_DIR}/ORFs_nucl/${db_path%/*}/${databases[$db_path]}" --threads 10
    fi
done

# Run Diamond/BLAST
ls "${MAGS_DIR}/ORFs_nucl"/*.fa > "${MAGS_DIR}/file_list2.txt"
while read -r query; do
    [[ -f "$query" ]] || continue
    SAMPLE=$(basename "${query%.*}")
    for db_type in "${!databases[@]}"; do
        db_name=${databases[$db_type]}
        out_dir="${MAGS_DIR}/ORFs_nucl/${db_type%/*}"
        if [[ "$db_name" == "MGEs" ]]; then
            blastn -query "$query" -db "${out_dir}/${db_name}" \
                -out "${out_dir}/${SAMPLE}-MGEs.txt" -evalue 1e-10 -perc_identity 80 -num_threads 10 -outfmt 6 -max_target_seqs 1
        else
            diamond blastx -d "${out_dir}/${db_name}" -q "$query" -o "${out_dir}/${SAMPLE}-${db_name%.nr}.txt" \
                --evalue 1e-10 --query-cover 70 --id 80 -k 1 --threads 10
        fi
    done
done < "${MAGS_DIR}/file_list2.txt"

# 4. Phage-Host Linkage (CRISPR)
echo "Running CRISPR analysis for phage-host linkage..."
cd "${MAGS_DIR}/04dRep"
CRT_JAR="CRT1.2-CLI.jar"
check_file "$CRT_JAR"
for genome in dereplicated_genomes/*.fa; do
    [[ -f "$genome" ]] || continue
    SAMPLE=$(basename "${genome%.*}")
    java -cp "$CRT_JAR" crt "$genome" "${MAGS_DIR}/01MAGs_crisprcas/${SAMPLE}.out"
    sed -i "s/k141\|k119/${SAMPLE}/g" "${MAGS_DIR}/01MAGs_crisprcas/${SAMPLE}.out"
done

cd "$MAGS_DIR"
check_file "MAGs_CRT_spacers.R"
Rscript MAGs_CRT_spacers.R
cd 03MAGs_spacers
cat *.txt > MAGs_spacers.fasta
seqtk seq -L 1 MAGs_spacers.fasta > "${VIRUS_DIR}/10host-link/01crispr/MAGs_spacers_nr.fasta"

# 5. Virus tRNA Analysis
echo "Running virus tRNA analysis..."
cd "${VIRUS_DIR}"
aragorn -t -fasta -wa -o 10host-link/02tRNA/AD_VCs_tRNA.fasta 01Clustergenome/AD_VCs.fasta
sed -e '/gene/d' -e '/      /d' 10host-link/02tRNA/AD_VCs_tRNA.fasta | \
    seqtk seq -N -L 1 | sed 's/ c\[/c\[/g; s/)\ \[/)\[/g' > 10host-link/02tRNA/AD_VCs_tRNA_nr.fasta

# 6. MAGs tRNA Database
echo "Building MAGs tRNA database..."
cd "${MAGS_DIR}/04dRep"
ensure_dir 02MAGs_database
for genome in dereplicated_genomes/*.fa; do
    [[ -f "$genome" ]] || continue
    SAMPLE=$(basename "${genome%.*}")
    sed "s/k141\|k119/${SAMPLE}/g" "$genome" > "02MAGs_database/${SAMPLE}.fa"
done

cat 02MAGs_database/*.fa > "${VIRUS_DIR}/10host-link/02tRNA/AD_MAGs_database.fa"
cd "${VIRUS_DIR}"
makeblastdb -dbtype nucl -in 10host-link/02tRNA/AD_MAGs_database.fa -input_type fasta \
    -title AD_MAGs -out 10host-link/02tRNA/AD_MAGs -blastdb_version 5

blastn -query 10host-link/02tRNA/AD_VCs_tRNA_nr.fasta -out 10host-link/02tRNA/AD_MAGs_hosts_tRNA.txt \
    -outfmt "6 qseqid sseqid length qlen slen qstart qend sstart send evalue pident staxids bitscore salltitles qcovs qcovhsp stitle" \
    -db 10host-link/02tRNA/AD_MAGs -dust no -perc_identity 100 -evalue 0.0001 -num_threads 20

makeblastdb -dbtype nucl -in 01Clustergenome/AD_VCs.fasta -input_type fasta \
    -title AD_VCs -out 10host-link/01crispr/AD_VCs -blastdb_version 5

blastn -task blastn-short -query 10host-link/01crispr/MAGs_spacers_nr.fasta -db 10host-link/01crispr/AD_VCs \
    -perc_identity 97 -out 10host-link/01crispr/AD_VCs_MAGs_criprcas_host.txt \
    -outfmt "6 qseqid sseqid length qlen slen qstart qend sstart send evalue pident bitscore mismatch qcovs qcovhsp" \
    -max_target_seqs 1 -num_threads 20

# 7. Homology Analysis
echo "Running homology analysis..."
cp 10host-link/02tRNA/AD_MAGs_database.fa 10host-link/03homology/
makeblastdb -dbtype nucl -in 10host-link/03homology/AD_MAGs_database.fa -input_type fasta \
    -title AD_MAGs -out 10host-link/03homology/AD_MAGs -blastdb_version 5

blastn -query 01Clustergenome/AD_VCs.fasta -out 10host-link/03homology/AD_Virus_homology.txt \
    -outfmt "6 qseqid sseqid length qlen slen qstart qend sstart send evalue pident staxids sscinames scomnames sblastnames bitscore salltitles qcovs qcovhsp stitle" \
    -db 10host-link/03homology/AD_MAGs -dust no -perc_identity 80 -evalue 0.00001 -num_threads 40

# 8. Antiviral Defense Systems
echo "Running antiviral defense system analysis..."
cd "${MAGS_DIR}/04dRep"
ensure_dir 03MAGs_prokka/{01err,02faa,03ffn,04fna,05fsa,06gbk,07gff,08log,09sqn,10tbl,11tsv,12txt}
for genome in dereplicated_genomes/*.fa; do
    [[ -f "$genome" ]] || continue
    SAMPLE=$(basename "${genome%.*}")
    TMPDIR="$TMP_DIR" prokka --outdir 03MAGs_prokka --cpus 40 --force --prefix "$SAMPLE" "$genome"
    for ext in err faa ffn fna fsa gbk gff log sqn tbl tsv txt; do
        [[ -f "03MAGs_prokka/${SAMPLE}.${ext}" ]] && mv "03MAGs_prokka/${SAMPLE}.${ext}" "03MAGs_prokka/0${ext#err}/${SAMPLE}.${ext}"
    done
done

# PADLOC
echo "Running PADLOC analysis..."
for genome in dereplicated_genomes/*.fa; do
    [[ -f "$genome" ]] || continue
    SAMPLE=$(basename "${genome%.*}")
    padloc --faa "03MAGs_prokka/02faa/${SAMPLE}.faa" --gff "03MAGs_prokka/07gff/${SAMPLE}.gff" \
        --outdir 04Defense_system --cpu 40
done

cd 04Defense_system
for file in *_padloc.csv; do
    [[ -f "$file" ]] || continue
    replace_str=${file%_padloc.csv}
    sed -i "s/k141/${replace_str}/g" "$file"
done
cat *_padloc.csv > MAGs_padloc_results.csv

# DefenseFinder
echo "Running DefenseFinder analysis..."
activate_conda defensefinder
for pro_file in "${MAGS_DIR}/ORFs_pro"/*.fa; do
    [[ -f "$pro_file" ]] || continue
    SAMPLE=$(basename "${pro_file%.*}")
    defense-finder run "$pro_file" --out-dir "${MAGS_DIR}/04dRep/04Defense_system" --threads 40
done

cat "${MAGS_DIR}/04dRep/04Defense_system"/*defense_finder_systems.tsv > \
    "${MAGS_DIR}/04dRep/04Defense_system/MAGs_defense_finder_results.tsv"

echo "MAGs module completed successfully!"