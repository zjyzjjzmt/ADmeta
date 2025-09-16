#!/bin/bash
#SBATCH --job-name=mags_analysis
#SBATCH --chdir=/work/$USER
#SBATCH --output=/work/%u/%x-%j.log
#SBATCH --time=0-01:00:00
#SBATCH --cpus-per-task=50
#SBATCH --mem=32G

# error
set -e
set -o pipefail

# pathway
WORK_DIR="/media/kunyu/data2/03DWTP"
READS_DIR="/media/kunyu/data7/01Riverine_left/01Reads_fq"
MAGS_DIR="${WORK_DIR}/06MAGs"
DB_DIR="/media/kunyu/data2/db"
TMP_DIR="${MAGS_DIR}/tmp"

# 检查工作目录是否存在
if [ ! -d "$WORK_DIR" ]; then
    echo "Error: Working directory $WORK_DIR does not exist"
    exit 1
fi

# 函数：激活conda环境
activate_conda() {
    local env=$1
    conda deactivate
    conda activate "$env" || { echo "Failed to activate conda env: $env"; exit 1; }
}

# 函数：检查文件是否存在
check_file() {
    local file=$1
    if [ ! -f "$file" ]; then
        echo "Error: File $file not found"
        exit 1
    }
}

# 初始化环境
source deactivate
conda activate base
cd "$WORK_DIR"

# 1. MAGs去重（dRep）
echo "Running dRep for MAGs dereplication..."
activate_conda dRep
mkdir -p "${MAGS_DIR}/03beforedRep" "${MAGS_DIR}/04dRep"

# 生成输入文件列表
ls "${READS_DIR}"/*.fastq.gz > "${MAGS_DIR}/file_list.txt"

while read -r forward_read; do
    reverse_read="${forward_read%_*}_2.fastq.gz"
    check_file "$forward_read"
    check_file "$reverse_read"
    
    basename=$(basename "$forward_read")
    sample=${basename%%_*}

    # 重命名bin文件
    bin_dir="${MAGS_DIR}/02refine_bins/${sample}/metawrap_50_20_bins"
    if [ -d "$bin_dir" ]; then
        cd "$bin_dir"
        rename "s/^/${sample}-/g" *.fa
        cd "$WORK_DIR"
        cp "${bin_dir}"/*.fa "${MAGS_DIR}/03beforedRep/"
    else
        echo "Warning: Directory $bin_dir not found for sample $sample"
    fi
done < "${MAGS_DIR}/file_list.txt"

dRep dereplicate "${MAGS_DIR}/04dRep" \
    -g "${MAGS_DIR}/03beforedRep"/*.fa \
    -comp 70 -con 10 -p 50

# 2. GTDB-Tk分类
echo "Running GTDB-Tk classification..."
activate_conda gtdbtk
mkdir -p "${MAGS_DIR}/06gtdbtk_classify"
gtdbtk classify_wf --skip_ani_screen \
    --genome_dir "${MAGS_DIR}/04dRep/dereplicated_genomes" \
    -x fa \
    --out_dir "${MAGS_DIR}/06gtdbtk_classify" \
    --cpus 50

# 3. MAGs ARGs分析
echo "Running ARGs analysis..."
activate_conda tools
mkdir -p "${MAGS_DIR}/ORFs_genes" "${MAGS_DIR}/ORFs_pro" "${MAGS_DIR}/ORFs_nucl"

# 生成基因预测
ls "${MAGS_DIR}/04dRep/dereplicated_genomes"/*.fa > "${MAGS_DIR}/file_list1.txt"
while read -r genome; do
    basename=$(basename "$genome")
    sample=${basename%%.*}
    
    prodigal -i "$genome" \
        -o "${MAGS_DIR}/ORFs_genes/${sample}.genes" \
        -a "${MAGS_DIR}/ORFs_pro/${sample}.pro.fa" \
        -d "${MAGS_DIR}/ORFs_nucl/${sample}.fa" \
        -p meta
done < "${MAGS_DIR}/file_list1.txt"

# 创建数据库
declare -A databases=(
    ["01SARG/SARG.2.2.fasta"]="SARG.2.2_nr"
    ["02card/CARD3.1.4.fasta"]="card3.1.4_nr"
    ["03bacmet/BacMet2.fasta"]="bacmet2_nr"
    ["04ICEs/ICE3.0.fasta"]="ICE3.0_nr"
    ["06victors/victors_pro.fasta"]="victors_nr"
)

for db_path in "${!databases[@]}"; do
    mkdir -p "${MAGS_DIR}/ORFs_nucl/$(dirname "$db_path")"
    diamond makedb --in "${MAGS_DIR}/ORFs_nucl/${db_path}" \
        -d "${MAGS_DIR}/ORFs_nucl/${db_path%/*}/${databases[$db_path]}"
done

# 运行DIAMOND比对
ls "${MAGS_DIR}/ORFs_nucl"/*.fa > "${MAGS_DIR}/file_list2.txt"
while read -r query; do
    basename=$(basename "$query")
    sample=${basename%%.*}
    
    declare -A blast_params=(
        ["SARG"]="--evalue 1e-10"
        ["card"]="--evalue 1e-5"
        ["bacmet2"]="--evalue 1e-10"
        ["ICE"]="--evalue 1e-10"
        ["victors"]="--evalue 1e-10"
    )
    
    for db_type in "${!blast_params[@]}"; do
        diamond blastx -d "${MAGS_DIR}/ORFs_nucl/${db_type%2*}/${db_type}_nr" \
            -q "$query" \
            -o "${MAGS_DIR}/ORFs_nucl/${db_type%2*}/${sample}-${db_type}.txt" \
            ${blast_params[$db_type]} --query-cover 70 --id 80 -k 1
    done
done < "${MAGS_DIR}/file_list2.txt"

# MGEs分析
makeblastdb -dbtype nucl \
    -in "${MAGS_DIR}/ORFs_nucl/05MGEs/MGEs.fasta" \
    -input_type fasta \
    -title MGEs \
    -out "${MAGS_DIR}/ORFs_nucl/05MGEs/MGEs"

while read -r query; do
    basename=$(basename "$query")
    sample=${basename%%.*}
    
    blastn -query "$query" \
        -db "${MAGS_DIR}/ORFs_nucl/05MGEs/MGEs" \
        -out "${MAGS_DIR}/ORFs_nucl/05MGEs/${sample}-MGEs.txt" \
        -evalue 1e-10 -perc_identity 80 -num_threads 10 -outfmt 6 -max_target_seqs 1
done < "${MAGS_DIR}/file_list2.txt"

# 4. Phage-host linkage分析
echo "Running phage-host linkage analysis..."
cd "${MAGS_DIR}/04dRep"

# CRISPR分析
mkdir -p "${MAGS_DIR}/01MAGs_crisprcas" "${WORK_DIR}/06VCs/10host-link/01crispr"
for genome in dereplicated_genomes/*.fa; do
    basename=$(basename "$genome")
    sample=${basename%%.*}
    
    java -cp CRT1.2-CLI.jar crt "$genome" "${MAGS_DIR}/01MAGs_crisprcas/${sample}.out"
    sed -i "s/k141\|k119/${sample}/g" "${MAGS_DIR}/01MAGs_crisprcas/${sample}.out"
done

cd "${MAGS_DIR}"
Rscript MAGs_CRT_spacers.R
cd 03MAGs_spacers
cat *.txt > MAGs_spacers.fasta
seqtk seq -L 1 MAGs_spacers.fasta > "${WORK_DIR}/06VCs/10host-link/01crispr/MAGs_spacers_nr.fasta"

# 5. Virus tRNA分析
echo "Running virus tRNA analysis..."
cd "${WORK_DIR}/06VCs"
mkdir -p 10host-link/02tRNA
aragorn -t -fasta -wa \
    -o 10host-link/02tRNA/DWTP_VCs_tRNA.fasta \
    01Clustergenome/DWTP_VCs.fasta

sed -e '/gene/d' 10host-link/02tRNA/DWTP_VCs_tRNA.fasta | \
    sed -e '/      /d' | \
    seqtk seq -N -L 1 > 10host-link/02tRNA/DWTP_VCs_tRNA_nr.fasta

sed -i 's/ c\[/c\[/g; s/)\ \[/)\[/g' 10host-link/02tRNA/DWTP_VCs_tRNA_nr.fasta

# 6. MAGs tRNA数据库
echo "Building MAGs tRNA database..."
cd "${MAGS_DIR}/04dRep"
mkdir -p 02MAGs_database
for genome in dereplicated_genomes/*.fa; do
    basename=$(basename "$genome")
    sample=${basename%%.*}
    sed "s/k141\|k119/${sample}/g" "$genome" > "02MAGs_database/${sample}.fa"
done

cat 02MAGs_database/*.fa > "${WORK_DIR}/06VCs/10host-link/02tRNA/DWTP_MAGs_database.fa"

cd "${WORK_DIR}/06VCs"
makeblastdb -dbtype nucl \
    -in 10host-link/02tRNA/DWTP_MAGs_database.fa \
    -input_type fasta \
    -title DWTP_MAGs \
    -out 10host-link/02tRNA/DWTP_MAGs \
    -blastdb_version 5

blastn -query 10host-link/02tRNA/DWTP_VCs_tRNA_nr.fasta \
    -out 10host-link/02tRNA/DWTP_MAGs_hosts_tRNA.txt \
    -outfmt "6 qseqid sseqid length qlen slen qstart qend sstart send evalue pident staxids bitscore salltitles qcovs qcovhsp stitle" \
    -db 10host-link/02tRNA/DWTP_MAGs \
    -dust no -perc_identity 100 -evalue 0.0001 -num_threads 20

makeblastdb -dbtype nucl \
    -in 01Clustergenome/DWTP_VCs.fasta \
    -input_type fasta \
    -title DWTP_VCs \
    -out 10host-link/01crispr/DWTP_VCs \
    -blastdb_version 5

blastn -task blastn-short \
    -query 10host-link/01crispr/MAGs_spacers_nr.fasta \
    -db 10host-link/01crispr/DWTP_VCs \
    -perc_identity 97 \
    -out 10host-link/01crispr/DWTP_VCs_MAGs_criprcas_host.txt \
    -outfmt "6 qseqid sseqid length qlen slen qstart qend sstart send evalue pident bitscore mismatch qcovs qcovhsp" \
    -max_target_seqs 1 -num_threads 20

# 7. Homology分析
echo "Running homology analysis..."
cp 10host-link/02tRNA/DWTP_MAGs_database.fa 10host-link/03homology/
makeblastdb -dbtype nucl \
    -in 10host-link/03homology/DWTP_MAGs_database.fa \
    -input_type fasta \
    -title DWTP_MAGs \
    -out 10host-link/03homology/DWTP_MAGs \
    -blastdb_version 5

blastn -query 01Clustergenome/DWTP_VCs.fasta \
    -out 10host-link/03homology/DWTP_Virus_homology.txt \
    -outfmt "6 qseqid sseqid length qlen slen qstart qend sstart send evalue pident staxids sscinames scomnames sblastnames bitscore salltitles qcovs qcovhsp stitle" \
    -db 10host-link/03homology/DWTP_MAGs \
    -dust no -perc_identity 80 -evalue 0.00001 -num_threads 40

# 8. 抗病毒防御系统分析
echo "Running antiviral defense system analysis..."
cd "${MAGS_DIR}/04dRep"
mkdir -p 03MAGs_prokka/{01err,02faa,03ffn,04fna,05fsa,06gbk,07gff,08log,09sqn,10tbl,11tsv,12txt}

for genome in dereplicated_genomes/*.fa; do
    basename=$(basename "$genome")
    sample=${basename%%.*}
    
    TMPDIR="$TMP_DIR" prokka --outdir 03MAGs_prokka \
        --cpus 40 --force --prefix "$sample" "$genome"
    
    mv 03MAGs_prokka/${sample}.* 03MAGs_prokka/$(echo ${sample}.* | cut -d. -f2)/
done

# PADLOC分析
mkdir -p 04Defense_system
for genome in dereplicated_genomes/*.fa; do
    basename=$(basename "$genome")
    sample=${basename%%.*}
    
    padloc --faa "03MAGs_prokka/02faa/${sample}.faa" \
        --gff "03MAGs_prokka/07gff/${sample}.gff" \
        --outdir 04Defense_system --cpu 40
done

cd 04Defense_system
for file in *_padloc.csv; do
    replace_str=${file%_padloc.csv}
    sed -i "s/k141/${replace_str}/g" "$file"
done
cat *_padloc.csv > MAGs_padloc_results.csv

# DefenseFinder分析
activate_conda defensefinder
for pro_file in "${MAGS_DIR}/ORFs_pro"/*.fa; do
    basename=$(basename "$pro_file")
    sample=${basename%%.*}
    
    defense-finder run "$pro_file" --out-dir "${MAGS_DIR}/04dRep/04Defense_system"
done

cat "${MAGS_DIR}/04dRep/04Defense_system"/*defense_finder_systems.tsv > \
    "${MAGS_DIR}/04dRep/04Defense_system/MAGs_defense_finder_results.tsv"

echo "Pipeline completed successfully!"