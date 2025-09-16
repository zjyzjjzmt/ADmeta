#!/bin/bash
#SBATCH --job-name=viral_analysis
#SBATCH --chdir=/work/$USER
#SBATCH --output=/work/%u/%x-%j.log
#SBATCH --time=0-00:30:00
#SBATCH --cpus-per-task=50
#SBATCH --mem=16G

# 设置错误处理
set -e
set -o pipefail

# 定义路径变量
WORK_DIR="/media/kunyu/data2/03DWTP"
READS_DIR="${WORK_DIR}/01Reads_fq"
CONTIGS_DIR="${WORK_DIR}/06VCs/00Contigs"
OUTPUT_DIR="${WORK_DIR}/06VCs"
DB_DIR="/media/kunyu/data2/db"

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
    fi
}

# 初始化环境
source deactivate
conda activate base
cd "$WORK_DIR"

# 1. 过滤大于5000bp的contigs
echo "Processing contigs > 5000bp..."
activate_conda tools
mkdir -p "${CONTIGS_DIR}"

for forward_read in "${READS_DIR}"/*_1.fastq; do
    reverse_read="${forward_read%_*}_2.fastq"
    check_file "$forward_read"
    check_file "$reverse_read"
    
    basename=$(basename "$forward_read")
    sample=${basename%%_*}

    sed "s/k141/${sample}/g" "05ASSEMBLY/${sample}/final_assembly.fasta" | \
        seqtk seq -L 5000 > "${CONTIGS_DIR}/${sample}.fa"
done

cat "${CONTIGS_DIR}"/*.fa > "${CONTIGS_DIR}/DWTP_contigs5000bp.fasta"
check_file "${CONTIGS_DIR}/DWTP_contigs5000bp.fasta"

# 2. VirFinder分析
echo "Running VirFinder analysis..."
activate_conda vs2
mkdir -p "${OUTPUT_DIR}/02Virfinder"

Rscript "${OUTPUT_DIR}/02Virfinder/virfinder.R" \
    "${CONTIGS_DIR}/Shahe_contigs5000bp.fasta" \
    "${OUTPUT_DIR}/02Virfinder/DWTP.virfinder.tsv"

awk '$4 >= 0.9 && $5 <= 0.01' "${OUTPUT_DIR}/02Virfinder/DWTP.virfinder.tsv" | \
    cut -f 2 | sed 's/"//g' > "${OUTPUT_DIR}/02Virfinder/DWTP.virfinder.txt"

seqtk subseq "${CONTIGS_DIR}/DWTP_contigs5000bp.fasta" \
    "${OUTPUT_DIR}/02Virfinder/DWTP.virfinder.txt" > \
    "${OUTPUT_DIR}/02Virfinder/DWTP.virf.fa"

checkv end_to_end "${OUTPUT_DIR}/02Virfinder/DWTP.virf.fa" \
    "${OUTPUT_DIR}/02Virfinder/checkv" \
    -t 28 -d "${DB_DIR}/checkv/checkv-db-v1.4"

awk '$6 == 0 && $7 > 1 { next } $6 == 0 && $7 == 1 && $2 < 10000 { next } { print }' \
    "${OUTPUT_DIR}/02Virfinder/checkv/quality_summary.tsv" > \
    "${OUTPUT_DIR}/02Virfinder/checkv/final_quality.tsv"

cut -f 1 "${OUTPUT_DIR}/02Virfinder/checkv/final_quality.tsv" > \
    "${OUTPUT_DIR}/02Virfinder/checkv/DWTP_quality_virf.txt"

seqtk subseq "${OUTPUT_DIR}/02Virfinder/DWTP.virf.fa" \
    "${OUTPUT_DIR}/02Virfinder/checkv/DWTP_quality_virf.txt" > \
    "${OUTPUT_DIR}/02Virfinder/DWTP_quality_virf.fasta"

# 3. VirSorter2分析
echo "Running VirSorter2 analysis..."
mkdir -p "${OUTPUT_DIR}/01Virsorter2"
virsorter run --keep-original-seq \
    -i "${CONTIGS_DIR}/DWTP_contigs5000bp.fasta" \
    -w "${OUTPUT_DIR}/01Virsorter2" \
    --include-groups dsDNAphage,ssDNA,NCLDV \
    --min-length 5000 \
    --min-score 0.5 \
    -j 28 all

checkv end_to_end "${OUTPUT_DIR}/01Virsorter2/final-viral-combined.fa" \
    "${OUTPUT_DIR}/01Virsorter2/checkv" \
    -t 28 -d "${DB_DIR}/checkv/checkv-db-v1.4"

awk '$6 == 0 && $7 >= 1 { next } { print }' \
    "${OUTPUT_DIR}/01Virsorter2/checkv/quality_summary.tsv" > \
    "${OUTPUT_DIR}/01Virsorter2/checkv/final_quality.tsv"

cut -f 1 "${OUTPUT_DIR}/01Virsorter2/checkv/final_quality.tsv" > \
    "${OUTPUT_DIR}/01Virsorter2/checkv/DWTP_quality_virs.txt"

seqtk subseq "${OUTPUT_DIR}/01Virsorter2/final-viral-combined.fa" \
    "${OUTPUT_DIR}/01Virsorter2/checkv/DWTP_quality_virs.txt" > \
    "${OUTPUT_DIR}/01Virsorter2/DWTP_quality_virs.fasta"

# 4. geNomad分析
echo "Running geNomad analysis..."
activate_conda genomad
mkdir -p "${OUTPUT_DIR}/01genomad"
genomad end-to-end --cleanup \
    "${CONTIGS_DIR}/DWTP_contigs5000bp.fasta" \
    "${OUTPUT_DIR}/01genomad" \
    "${DB_DIR}/genomad_db/genomad_db"

# 5. 合并结果并去重
echo "Merging and deduplicating results..."
cat "${OUTPUT_DIR}/02Virfinder/DWTP_quality_virf.fasta" \
    "${OUTPUT_DIR}/01Virsorter2/DWTP_quality_virs.fasta" \
    "${OUTPUT_DIR}/01genomad/DWTP_contigs5000bp_summary/DWTP_contigs5000bp_virus.fna" > \
    "${OUTPUT_DIR}/DWTP_quality_virus.fasta"

cd-hit-est -i "${OUTPUT_DIR}/DWTP_quality_virus.fasta" \
    -o "${OUTPUT_DIR}/DWTP_virus_nr.fasta" \
    -c 1.0 -n 10 -M 16000 -T 8

# 6. 聚类分析
echo "Performing clustering analysis..."
cd "${OUTPUT_DIR}"
mkdir -p 01Clustergenome
activate_conda tools

makeblastdb -in DWTP_virus_nr.fasta -dbtype nucl -out 01Clustergenome/DWTP_virus_nr
blastn -query DWTP_virus_nr.fasta \
    -db 01Clustergenome/DWTP_virus_nr \
    -outfmt '6 std qlen slen' \
    -max_target_seqs 10000 \
    -out 01Clustergenome/DWTP_virus_nr.tsv \
    -num_threads 32

python 01Clustergenome/anicalc.py \
    -i 01Clustergenome/DWTP_virus_nr.tsv \
    -o 01Clustergenome/DWTP_virus_ani.tsv

python 01Clustergenome/aniclust.py \
    --fna DWTP_virus_nr.fasta \
    --ani 01Clustergenome/DWTP_virus_ani.tsv \
    --out 01Clustergenome/DWTP_virus_clusters.tsv \
    --min_ani 95 --min_tcov 85 --min_qcov 0

cut -f 1 01Clustergenome/DWTP_virus_clusters.tsv > 01Clustergenome/DWTP_virus_clusters.txt
seqtk subseq DWTP_virus_nr.fasta 01Clustergenome/DWTP_virus_clusters.txt > \
    01Clustergenome/DWTP_VCs.fasta

# 7. 最终CheckV验证
echo "Running final CheckV validation..."
activate_conda vs2
checkv end_to_end 01Clustergenome/DWTP_VCs.fasta \
    03checkv -t 50 -d "${DB_DIR}/checkv/checkv-db-v1.4"

# 8. CoverM分析
echo "Running CoverM analysis..."
cd "${WORK_DIR}"
mkdir -p "${OUTPUT_DIR}/11coverM"
activate_conda tools

for forward_read in "${READS_DIR}"/*_1.fastq; do
    reverse_read="${forward_read%_*}_2.fastq"
    basename=$(basename "$forward_read")
    sample=${basename%%_*}

    coverm make -r "${OUTPUT_DIR}/01Clustergenome/DWTP_VCs.fasta" \
        -1 "$forward_read" -2 "$reverse_read" \
        -o "${OUTPUT_DIR}/11coverM" -t 50

    coverm filter -b "${OUTPUT_DIR}/11coverM/DWTP_VCs.fasta.${basename}.bam" \
        -o "${OUTPUT_DIR}/11coverM/${sample}.bam" \
        --min-read-percent-identity 95 \
        --min-read-aligned-percent 75 \
        --threads 50

    coverm contig --bam-files "${OUTPUT_DIR}/11coverM/${sample}.bam" \
        -o "${OUTPUT_DIR}/11coverM/${sample}.rpkm.txt" \
        --trim-min 10 --trim-max 90 \
        --min-read-percent-identity 95 \
        --min-read-aligned-percent 75 \
        -m rpkm -t 50

    coverm contig --bam-files "${OUTPUT_DIR}/11coverM/${sample}.bam" \
        -o "${OUTPUT_DIR}/11coverM/${sample}.count.txt" \
        --min-covered-fraction 0 -m count -t 50

    coverm contig --bam-files "${OUTPUT_DIR}/11coverM/${sample}.bam" \
        -o "${OUTPUT_DIR}/11coverM/${sample}.bases.txt" \
        --min-covered-fraction 0 -m covered_bases -t 50
done

cd "${OUTPUT_DIR}/11coverM"
paste *.rpkm.txt > DWTP_virus.rpkm.txt
paste *.count.txt > DWTP_virus.count.txt
paste *.bases.txt > DWTP_virus.bases.txt

# 9. VIBRANT分析
echo "Running VIBRANT analysis..."
cd "${OUTPUT_DIR}"
activate_conda vibrant
VIBRANT_run.py -i 01Clustergenome/DWTP_VCs.fasta -f nucl -t 50 -virome

# 10. ARGs分析
echo "Running ARGs analysis..."
activate_conda tools
mkdir -p "${OUTPUT_DIR}/06ARGs"

prodigal -i 01Clustergenome/DWTP_VCs.fasta \
    -o 01Clustergenome/DWTP_VCs.genes \
    -a 01Clustergenome/DWTP_VCs.faa \
    -d 01Clustergenome/DWTP_VCs.fna -p meta

diamond makedb --in "${OUTPUT_DIR}/06ARGs/SARG.2.2.fasta" \
    -d "${OUTPUT_DIR}/06ARGs/SARG.2.2_nr"

diamond blastx -d "${OUTPUT_DIR}/06ARGs/SARG.2.2_nr" \
    -q 01Clustergenome/DWTP_VCs.fna \
    -o "${OUTPUT_DIR}/06ARGs/DWTP_VCs-SARG.txt" \
    --evalue 1e-10 --query-cover 70 --id 80 -k 1

diamond makedb --in "${OUTPUT_DIR}/06ARGs/CARD3.1.4.fasta" \
    -d "${OUTPUT_DIR}/06ARGs/card3.1.4_nr"

diamond blastx -d "${OUTPUT_DIR}/06ARGs/card3.1.4_nr" \
    -q 01Clustergenome/DWTP_VCs.fna \
    -o "${OUTPUT_DIR}/06ARGs/DWTP_VCs-card.txt" \
    --evalue 1e-10 --query-cover 70 --id 80 -k 1

# 11. 后续分析（geNomad分类、CAT、vContact2等）
echo "Running additional analyses..."
# geNomad分类
activate_conda genomad
mkdir -p "${OUTPUT_DIR}/02genomad"
genomad end-to-end --cleanup \
    01Clustergenome/DWTP_VCs.fasta \
    "${OUTPUT_DIR}/02genomad" \
    "${DB_DIR}/genomad_db/genomad_db"

# CAT分类
echo "Running CAT classification..."
activate_conda tools
mkdir -p "${OUTPUT_DIR}/02CAT"
CAT contigs -c 01Clustergenome/DWTP_VCs.fasta \
    -d "${DB_DIR}/CAT_prepare_20210107/2021-01-07_CAT_database" \
    -t "${DB_DIR}/CAT_prepare_20210107/2021-01-07_taxonomy"

CAT add_names -i out.CAT.contig2classification.txt \
    -o "${OUTPUT_DIR}/02CAT/CAT_contigs_classification.txt" \
    -t "${DB_DIR}/CAT_prepare_20210107/2021-01-07_taxonomy"

CAT add_names -i out.CAT.ORF2LCA.txt \
    -o "${OUTPUT_DIR}/02CAT/CAT_ORFs_classification.txt" \
    -t "${DB_DIR}/CAT_prepare_20210107/2021-01-07_taxonomy"

mv out.CAT.* "${OUTPUT_DIR}/02CAT/"
rm -rf out.CAT.alignment.diamond

# vContact2分析
echo "Running vContact2 analysis..."
activate_conda vcontact2
mkdir -p "${OUTPUT_DIR}/04vcontact2"
seqtk seq -L 10000 01Clustergenome/DWTP_VCs.fasta > \
    "${OUTPUT_DIR}/04vcontact2/DWTP_VCs_10kb.fasta"

prodigal -i "${OUTPUT_DIR}/04vcontact2/DWTP_VCs_10kb.fasta" \
    -o "${OUTPUT_DIR}/04vcontact2/DWTP_VCs_10kb.genes" \
    -a "${OUTPUT_DIR}/04vcontact2/DWTP_VCs_10kb.faa" -p meta

vcontact2_gene2genome -p "${OUTPUT_DIR}/04vcontact2/DWTP_VCs_10kb.faa" \
    -o "${OUTPUT_DIR}/04vcontact2/DWTP_VCs_10kb_g2g.csv" \
    -s 'Prodigal-FAA'

vcontact2 -r "${OUTPUT_DIR}/04vcontact2/DWTP_VCs_10kb.faa" \
    --rel-mode 'Diamond' \
    --proteins-fp "${OUTPUT_DIR}/04vcontact2/DWTP_VCs_10kb_g2g.csv" \
    --db 'ProkaryoticViralRefSeq201-Merged' \
    --pcs-mode MCL \
    --vcs-mode ClusterONE \
    --c1-bin /home/kunyu/miniconda3/envs/vcontact2/bin/cluster_one-1.0.jar \
    --output-dir "${OUTPUT_DIR}/04vcontact2/vContactOut/"

# 12. RefSeq和IMG/VR分析
echo "Running RefSeq and IMG/VR analyses..."
activate_conda tools
mkdir -p "${OUTPUT_DIR}/08majority_rules" "${OUTPUT_DIR}/09blast2IMGVR" "${OUTPUT_DIR}/10novelty"

# RefSeq分析
diamond makedb --in "${DB_DIR}/viral_refseq/viral.1.protein.faa" \
    -d "${OUTPUT_DIR}/08majority_rules/viral_refseq"

diamond blastp -q 01Clustergenome/DWTP_VCs.faa \
    -d "${OUTPUT_DIR}/08majority_rules/viral_refseq" \
    -o "${OUTPUT_DIR}/08majority_rules/DWTP_VCs_blastp.txt" \
    --query-cover 50 --subject-cover 50 --evalue 1e-5 -k 1

csvtk filter -t -f "12>=50" "${OUTPUT_DIR}/08majority_rules/DWTP_VCs_blastp.txt" > \
    "${OUTPUT_DIR}/08majority_rules/DWTP_VCs_blastp_50score.txt"

cut -f 2 "${OUTPUT_DIR}/08majority_rules/DWTP_VCs_blastp_50score.txt" > \
    "${OUTPUT_DIR}/08majority_rules/DWTP_VCs_blastp_accession.txt"

mkdir -p "${OUTPUT_DIR}/08majority_rules/02blastp_accession"
split -l 1000 "${OUTPUT_DIR}/08majority_rules/DWTP_VCs_blastp_accession.txt" \
    "${OUTPUT_DIR}/08majority_rules/02blastp_accession/part_"

for file in "${OUTPUT_DIR}/08majority_rules/02blastp_accession/part_"*; do
    rg -f "$file" "${DB_DIR}/protein_taxID/prot.accession2taxid.FULL" \
        --no-line-number >> "${OUTPUT_DIR}/08majority_rules/DWTP_VCs_refseq_accession2taxid.txt"
done

cut -f 2 "${OUTPUT_DIR}/08majority_rules/DWTP_VCs_refseq_accession2taxid.txt" > \
    "${OUTPUT_DIR}/08majority_rules/DWTP_VCs_blastp_taxid.txt"

taxonkit lineage "${OUTPUT_DIR}/08majority_rules/DWTP_VCs_blastp_taxid.txt" \
    --data-dir "${DB_DIR}/taxdump" | \
    taxonkit reformat -F --data-dir "${DB_DIR}/taxdump" | \
    cut -f 1,3 > "${OUTPUT_DIR}/08majority_rules/DWTP_VCs_taxid_taxonomy.txt"

paste "${OUTPUT_DIR}/08majority_rules/DWTP_VCs_refseq_accession2taxid.txt" \
    "${OUTPUT_DIR}/08majority_rules/DWTP_VCs_taxid_taxonomy.txt" | \
    cut -f 1,4 > "${OUTPUT_DIR}/08majority_rules/DWTP_VCs_accession_taxonomy.txt"

awk -F'\t' 'FNR==NR{a[$1]=$2; next}; {if($2 in a) {print $0, "\t"a[$2];} else {print $0, "\tNA"}}' \
    "${OUTPUT_DIR}/08majority_rules/DWTP_VCs_accession_taxonomy.txt" \
    "${OUTPUT_DIR}/08majority_rules/DWTP_VCs_blastp_50score.txt" > \
    "${OUTPUT_DIR}/08majority_rules/DWTP_VCs_blastp_viref_50_tax.txt"

awk -F'\t' -v OFS='\t' '{sub(/_[0-9]+$/, "", $1); split($13, a, ";"); $13=a[1] ";" a[2] ";" a[3] ";" a[4] ";" a[5]; print $1, $13, 1}' \
    "${OUTPUT_DIR}/08majority_rules/DWTP_VCs_blastp_viref_50_tax.txt" > \
    "${OUTPUT_DIR}/08majority_rules/DWTP_VCs_family.txt"

awk -F'\t' '{sums[$1"\t"$2] += $3} END {for (key in sums) print key"\t"sums[key]}' \
    "${OUTPUT_DIR}/08majority_rules/DWTP_VCs_family.txt" > \
    "${OUTPUT_DIR}/08majority_rules/DWTP_VCs_family_stax.txt"

cut -f 1,3 "${OUTPUT_DIR}/08majority_rules/DWTP_VCs_family_stax.txt" | \
    datamash -sW -g1 sum 2 > "${OUTPUT_DIR}/08majority_rules/DWTP_VCs_contigs_stax.txt"

awk 'FNR==NR{a[$1]=$2; next}; {if($1 in a) {print $0, "\t"a[$1];} else {print $0, "\tNA"}}' \
    "${OUTPUT_DIR}/08majority_rules/DWTP_VCs_contigs_stax.txt" \
    "${OUTPUT_DIR}/08majority_rules/DWTP_VCs_family_stax.txt" > \
    "${OUTPUT_DIR}/08majority_rules/DWTP_VCs_family_contigs.txt"

# IMG/VR分析
makeblastdb -dbtype nucl \
    -in "${DB_DIR}/viral_refseq/viral.1.1.genomic.fna" \
    -input_type fasta \
    -title virus_genome \
    -out "${OUTPUT_DIR}/09blast2IMGVR/virus_genome" \
    -blastdb_version 5

blastn -query 01Clustergenome/DWTP_VCs.fasta \
    -out "${OUTPUT_DIR}/09blast2IMGVR/DWTP_VCs_blastnVG.txt" \
    -outfmt "6 qseqid sseqid length qlen slen qstart qend sstart send evalue pident staxids sscinames scomnames sblastnames bitscore salltitles qcovs qcovhsp stitle" \
    -db "${OUTPUT_DIR}/09blast2IMGVR/virus_genome" \
    -dust no -max_target_seqs 1 -perc_identity 95 -evalue 0.00001 -num_threads 20

blastn -query 01Clustergenome/DWTP_VCs.fasta \
    -out "${OUTPUT_DIR}/09blast2IMGVR/DWTP_VCs_blastnIMGVR.txt" \
    -outfmt "6 qseqid sseqid length qlen slen qstart qend sstart send evalue pident staxids sscinames scomnames sblastnames bitscore salltitles qcovs qcovhsp stitle" \
    -db "${DB_DIR}/IMG_VR/IMGVR" \
    -dust no -max_target_seqs 1 -perc_identity 95 -evalue 0.00001 -num_threads 20

cut -f 2 "${OUTPUT_DIR}/09blast2IMGVR/DWTP_VCs_blastnIMGVR.txt" > \
    "${OUTPUT_DIR}/09blast2IMGVR/DWTP_IMGVR_VCblastn.txt"

cat "${DB_DIR}/IMG_VR/IMGVR_all_Host_information-high_confidence.tsv" | \
    csvtk grep -t -P "${OUTPUT_DIR}/09blast2IMGVR/DWTP_IMGVR_VCblastn.txt" > \
    "${OUTPUT_DIR}/09blast2IMGVR/IMGVR_host_results.txt"

# 13. Novelty分析
echo "Running novelty analysis..."
diamond blastp -q 01Clustergenome/DWTP_VCs.faa \
    -d "${DB_DIR}/IMG_VR/IMGVR_proteins" \
    -o "${OUTPUT_DIR}/10novelty/01DWTP_VCs_IMGVR.txt" \
    --id 30 --query-cover 50 --evalue 1e-5 -k 1

diamond blastp -q 01Clustergenome/DWTP_VCs.faa \
    -d "${OUTPUT_DIR}/08majority_rules/viral_refseq" \
    -o "${OUTPUT_DIR}/10novelty/02DWTP_VCs_Refseq.txt" \
    --id 30 --query-cover 50 --evalue 1e-5 -k 1

echo "Pipeline completed successfully!"
