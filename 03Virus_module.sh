#!/bin/bash
#SBATCH --job-name=AD_virus_analysis
#SBATCH --output=AD_virus_analysis-%j.log
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=50
#SBATCH --mem=64G

# Set strict error handling
set -euo pipefail

# ---------------------- Path Settings ----------------------
WORKDIR="/media/kunyu/data2/03AD"
READS_DIR="${WORKDIR}/01Reads_fq"
CONTIGS_DIR="${WORKDIR}/06VCs/00Contigs"
OUTPUT_DIR="${WORKDIR}/06VCs"
DB_DIR="/media/kunyu/data2/db"

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
ensure_dir "$CONTIGS_DIR" "${OUTPUT_DIR}/01Virsorter2" "${OUTPUT_DIR}/02Virfinder" "${OUTPUT_DIR}/01genomad" \
           "${OUTPUT_DIR}/02CAT" "${OUTPUT_DIR}/04vcontact2" "${OUTPUT_DIR}/06ARGs" "${OUTPUT_DIR}/08majority_rules" \
           "${OUTPUT_DIR}/09blast2IMGVR" "${OUTPUT_DIR}/10novelty"
activate_conda base
cd "$WORKDIR"

# 1. Filter Contigs > 5000bp
echo "Filtering contigs > 5000bp..."
activate_conda tools
for F in "${READS_DIR}"/*_1.fastq; do
    [[ -f "$F" ]] || continue
    R="${F%_*}_2.fastq"
    check_file "$F"
    check_file "$R"
    SAMPLE=$(basename "${F%_*}")

    ASSEMBLY="${WORKDIR}/05ASSEMBLY/${SAMPLE}/final_assembly.fasta"
    check_file "$ASSEMBLY"

    sed "s/k141/${SAMPLE}/g" "$ASSEMBLY" | seqtk seq -L 5000 > "${CONTIGS_DIR}/${SAMPLE}.fa"
done

cat "${CONTIGS_DIR}"/*.fa > "${CONTIGS_DIR}/AD_contigs5000bp.fasta"
check_file "${CONTIGS_DIR}/AD_contigs5000bp.fasta"

# 2. VirFinder Analysis
echo "Running VirFinder analysis..."
activate_conda vs2
check_file "${CONTIGS_DIR}/AD_contigs5000bp.fasta" # Fixed Shahe to AD
Rscript "${OUTPUT_DIR}/02Virfinder/virfinder.R" \
    "${CONTIGS_DIR}/AD_contigs5000bp.fasta" \
    "${OUTPUT_DIR}/02Virfinder/AD.virfinder.tsv"

awk '$4 >= 0.9 && $5 <= 0.01' "${OUTPUT_DIR}/02Virfinder/AD.virfinder.tsv" | \
    cut -f 2 | sed 's/"//g' > "${OUTPUT_DIR}/02Virfinder/AD.virfinder.txt"

seqtk subseq "${CONTIGS_DIR}/AD_contigs5000bp.fasta" \
    "${OUTPUT_DIR}/02Virfinder/AD.virfinder.txt" > "${OUTPUT_DIR}/02Virfinder/AD.virf.fa"

checkv end_to_end "${OUTPUT_DIR}/02Virfinder/AD.virf.fa" "${OUTPUT_DIR}/02Virfinder/checkv" \
    -t 28 -d "${DB_DIR}/checkv/checkv-db-v1.4"

awk '$6 == 0 && $7 > 1 { next } $6 == 0 && $7 == 1 && $2 < 10000 { next } { print }' \
    "${OUTPUT_DIR}/02Virfinder/checkv/quality_summary.tsv" > \
    "${OUTPUT_DIR}/02Virfinder/checkv/final_quality.tsv"

cut -f 1 "${OUTPUT_DIR}/02Virfinder/checkv/final_quality.tsv" > \
    "${OUTPUT_DIR}/02Virfinder/checkv/AD_quality_virf.txt"

seqtk subseq "${OUTPUT_DIR}/02Virfinder/AD.virf.fa" \
    "${OUTPUT_DIR}/02Virfinder/checkv/AD_quality_virf.txt" > "${OUTPUT_DIR}/02Virfinder/AD_quality_virf.fasta"

# 3. VirSorter2 Analysis
echo "Running VirSorter2 analysis..."
virsorter run --keep-original-seq -i "${CONTIGS_DIR}/AD_contigs5000bp.fasta" \
    -w "${OUTPUT_DIR}/01Virsorter2" --include-groups dsDNAphage,ssDNA,NCLDV \
    --min-length 5000 --min-score 0.5 -j 28 all

checkv end_to_end "${OUTPUT_DIR}/01Virsorter2/final-viral-combined.fa" \
    "${OUTPUT_DIR}/01Virsorter2/checkv" -t 28 -d "${DB_DIR}/checkv/checkv-db-v1.4"

awk '$6 == 0 && $7 >= 1 { next } { print }' "${OUTPUT_DIR}/01Virsorter2/checkv/quality_summary.tsv" > \
    "${OUTPUT_DIR}/01Virsorter2/checkv/final_quality.tsv"

cut -f 1 "${OUTPUT_DIR}/01Virsorter2/checkv/final_quality.tsv" > \
    "${OUTPUT_DIR}/01Virsorter2/checkv/AD_quality_virs.txt"

seqtk subseq "${OUTPUT_DIR}/01Virsorter2/final-viral-combined.fa" \
    "${OUTPUT_DIR}/01Virsorter2/checkv/AD_quality_virs.txt" > "${OUTPUT_DIR}/01Virsorter2/AD_quality_virs.fasta"

# 4. geNomad Analysis
echo "Running geNomad analysis..."
activate_conda genomad
genomad end-to-end --cleanup "${CONTIGS_DIR}/AD_contigs5000bp.fasta" \
    "${OUTPUT_DIR}/01genomad" "${DB_DIR}/genomad_db/genomad_db"

# 5. Merge and Deduplicate Viral Sequences
echo "Merging and deduplicating viral sequences..."
cat "${OUTPUT_DIR}/02Virfinder/AD_quality_virf.fasta" \
    "${OUTPUT_DIR}/01Virsorter2/AD_quality_virs.fasta" \
    "${OUTPUT_DIR}/01genomad/AD_contigs5000bp_summary/AD_contigs5000bp_virus.fna" > \
    "${OUTPUT_DIR}/AD_quality_virus.fasta"

cd-hit-est -i "${OUTPUT_DIR}/AD_quality_virus.fasta" -o "${OUTPUT_DIR}/AD_virus_nr.fasta" \
    -c 1.0 -n 10 -M 16000 -T 8

# 6. Viral Clustering
echo "Performing viral clustering..."
cd "${OUTPUT_DIR}"
ensure_dir 01Clustergenome
makeblastdb -in AD_virus_nr.fasta -dbtype nucl -out 01Clustergenome/AD_virus_nr
blastn -query AD_virus_nr.fasta -db 01Clustergenome/AD_virus_nr \
    -outfmt '6 std qlen slen' -max_target_seqs 10000 -out 01Clustergenome/AD_virus_nr.tsv -num_threads 32

python 01Clustergenome/anicalc.py -i 01Clustergenome/AD_virus_nr.tsv -o 01Clustergenome/AD_virus_ani.tsv
python 01Clustergenome/aniclust.py --fna AD_virus_nr.fasta --ani 01Clustergenome/AD_virus_ani.tsv \
    --out 01Clustergenome/AD_virus_clusters.tsv --min_ani 95 --min_tcov 85 --min_qcov 0

cut -f 1 01Clustergenome/AD_virus_clusters.tsv > 01Clustergenome/AD_virus_clusters.txt
seqtk subseq AD_virus_nr.fasta 01Clustergenome/AD_virus_clusters.txt > 01Clustergenome/AD_VCs.fasta

# 7. Final CheckV Validation
echo "Running final CheckV validation..."
activate_conda vs2
checkv end_to_end 01Clustergenome/AD_VCs.fasta "${OUTPUT_DIR}/03checkv" -t 50 -d "${DB_DIR}/checkv/checkv-db-v1.4"

# 8. CoverM Analysis
echo "Running CoverM analysis..."
activate_conda tools
ensure_dir "${OUTPUT_DIR}/11coverM"
for F in "${READS_DIR}"/*_1.fastq; do
    [[ -f "$F" ]] || continue
    R="${F%_*}_2.fastq"
    check_file "$R"
    SAMPLE=$(basename "${F%_*}")

    echo "Mapping reads for sample: $SAMPLE"
    coverm make -r "${OUTPUT_DIR}/01Clustergenome/AD_VCs.fasta" -1 "$F" -2 "$R" \
        -o "${OUTPUT_DIR}/11coverM" -t 50

    coverm filter -b "${OUTPUT_DIR}/11coverM/AD_VCs.fasta.${SAMPLE}_1.fastq.bam" \
        -o "${OUTPUT_DIR}/11coverM/${SAMPLE}.bam" --min-read-percent-identity 95 \
        --min-read-aligned-percent 75 --threads 50

    coverm contig --bam-files "${OUTPUT_DIR}/11coverM/${SAMPLE}.bam" \
        -o "${OUTPUT_DIR}/11coverM/${SAMPLE}.rpkm.txt" --trim-min 10 --trim-max 90 \
        --min-read-percent-identity 95 --min-read-aligned-percent 75 -m rpkm -t 50
    coverm contig --bam-files "${OUTPUT_DIR}/11coverM/${SAMPLE}.bam" \
        -o "${OUTPUT_DIR}/11coverM/${SAMPLE}.count.txt" --min-covered-fraction 0 -m count -t 50
    coverm contig --bam-files "${OUTPUT_DIR}/11coverM/${SAMPLE}.bam" \
        -o "${OUTPUT_DIR}/11coverM/${SAMPLE}.bases.txt" --min-covered-fraction 0 -m covered_bases -t 50
done

cd "${OUTPUT_DIR}/11coverM"
paste *.rpkm.txt > AD_virus.rpkm.txt
paste *.count.txt > AD_virus.count.txt
paste *.bases.txt > AD_virus.bases.txt

# 9. VIBRANT Analysis
echo "Running VIBRANT analysis..."
activate_conda vibrant
VIBRANT_run.py -i 01Clustergenome/AD_VCs.fasta -f nucl -t 50 -virome

# 10. ARGs Analysis
echo "Running ARGs analysis on viral sequences..."
activate_conda tools
ensure_dir "${OUTPUT_DIR}/06ARGs"
prodigal -i 01Clustergenome/AD_VCs.fasta -o 01Clustergenome/AD_VCs.genes \
    -a 01Clustergenome/AD_VCs.faa -d 01Clustergenome/AD_VCs.fna -p meta

for db in SARG.2.2 card3.1.4; do
    check_file "${OUTPUT_DIR}/06ARGs/${db}.fasta"
    diamond makedb --in "${OUTPUT_DIR}/06ARGs/${db}.fasta" -d "${OUTPUT_DIR}/06ARGs/${db}_nr" --threads 10
    diamond blastx -d "${OUTPUT_DIR}/06ARGs/${db}_nr" -q 01Clustergenome/AD_VCs.fna \
        -o "${OUTPUT_DIR}/06ARGs/AD_VCs-${db}.txt" --evalue 1e-10 --query-cover 70 --id 80 -k 1 --threads 10
done

# 11. Additional Analyses (geNomad, CAT, vContact2)
echo "Running additional analyses..."
# geNomad
activate_conda genomad
genomad end-to-end --cleanup 01Clustergenome/AD_VCs.fasta "${OUTPUT_DIR}/02genomad" "${DB_DIR}/genomad_db/genomad_db"

# CAT
echo "Running CAT classification..."
activate_conda tools
CAT contigs -c 01Clustergenome/AD_VCs.fasta \
    -d "${DB_DIR}/CAT_prepare_20210107/2021-01-07_CAT_database" \
    -t "${DB_DIR}/CAT_prepare_20210107/2021-01-07_taxonomy"

CAT add_names -i out.CAT.contig2classification.txt -o "${OUTPUT_DIR}/02CAT/CAT_contigs_classification.txt" \
    -t "${DB_DIR}/CAT_prepare_20210107/2021-01-07_taxonomy"
CAT add_names -i out.CAT.ORF2LCA.txt -o "${OUTPUT_DIR}/02CAT/CAT_ORFs_classification.txt" \
    -t "${DB_DIR}/CAT_prepare_20210107/2021-01-07_taxonomy"

mv out.CAT.* "${OUTPUT_DIR}/02CAT/" || true
rm -rf out.CAT.alignment.diamond || true

# vContact2
echo "Running vContact2 analysis..."
activate_conda vcontact2
seqtk seq -L 10000 01Clustergenome/AD_VCs.fasta > "${OUTPUT_DIR}/04vcontact2/AD_VCs_10kb.fasta"
prodigal -i "${OUTPUT_DIR}/04vcontact2/AD_VCs_10kb.fasta" -o "${OUTPUT_DIR}/04vcontact2/AD_VCs_10kb.genes" \
    -a "${OUTPUT_DIR}/04vcontact2/AD_VCs_10kb.faa" -p meta

vcontact2_gene2genome -p "${OUTPUT_DIR}/04vcontact2/AD_VCs_10kb.faa" \
    -o "${OUTPUT_DIR}/04vcontact2/AD_VCs_10kb_g2g.csv" -s 'Prodigal-FAA'

vcontact2 -r "${OUTPUT_DIR}/04vcontact2/AD_VCs_10kb.faa" --rel-mode 'Diamond' \
    --proteins-fp "${OUTPUT_DIR}/04vcontact2/AD_VCs_10kb_g2g.csv" --db 'ProkaryoticViralRefSeq201-Merged' \
    --pcs-mode MCL --vcs-mode ClusterONE --c1-bin /home/kunyu/miniconda3/envs/vcontact2/bin/cluster_one-1.0.jar \
    --output-dir "${OUTPUT_DIR}/04vcontact2/vContactOut/"

# 12. RefSeq and IMG/VR Analysis
echo "Running RefSeq and IMG/VR analyses..."
activate_conda tools
check_file "${DB_DIR}/viral_refseq/viral.1.protein.faa"
diamond makedb --in "${DB_DIR}/viral_refseq/viral.1.protein.faa" -d "${OUTPUT_DIR}/08majority_rules/viral_refseq" --threads 10

diamond blastp -q 01Clustergenome/AD_VCs.faa -d "${OUTPUT_DIR}/08majority_rules/viral_refseq" \
    -o "${OUTPUT_DIR}/08majority_rules/AD_VCs_blastp.txt" --query-cover 50 --subject-cover 50 --evalue 1e-5 -k 1 --threads 10

csvtk filter -t -f "12>=50" "${OUTPUT_DIR}/08majority_rules/AD_VCs_blastp.txt" > \
    "${OUTPUT_DIR}/08majority_rules/AD_VCs_blastp_50score.txt"

cut -f 2 "${OUTPUT_DIR}/08majority_rules/AD_VCs_blastp_50score.txt" > \
    "${OUTPUT_DIR}/08majority_rules/AD_VCs_blastp_accession.txt"

ensure_dir "${OUTPUT_DIR}/08majority_rules/02blastp_accession"
split -l 1000 "${OUTPUT_DIR}/08majority_rules/AD_VCs_blastp_accession.txt" \
    "${OUTPUT_DIR}/08majority_rules/02blastp_accession/part_"

for file in "${OUTPUT_DIR}/08majority_rules/02blastp_accession/part_"*; do
    rg -f "$file" "${DB_DIR}/protein_taxID/prot.accession2taxid.FULL" --no-line-number >> \
        "${OUTPUT_DIR}/08majority_rules/AD_VCs_refseq_accession2taxid.txt"
done

cut -f 2 "${OUTPUT_DIR}/08majority_rules/AD_VCs_refseq_accession2taxid.txt" > \
    "${OUTPUT_DIR}/08majority_rules/AD_VCs_blastp_taxid.txt"

taxonkit lineage "${OUTPUT_DIR}/08majority_rules/AD_VCs_blastp_taxid.txt" --data-dir "${DB_DIR}/taxdump" | \
    taxonkit reformat -F --data-dir "${DB_DIR}/taxdump" | \
    cut -f 1,3 > "${OUTPUT_DIR}/08majority_rules/AD_VCs_taxid_taxonomy.txt"

paste "${OUTPUT_DIR}/08majority_rules/AD_VCs_refseq_accession2taxid.txt" \
    "${OUTPUT_DIR}/08majority_rules/AD_VCs_taxid_taxonomy.txt" | \
    cut -f 1,4 > "${OUTPUT_DIR}/08majority_rules/AD_VCs_accession_taxonomy.txt"

awk -F'\t' 'FNR==NR{a[$1]=$2; next}; {if($2 in a) print $0 "\t" a[$2]; else print $0 "\tNA"}' \
    "${OUTPUT_DIR}/08majority_rules/AD_VCs_accession_taxonomy.txt" \
    "${OUTPUT_DIR}/08majority_rules/AD_VCs_blastp_50score.txt" > \
    "${OUTPUT_DIR}/08majority_rules/AD_VCs_blastp_viref_50_tax.txt"

awk -F'\t' -v OFS='\t' '{sub(/_[0-9]+$/, "", $1); split($13, a, ";"); $13=a[1] ";" a[2] ";" a[3] ";" a[4] ";" a[5]; print $1, $13, 1}' \
    "${OUTPUT_DIR}/08majority_rules/AD_VCs_blastp_viref_50_tax.txt" > \
    "${OUTPUT_DIR}/08majority_rules/AD_VCs_family.txt"

awk -F'\t' '{sums[$1"\t"$2] += $3} END {for (key in sums) print key"\t"sums[key]}' \
    "${OUTPUT_DIR}/08majority_rules/AD_VCs_family.txt" > \
    "${OUTPUT_DIR}/08majority_rules/AD_VCs_family_stax.txt"

cut -f 1,3 "${OUTPUT_DIR}/08majority_rules/AD_VCs_family_stax.txt" | \
    datamash -sW -g1 sum 2 > "${OUTPUT_DIR}/08majority_rules/AD_VCs_contigs_stax.txt"

awk 'FNR==NR{a[$1]=$2; next}; {if($1 in a) print $0 "\t" a[$1]; else print $0 "\tNA"}' \
    "${OUTPUT_DIR}/08majority_rules/AD_VCs_contigs_stax.txt" \
    "${OUTPUT_DIR}/08majority_rules/AD_VCs_family_stax.txt" > \
    "${OUTPUT_DIR}/08majority_rules/AD_VCs_family_contigs.txt"

# IMG/VR
check_file "${DB_DIR}/viral_refseq/viral.1.1.genomic.fna"
makeblastdb -dbtype nucl -in "${DB_DIR}/viral_refseq/viral.1.1.genomic.fna" -input_type fasta \
    -title virus_genome -out "${OUTPUT_DIR}/09blast2IMGVR/virus_genome" -blastdb_version 5

blastn -query 01Clustergenome/AD_VCs.fasta -out "${OUTPUT_DIR}/09blast2IMGVR/AD_VCs_blastnVG.txt" \
    -outfmt "6 qseqid sseqid length qlen slen qstart qend sstart send evalue pident staxids sscinames scomnames sblastnames bitscore salltitles qcovs qcovhsp stitle" \
    -db "${OUTPUT_DIR}/09blast2IMGVR/virus_genome" -dust no -max_target_seqs 1 -perc_identity 95 -evalue 0.00001 -num_threads 20

blastn -query 01Clustergenome/AD_VCs.fasta -out "${OUTPUT_DIR}/09blast2IMGVR/AD_VCs_blastnIMGVR.txt" \
    -outfmt "6 qseqid sseqid length qlen slen qstart qend sstart send evalue pident staxids sscinames scomnames sblastnames bitscore salltitles qcovs qcovhsp stitle" \
    -db "${DB_DIR}/IMG_VR/IMGVR" -dust no -max_target_seqs 1 -perc_identity 95 -evalue 0.00001 -num_threads 20

cut -f 2 "${OUTPUT_DIR}/09blast2IMGVR/AD_VCs_blastnIMGVR.txt" > "${OUTPUT_DIR}/09blast2IMGVR/AD_IMGVR_VCblastn.txt"
csvtk grep -t -P "${OUTPUT_DIR}/09blast2IMGVR/AD_IMGVR_VCblastn.txt" "${DB_DIR}/IMG_VR/IMGVR_all_Host_information-high_confidence.tsv" \
    > "${OUTPUT_DIR}/09blast2IMGVR/IMGVR_host_results.txt"

# 13. Novelty Analysis
echo "Running novelty analysis..."
diamond blastp -q 01Clustergenome/AD_VCs.faa -d "${DB_DIR}/IMG_VR/IMGVR_proteins" \
    -o "${OUTPUT_DIR}/10novelty/01AD_VCs_IMGVR.txt" --id 30 --query-cover 50 --evalue 1e-5 -k 1 --threads 10
diamond blastp -q 01Clustergenome/AD_VCs.faa -d "${OUTPUT_DIR}/08majority_rules/viral_refseq" \
    -o "${OUTPUT_DIR}/10novelty/02AD_VCs_Refseq.txt" --id 30 --query-cover 50 --evalue 1e-5 -k 1 --threads 10

echo "Virus module completed successfully!"