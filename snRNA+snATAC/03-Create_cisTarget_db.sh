#!bin/bash
REGION_BED="/home3/ciervo/scMULTIOME/Analisi/scRNA_scATAC/malignant_consensusPeaks.bed"
GENOME_FASTA="/storage/qnap_vol1/SHARED/NGSTOOLS/BCBIO2023/genomes/Hsapiens/hg38/seq/hg38.fa"
CHROMSIZES="/home3/ciervo/scMULTIOME/Analisi/scRNA_scATAC/hg38.chrom.sizes"
SCRIPT_DIR="/home3/ciervo/scMULTIOME/Analisi/scRNA_scATAC/SCENICplus/create_cisTarget_databases"

# ${SCRIPT_DIR}/create_fasta_with_padded_bg_from_bed.sh \
#         ${GENOME_FASTA} \
#         ${CHROMSIZES} \
#         ${REGION_BED} \
#         hg38.10x_NIBIT.with_1kb_bg_padding.fa \
#         1000 \
#         yes

OUT_DIR="/home3/ciervo/scMULTIOME/Analisi/scRNA_scATAC/SCENICplus/cisTarget_db_output/"
CBDIR="/home3/ciervo/scMULTIOME/Analisi/scRNA_scATAC/SCENICplus/aertslab_motif_colleciton/v10nr_clust_public/singletons"
FASTA_FILE="/home3/ciervo/scMULTIOME/Analisi/scRNA_scATAC/SCENICplus/hg38.10x_NIBIT.with_1kb_bg_padding.fa"
MOTIF_LIST="/home3/ciervo/scMULTIOME/Analisi/scRNA_scATAC/SCENICplus/motifs.txt"
DATABASE_PREFIX='NIBIT_M4'

"${SCRIPT_DIR}/create_cistarget_motif_databases.py" \
    -f ${FASTA_FILE} \
    -M ${CBDIR} \
    -m ${MOTIF_LIST} \
    -o ${OUT_DIR}/${DATABASE_PREFIX} \
    --bgpadding 1000 \
    -t 75