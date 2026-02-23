#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=8G
#SBATCH --job-name=pySCENIC
#SBATCH --output=pySCENIC.%j.out

# inizialize conda
source /home/ciervo/miniconda3/etc/profile.d/conda.sh
conda activate pyscenic

# setwd
cd /home3/ciervo/scMULTIOME/Analisi/scRNA/Pozniak_SCENIC/
  export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

arboreto_with_multiprocessing.py \
pySCENIC_input_malignant.loom \
/home/ciervo/HeadAndNeck/Analisi_03aprile24/Merged/SCENIC/cisTarget_databases/allTFs_hg38.txt \
--method grnboost2 \
--output expr_mat.adjacencies.tsv \
--num_workers 50 \
--seed 1234

pyscenic ctx \
expr_mat.adjacencies.tsv \
/home/ciervo/HeadAndNeck/Analisi_03aprile24/Merged/SCENIC/cisTarget_databases/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather \
/home/ciervo/HeadAndNeck/Analisi_03aprile24/Merged/SCENIC/cisTarget_databases/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather \
--annotations_fname /home/ciervo/HeadAndNeck/Analisi_03aprile24/Merged/SCENIC/cisTarget_databases/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname pySCENIC_input_malignant.loom \
--mode "custom_multiprocessing" \
--output regulons.csv \
--num_workers 50

pyscenic aucell \
pySCENIC_input_malignant.loom \
regulons.csv \
-o auc_mtx.csv \
--num_workers 50