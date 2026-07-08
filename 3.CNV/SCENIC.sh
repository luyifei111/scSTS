#!/bin/bash

#SBATCH -J scenic-ALL_test
#SBATCH -p protein
#SBATCH -n 40
#SBATCH -N 1
#SBATCH --ntasks-per-node=40
#SBATCH -o /home/u22211240018/STS/scenic/logall/%j.out
#SBATCH -e /home/u22211240018/STS/scenic/logall/%j.err                 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=22211240018@m.fudan.edu.cn

source activate
conda activate pyscenic310

f_loom_path_scenic='/home/u22211240018/STS/tumor_loom/STS_tumor_95samples_nomiro_harmony_nodoublet_20240507.loom'
F_TFS='/home/u22211240018/STS/seperated/T1620/scenic/allTFs_hg38.txt'
ADJACENCIES_FNAME='/home/u22211240018/STS/scenic/adj/all_adj.csv'

pyscenic grn f_loom_path_scenic F_TFS -o ADJACENCIES_FNAME --num_workers 40

MOTIFS_FNAME='/home/u22211240018/STS/scenic/reg/all_reg.csv'
DBS_PARAM='/home/u22211240018/STS/seperated/T1620/scenic/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather'
f_motif_path='/home/u22211240018/STS/seperated/T1620/scenic/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl'

pyscenic ctx ADJACENCIES_FNAME DBS_PARAM --annotations_fname f_motif_path --expression_mtx_fname f_loom_path_scenic --output MOTIFS_FNAME --mask_dropouts --num_workers 40

f_pyscenic_output='/home/u22211240018/STS/scenic/output_loom/all_pyscenic_output.loom'
pyscenic aucell f_loom_path_scenic MOTIFS_FNAME --output f_pyscenic_output --num_workers 40