#!/bin/bash
#SBATCH -o job.%j.out
#SBATCH -e job.%j.err
#SBATCH -p node01
#SBATCH -c 20
#SBATCH --mem 150G

export LD_LIBRARY_PATH=/public/work/Personal/renxiaozhen/software/miniconda3/envs/scDRS/lib/:$LD_LIBRARY_PATH
source /public/work/Personal/renxiaozhen/software/miniconda3/bin/activate /public/work/Personal/renxiaozhen/software/miniconda3/envs/scDRS

# 参数 -----
##文件路径
h5ad_path=/public/work/Project/Single_cell/kangxuhe_mmy/07.GWAS/scDRS/data/test.h5ad
celltype_colname=general_pred_celltype
species=mouse
cov_path=
gs_path=/public/work/Project/Single_cell/kangxuhe_mmy/07.GWAS/scDRS/data/T2D.gs
outdir=/public/work/Project/Single_cell/kangxuhe_mmy/07.GWAS/scDRS/result

# run scDRS.py
python /public/work/Personal/renxiaozhen/pipline/scDRS/scDRS.py \
    --h5ad_path $h5ad_path \
    --celltype_colname $celltype_colname \
    --species $species \
    --cov_path $cov_path \
    --gs_path $gs_path \
    --outdir $outdir

