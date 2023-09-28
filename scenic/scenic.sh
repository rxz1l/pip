#!/bin/bash
#SBATCH -o job.%j.out
#SBATCH -e job.%j.err
#SBATCH -p bnode
#SBATCH -c 20
#SBATCH --mem 150G
export LD_LIBRARY_PATH=/public/work/Personal/renxiaozhen/software/miniconda3/envs/pyscenic/lib/:$LD_LIBRARY_PATH
source /public/work/Personal/renxiaozhen/software/miniconda3/bin/activate /public/work/Personal/renxiaozhen/software/miniconda3/envs/pyscenic

# dir=/data/index_genome/cisTarget_databases/ #改成自己的目录
# tfs=$dir/hs_hgnc_tfs.txt
# feather=$dir/hg19-tss-centered-10kb-10species.mc9nr.genes_vs_motifs.rankings.feather
# tbl=$dir/motifs-v9-nr.hgnc-m0.001-o0.0.tbl 
# # 一定要保证上面的数据库文件完整无误哦 
# input_loom=scenic_data.loom
# ls $tfs  $feather  $tbl  

# 参数 -----
##rds文件路径
input_path=  
##随机抽取数据比例
percentage=
##过滤后数据输出路径
file_data_outdir=
##结果文件路径
result_dir=
##转录因子文件可以在https://link.zhihu.com/?target=https%3A//github.com/aertslab/pySCENIC/tree/master/resources中下载
feather=
##ctx数据库文件，根据自己的物种下载https://resources.aertslab.org/cistarget
tbl=
tfs=
##线程
num_workers=

# 过滤数据
/public/work/Personal/renxiaozhen/software/miniconda3/envs/R/bin/Rscript /public/work/Personal/renxiaozhen/pipline/scenic/data.r \
    --input_path $input_path --percentage $percentage --out_dir $file_data_outdir

# 运行change.py
/public/work/Personal/renxiaozhen/software/miniconda3/envs/pyscenic/bin/python /public/work/Personal/renxiaozhen/pipline/scenic/change.py \
    --path $file_data_outdir

# 运行pySCENIC
## 2.1 grn
pyscenic grn \
    --num_workers $num_workers \
    --output $result_dir/grn_result.csv \
    --method grnboost2 \
    $file_data_outdir/scenic.loom \
    $tfs #转录因子文件，1839个基因的名字列表

## 2.2 cistarget
pyscenic ctx \
    $result_dir/grn_result.csv\
    $feather \
    --annotations_fname $tbl \
    --expression_mtx_fname $file_data_outdir/scenic.loom  \
    --mode "dask_multiprocessing" \
    --output $result_dir/ctx.csv \
    --num_workers $num_workers  \
    --mask_dropouts

## 2.3 AUCell
pyscenic aucell \
    $file_data_outdir/scenic.loom \
    $result_dir/ctx.csv \
    --output $result_dir/out_SCENIC.loom \
    --num_workers  $num_workers