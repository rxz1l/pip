#!/usr/bin/bash
#SBATCH -o job.%j.out
#SBATCH -e job.%j.err
#SBATCH -p bnode,node01,node02,node03
#SBATCH -c 20
#SBATCH --mem 250G

unset PYTHONPATH 
export LD_LIBRARY_PATH=/public/work/Personal/liujinxu/software/Anaconda/Anaconda3-2021.11/envs/report_python/lib/:$LD_LIBRARY_PATH
source /public/work/Personal/liujinxu/software/Anaconda/Anaconda3-2021.11/bin/activate /public/work/Personal/liujinxu/software/Anaconda/Anaconda3-2021.11/envs/report_python

/public/work/Personal/liujinxu/software/R/R-4.1.2/bin/Rscript /public/work/Project/Single_cell/SA2022041401_jiuzhitang_10x_wkl/04.Trajectory/monocle_barV2.R \
                                                            --seurat_path /public/work/Project/Single_cell/atri_xbc_0921_rxz/03.sample_combind/harmony.rds \
                                                            --outdir  /public/work/Project/Single_cell/atri_xbc_0921_rxz/04.monocle/result \
                                                            --cluster_label /public/work/Project/Single_cell/atri_xbc_0921_rxz/04.monocle/monocle_cluster.csv \
                                                            --ordering_gene dispersion \
