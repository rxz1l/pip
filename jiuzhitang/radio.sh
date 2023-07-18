export LD_LIBRARY_PATH=/public/work/Personal/liujinxu/software/Anaconda/Anaconda3-2021.11/envs/report_python/lib/:$LD_LIBRARY_PATH
source /public/work/Personal/liujinxu/software/Anaconda/Anaconda3-2021.11/bin/activate /public/work/Personal/liujinxu/software/Anaconda/Anaconda3-2021.11/envs/report_python

sample=(S10-TC S1-CB GM-S3-A7A)

for i in ${sample[*]}
do
    {
    /public/work/Personal/liujinxu/software/R/R-4.1.2/bin/Rscript /public/work/Project/Single_cell/SA2023041801-03_peiyangji_renxiaozhen/c4.R \
        --RDSor10X RDS \
        --cellranger_path /public/work/Project/Single_cell/SA2023041801-03_peiyangji_renxiaozhen/08.seurat_qc/$i/QC/${i}_QC.rds \
        --sample $i \
        --species Human \
        --mark /public/work/Pipline/Single_RNA/standard_analysis/C4/tongji_ratio/data/mark_type/mark_type_jiuzhitang.xls \
        --outdir  /public/work/Project/Single_cell/SA2023041801-03_peiyangji_renxiaozhen/08.seurat_qc/radio/$i
    }&
done