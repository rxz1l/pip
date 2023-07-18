export LD_LIBRARY_PATH=/public/work/Personal/liujinxu/software/Anaconda/Anaconda3-2021.11/envs/report_python/lib/:$LD_LIBRARY_PATH
source /public/work/Personal/liujinxu/software/Anaconda/Anaconda3-2021.11/bin/activate /public/work/Personal/liujinxu/software/Anaconda/Anaconda3-2021.11/envs/report_python

sample=(S10-TC S1-CB GM-S3-A7A)

for i in ${sample[*]}
do
    {
    /public/work/Personal/liujinxu/software/R/R-4.1.2/bin/Rscript /public/work/Project/Single_cell/SA2023041801-03_peiyangji_renxiaozhen/c4_average1.R \
        --RDSor10X RDS \
        --cellranger_path /public/work/Project/Single_cell/SA2023041801-03_peiyangji_renxiaozhen/08.seurat_qc/$i/QC/${i}_QC.rds \
        --sample $i \
        --species Human \
        --outdir  /public/work/Project/Single_cell/SA2023041801-03_peiyangji_renxiaozhen/08.seurat_qc/average/$i
    }&

done