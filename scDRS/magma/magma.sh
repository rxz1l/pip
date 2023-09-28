#! /bin/bash

export LD_LIBRARY_PATH=/public/work/Personal/renxiaozhen/software/magma/lib64
/public/work/Personal/renxiaozhen/software/magma/magma --annotate --snp-loc /public/work/Personal/renxiaozhen/software/magma/g1000/g1000_eur.bim \
     --gene-loc /public/work/Personal/renxiaozhen/software/magma/NCBI37/NCBI37.3.gene.loc --out g1000_eur
