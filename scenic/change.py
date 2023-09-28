import os, sys 
os.getcwd()
os.listdir(os.getcwd()) 

import argparse
import loompy as lp;
import numpy as np;
import scanpy as sc;

# 传入参数 ----
parser = argparse.ArgumentParser(description='')
parser.add_argument("--path", type = str)
args = parser.parse_args()

x=sc.read_csv(f"{args.path}/scenic.data.csv");
row_attrs = {"Gene": np.array(x.var_names),};
col_attrs = {"CellID": np.array(x.obs_names)};
lp.create(f"{args.path}/scenic.loom",x.X.transpose(),row_attrs,col_attrs);