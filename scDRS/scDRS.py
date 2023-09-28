# load packages ---- 
import argparse
import scdrs
import scanpy as sc
from anndata import AnnData
from scipy import stats
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import warnings


warnings.filterwarnings("ignore")
sc.set_figure_params(dpi=125)

# 传入参数 ----
parser = argparse.ArgumentParser(description='')
parser.add_argument("--h5ad_path", type = str)
parser.add_argument("--celltype_colname", type = str)
parser.add_argument("--species", type = str)
parser.add_argument("--cov_path", type = str)
parser.add_argument("--gs_path", type = str)
parser.add_argument("--outdir", type = str)

args = parser.parse_args()

# 创建文件夹 ----
os.makedirs(args.outdir, exist_ok = True)

# load data ----
## load adata
adata = sc.read_h5ad(args.h5ad_path)

## subset gene sets
df_gs = pd.read_csv(args.gs_path, sep  = "\t", index_col = 0)

# df_gs = df_gs.loc[
#     [
#         "PASS_Schizophrenia_Pardinas2018",
#         "spatial_dorsal",
#         "UKB_460K.body_HEIGHTz",
#     ],
#     :,
# ].rename(
#     {
#         "PASS_Schizophrenia_Pardinas2018": "SCZ",
#         "spatial_dorsal": "Dorsal",
#         "UKB_460K.body_HEIGHTz": "Height",
#     }
# )
# display(df_gs)

# df_gs.to_csv("/public/work/Personal/renxiaozhen/pipline/scDRS/data/processed_geneset.gs", sep="\t")

# scDRS analysis of disease enrichment for individual cells ----
os.system(f"""scdrs compute-score \\
    --h5ad-file {args.h5ad_path} \\
    --h5ad-species {args.species} \\
    --gs-file {args.gs_path} \\
    --gs-species {args.species} \\
    --cov-file {args.cov_path} \\
    --flag-filter-data True \\
    --flag-raw-count True \\
    --flag-return-ctrl-raw-score False \\
    --flag-return-ctrl-norm-score True \\
    --out-folder {args.outdir}
""")

# scDRS score umap plot ----
dict_score = {
    trait: pd.read_csv(f"{args.outdir}/{trait}.full_score.gz", sep="\t", index_col=0)
    for trait in df_gs.index
}
for trait in dict_score:
    adata.obs[trait] = dict_score[trait]["norm_score"]

sc.set_figure_params(figsize =[2.5, 2.5], dpi=300)
sc.pl.umap(
    adata,
    color=args.celltype_colname,
    ncols=1,
    color_map="RdBu_r",
    vmin=-5,
    vmax=5,
    show = False
)
plt.savefig(f"{args.outdir}/celltype.png")
plt.savefig(f"{args.outdir}/celltype.pdf")

sc.pl.umap(
    adata,
    color=dict_score.keys(),
    color_map="RdBu_r",
    vmin=-5,
    vmax=5,
    s=20,
    show=False
)
plt.savefig(f"{args.outdir}/scDRS.png")
plt.savefig(f"{args.outdir}/scDRS.pdf")

# scDRS test of group level statistics ----

for trait in dict_score.keys():
    os.system(f"""scdrs perform-downstream \\
        --h5ad-file {args.h5ad_path} \\
        --score-file {args.outdir}/{trait}.full_score.gz \\
        --out-folder {args.outdir} \\
        --group-analysis {args.celltype_colname} \\
        --flag-filter-data True \\
        --flag-raw-count True""")

dict_df_stats = {
    trait: pd.read_csv(f"{args.outdir}/{trait}.scdrs_group.{args.celltype_colname}", sep="\t", index_col=0)
    for trait in dict_score.keys()
}
# dict_celltype_display_name = {
#     "pyramidal_CA1": "Pyramidal CA1",
#     "oligodendrocytes": "Oligodendrocyte",
#     "pyramidal_SS": "Pyramidal SS",
#     "interneurons": "Interneuron",
#     "endothelial-mural": "Endothelial",
#     "astrocytes_ependymal": "Astrocyte",
#     "microglia": "Microglia",
# }

fig, ax = scdrs.util.plot_group_stats(
    dict_df_stats,
    plot_kws={
        "vmax": 0.2,
        "cb_fraction":0.12
    }
)
plt.savefig(f"{args.outdir}/cor_heatmap.png", bbox_inches = 'tight')
plt.savefig(f"{args.outdir}/cor_heatmap.pdf", bbox_inches = 'tight')



# # Further analysis of cell subsets ----
# adata_ca1 = adata[adata.obs["level2class"].isin(["CA1Pyr1", "CA1Pyr2"])].copy()
# sc.pp.filter_cells(adata_ca1, min_genes=0)
# sc.pp.filter_genes(adata_ca1, min_cells=1)
# sc.pp.normalize_total(adata_ca1, target_sum=1e4)
# sc.pp.log1p(adata_ca1)

# sc.pp.highly_variable_genes(adata_ca1, min_mean=0.0125, max_mean=3, min_disp=0.5)
# adata_ca1 = adata_ca1[:, adata_ca1.var.highly_variable]
# sc.pp.scale(adata_ca1, max_value=10)
# sc.tl.pca(adata_ca1, svd_solver="arpack")

# sc.pp.neighbors(adata_ca1, n_neighbors=10, n_pcs=40)
# sc.tl.umap(adata_ca1, n_components=2)

# # assign scDRS score
# for trait in dict_score:
#     adata_ca1.obs[trait] = dict_score[trait]["norm_score"]

# sc.pl.umap(
#     adata_ca1,
#     color=dict_score.keys(),
#     color_map="RdBu_r",
#     vmin=-5,
#     vmax=5,
#     s=20,
# )

# ## line chart
# df_plot = adata_ca1.obs[["Dorsal", "SCZ", "Height"]].copy()
# df_plot["Dorsal quintile"] = pd.qcut(df_plot["Dorsal"], 5, labels=np.arange(5))

# fig, ax = plt.subplots(figsize=(3.5, 3.5))
# for trait in ["SCZ", "Height"]:
#     sns.lineplot(
#         data=df_plot,
#         x="Dorsal quintile",
#         y=trait,
#         label=trait,
#         err_style="bars",
#         marker="o",
#         ax=ax,
#     )
# ax.set_xticks(np.arange(5))
# ax.set_xlabel("Dorsal quintile")
# ax.set_ylabel("Mean scDRS disease score")
