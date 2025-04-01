import scanpy as sc
import json
import os
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
import shutil

def compare_dge_between_clusters(adata, n_genes=20):
    """Compare differential gene expression between UMAP clusters using Wilcoxon test."""
    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")
    sc.pl.rank_genes_groups(adata, n_genes=n_genes, sharey=False)
    # Save results to CSV
    for cluster in adata.obs["leiden"].cat.categories:
        df = sc.get.rank_genes_groups_df(adata, group=cluster)
        df.to_csv(f"data/degs_cluster_{cluster}.csv", index=False)
    print("✅ Differential expression analysis completed and results saved.")
    return adata
