import scanpy as sc
import json
import os
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
import shutil

# Function to perform clustering & generate UMAP
def generate_umap(adata, config):
    print("🔄 Running Leiden clustering...")
    sc.tl.leiden(adata, resolution=config.get("differential_expression", {}).get("resolution", 0.5))
    print("✅ Leiden clustering complete!")

    print("🔄 Running PCA & Neighbors...")
    sc.tl.pca(adata, svd_solver='arpack')  # Ensure PCA is computed
    sc.pp.neighbors(adata, n_neighbors=config.get("preprocessing_params", {}).get("n_neighbors", 10),
                    n_pcs=config.get("preprocessing_params", {}).get("n_pcs", 40))
    print("✅ PCA & neighbors computed!")

    print("🔄 Running UMAP...")
    sc.tl.umap(adata)
    print("✅ UMAP computation complete!")

    # Generate and save UMAP plots
    figures_dir = "figures"

    # Check if the directory exists
    if os.path.exists(figures_dir):
        # Delete the directory and all its contents
        shutil.rmtree(figures_dir)

    # Recreate the directory
    os.makedirs(figures_dir, exist_ok=True)

    # Default UMAP plots
    sc.pl.umap(adata, color=["leiden", "total_counts", "n_genes_by_counts"], wspace=0.4, save=f"umap_qc.png")

    color_genes = list(adata.var_names[:5])
    with rc_context({"figure.figsize": (5, 1)}):
        sc.pl.umap(adata, color=color_genes, s=50, frameon=False, ncols=4, vmax="p99")
    #sc.pl.umap(adata, list(adata.var_names[:5]), save=f"genes_umap.png")  # Add first 5 genes for visualization

    #for feature in umap_features:
        #sc.pl.umap(adata, color=feature, save=f"_umap_{feature}.png")

    print(f"✅ UMAP plots saved")

# Function to perform differential gene expression analysis
def perform_differential_expression(adata, config):
    num_degs = config.get("differential_expression", {}).get("num_degs", 20)
    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")
    sc.pl.rank_genes_groups(adata, n_genes=num_degs, sharey=False)

    # Save DEGs to CSV

    os.makedirs("data", exist_ok=True)   
    for cluster in adata.obs["leiden"].cat.categories:
        df = sc.get.rank_genes_groups_df(adata, group=cluster)
        df.to_csv(f"data/degs_cluster_{cluster}.csv", index=False)

    print("✅ Differential expression analysis completed and results saved.")





