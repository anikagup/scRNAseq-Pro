import scanpy as sc
import json
import os
import matplotlib.pyplot as plt
import pandas as pd 
import warnings
import shutil

warnings.filterwarnings("ignore", category=RuntimeWarning)

# Function to perform clustering & generate UMAP
def generate_umap(adata, config):

    print("🔄 Running PCA & Neighbors...")
    sc.tl.pca(adata, svd_solver='arpack')  # Ensure PCA is computed
    sc.pp.neighbors(adata, n_neighbors=config.get("preprocessing_params", {}).get("n_neighbors", 10),
                    n_pcs=config.get("preprocessing_params", {}).get("n_pcs", 40))
    
    print("✅ PCA & neighbors computed!")
    print("🔄 Running Leiden clustering...")
    
    sc.tl.leiden(
        adata,
        resolution=config.get("differential_expression", {}).get("resolution", 0.5),
        flavor="igraph",
        directed=False,
        n_iterations=2
    )

    print("✅ Leiden clustering complete!")

    print("🔄 Running UMAP...")
    sc.tl.umap(adata)
    print("✅ UMAP computation complete!")

    # Generate and save UMAP plots
    figures_dir = "figures"

    # Recreate the directory
    os.makedirs(figures_dir, exist_ok=True)

    # Default UMAP plots
    sc.pl.umap(adata, color=["leiden", "total_counts", "n_genes_by_counts"], wspace=0.4, save=f"_qc.png", show=False)
    sc.pl.umap(adata, color=adata.var_names[:5], wspace=0.4, save=f"_top5.png", show=False)

    print(f"✅ UMAP plots saved")

def user_umap(adata, params):
    requested_genes = params.get("custom_genes", ["IL-1"])
    print("received gene input")
    # Always treat as a list
    if isinstance(requested_genes, str):
        requested_genes = [requested_genes]

    # Strip and filter by genes that exist
    valid_genes = [gene.strip() for gene in requested_genes if gene.strip() in adata.var_names]
    if not valid_genes:
        print("❌ No valid genes found in the dataset for custom UMAP.")
        return
    # Plot valid genes
    sc.pl.umap(adata, color=valid_genes, save="_custom_gene.png", show=False)
    print(f"✅ User-defined UMAP generated for: {valid_genes}")

def prediction_umap(adata, config):
    sc.pl.umap(adata, color="Classifier_predictions", save=f"ML_umap.png", show=False)
    print(f"ML-classified plot is done!")

# Function to perform differential gene expression analysis
def perform_differential_expression(adata, config):
    num_degs = config.get("differential_expression", {}).get("num_degs", 20)
    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon", use_raw=False)
    sc.pl.rank_genes_groups(adata, n_genes=num_degs, sharey=False, save =f".png", show=False )

    # Save DEGs to CSV    
    all_degs = []
    top10_degs = []

    for cluster in adata.obs["leiden"].cat.categories:
        df = sc.get.rank_genes_groups_df(adata, group=cluster)
        df["cluster"] = cluster
        all_degs.append(df)

        # Keep only top 10 for summary CSV
        top10_degs.append(df.head(10))

    # Concatenate all DEGs across clusters
    all_degs_df = pd.concat(all_degs)
    top10_degs_df = pd.concat(top10_degs)

    # Save full and top-10 DEG tables
    all_degs_df.to_csv("processed_data/all_degs.csv", index=False)
    top10_degs_df.to_csv("processed_data/top10_degs_by_cluster.csv", index=False)


    print("✅ Differential expression analysis completed and results saved.")





