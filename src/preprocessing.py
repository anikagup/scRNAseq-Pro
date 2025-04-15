import scanpy as sc
import pandas as pd
import numpy as np
import os
import json
import matplotlib.pyplot as plt
import scrublet as scr


# Load data based on file type
def load_data(file, file_type):
    if file_type == "10x" or (file_type == "auto" and file.endswith(".h5")):
        adata = sc.read_10x_h5(file)
    elif file_type == "loom" or (file_type == "auto" and file.endswith(".loom")):
        adata = sc.read_loom(file)
    elif file_type == "h5ad" or (file_type == "auto" and file.endswith(".h5ad")):
        adata = sc.read_h5ad(file)
    elif file_type in ["csv", "txt"] or (file_type == "auto" and file.endswith(".csv")):
        df = pd.read_csv(file, index_col=0)
        adata = sc.AnnData(df.T)  # Transpose to match scRNA-seq format
    else:
        raise ValueError("Unsupported file type.")
    return adata

# Run Scrublet
def detect_doublets(adata):
    # Ensure raw count matrix is used
    counts_matrix = adata.raw.X if adata.raw is not None else adata.X
    if hasattr(counts_matrix, "toarray"):
        counts_matrix = counts_matrix.toarray()

    scrub = scr.Scrublet(counts_matrix)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()

    adata.obs["doublet_score"] = doublet_scores
    adata.obs["predicted_doublet"] = predicted_doublets

    # Optional: save histogram
    import matplotlib.pyplot as plt
    scrub.plot_histogram()
    plt.savefig("figures/doublet_score_histogram.png")

    # Filter out doublets
    adata = adata[~predicted_doublets].copy()
    print(f"✅ Removed {np.sum(predicted_doublets)} predicted doublets.")
    return adata
    
# Function to calculate QC metrics and generate plots
def generate_qc_metrics(adata):
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # Identify mitochondrial genes
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    # Generate violin plots for QC metrics
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True, save="_qc_metrics.png", show = False)
    ## Maybe make this three separate figs 
    print("QC metrics calculated and violin plots saved as '_qc_metrics.png'")

# Preprocessing function
def preprocess_data(adata, params):
    generate_qc_metrics(adata)  # Calculate QC metrics before preprocessing
    
    sc.pp.filter_cells(adata, min_genes=params.get("min_genes", 200)) #Changeable
    sc.pp.filter_genes(adata, min_cells=params.get("min_cells", 3)) #Changeable

    print("Filtering is done")

    # Identify mitochondrial genes
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")

    # Calculate QC metrics including mitochondrial percentage
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    # Filter out cells with high mitochondrial content (usually > 5–10%)
    adata = adata[adata.obs["pct_counts_mt"] < 10].copy()

    print("Successfully filtered out cells w/>10% mitochondrial content")

    
    adata=detect_doublets(adata)
    print("Successfully filtered out doublets")

    sc.pp.normalize_total(adata, target_sum=params.get("target_sum", 1e4))
    sc.pp.log1p(adata)
    print("Normalization is done")
    sc.pp.highly_variable_genes(adata, n_top_genes=params.get("n_top_genes", 2000)) #Changeable
    print("Highly variable genes have been found")
    #make a copy of data before subsetting highly variable genes to feed to ML model for labeling
    ML_data = adata.copy()
    adata = adata[:, adata.var["highly_variable"]]

    sc.pp.scale(adata, max_value=params.get("max_value", 10)) #Download CSV of count matrix
    sc.pp.scale(ML_data, max_value=params.get("max_value", 10)) #Download CSV of count matrix

    return adata, ML_data


def export_adata_to_csv(adata, output_path="processed_data/processed_matrix.csv"):
    dense_matrix = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X
    df = pd.DataFrame(dense_matrix, index=adata.obs_names, columns=adata.var_names)
    df.to_csv(output_path)
    print(f"✅ Processed expression matrix exported to {output_path}")
