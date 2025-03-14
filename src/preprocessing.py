import scanpy as sc
import pandas as pd
import numpy as np
import os
import json
import matplotlib.pyplot as plt

# Load configuration from config.json
with open("config.json", "r") as config_file:
    config = json.load(config_file)

# Define file paths
input_file = config.get("input_file", "data/input_file")
file_type = config.get("file_type", "auto")

# Load data based on file type
def load_data(file, file_type):
    if file_type == "10x" or (file_type == "auto" and file.endswith(".h5")):
        adata = sc.read_10x_h5(file)
    elif file_type == "loom" or (file_type == "auto" and file.endswith(".loom")):
        adata = sc.read_loom(file)
    elif file_type == "h5ad" or (file_type == "auto" and file.endswith(".h5ad")):
        adata = sc.read_h5ad(file)
    elif file_type in ["csv", "txt"] or (file_type == "auto" and (file.endswith(".csv") or file.endswith(".txt"))):
        adata = sc.read_text(file)
    else:
        raise ValueError("Unsupported file type.")
    return adata

# Function to calculate QC metrics and generate plots
def generate_qc_metrics(adata):
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # Identify mitochondrial genes
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    # Generate violin plots for QC metrics
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True, save="_qc_metrics.png")
    
    print("QC metrics calculated and violin plots saved as '_qc_metrics.png'")

# Preprocessing function
def preprocess_data(adata, params):
    generate_qc_metrics(adata)  # Calculate QC metrics before preprocessing
    
    sc.pp.filter_cells(adata, min_genes=params.get("min_genes", 200))
    sc.pp.filter_genes(adata, min_cells=params.get("min_cells", 3))
    sc.pp.normalize_total(adata, target_sum=params.get("target_sum", 1e4))
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=params.get("n_top_genes", 2000))
    adata = adata[:, adata.var["highly_variable"]]
    sc.pp.scale(adata, max_value=params.get("max_value", 10))
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=params.get("n_neighbors", 10), n_pcs=params.get("n_pcs", 40))
    sc.tl.umap(adata)
    return adata

# Main execution
params = config.get("preprocessing_params", {})
adata = load_data(input_file, file_type)
if adata is not None:
    adata = preprocess_data(adata, params)
    sc.pl.umap(adata, color='CST3')
    adata.write("data/processed_data.h5ad")
else:
    print("No valid data loaded.")