import scanpy as sc
import pandas as pd
import numpy as np
import os

# Define file paths
input_file = "data/input_file"
file_type = "auto"  # Options: '10x', 'loom', 'h5ad', 'csv', 'txt', 'fastq', or 'auto'

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
    elif file_type == "fastq":
        print("FASTQ files require alignment before Scanpy. Please preprocess first.")
        return None
    else:
        raise ValueError("Unsupported file type.")
    return adata

# Preprocessing function
def preprocess_data(adata):
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    adata = adata[:, adata.var["highly_variable"]]
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    return adata

# Main execution
adata = load_data(input_file, file_type)
if adata is not None:
    adata = preprocess_data(adata)
    sc.pl.umap(adata, color='CST3')  # Example gene for visualization
    adata.write("data/processed_data.h5ad")
else:
    print("No valid data loaded.")
