import json
import scanpy as sc
import os

from preprocessing import preprocess_data, load_data

import os
import json

# Get the absolute path of the directory containing main.py (which is src/)
SRC_DIR = os.path.dirname(os.path.abspath(__file__))

# Define the correct path to config.json (always inside src/)
CONFIG_PATH = os.path.join(SRC_DIR, "config.json")

print(f"📂 Looking for config.json at: {CONFIG_PATH}")  # Debugging print

# Load config.json safely
if not os.path.exists(CONFIG_PATH):
    raise FileNotFoundError(f"❌ ERROR: config.json not found at {CONFIG_PATH}")

with open(CONFIG_PATH, "r") as config_file:
    config = json.load(config_file)

print(f"✅ Loaded config from: {CONFIG_PATH}")


# Define file paths
input_file = config.get("input_file", "data/input_file")
file_type = config.get("file_type", "auto")


input_file = config.get("input_file", "data/pbmc3k.h5ad")
file_type = config.get("file_type", "h5ad")
params = config.get("preprocessing_params", {})

print("📂 Loading dataset...")
adata = load_data(input_file, file_type)

if adata is None:
    print("❌ ERROR: No valid data loaded. Exiting.")
    exit(1)

print(f"✅ Data loaded! Shape: {adata.shape}")  # Debugging step

# Preprocess data
print("🔄 Running preprocessing...")
adata = preprocess_data(adata, params)
print("✅ Preprocessing complete!")

# Perform clustering & UMAP
print("🔄 Running Leiden clustering...")
sc.tl.leiden(adata, resolution=config.get("differential_expression", {}).get("resolution", 0.5))
print("✅ Leiden clustering complete!")

print("🔄 Running PCA & Neighbors...")
sc.tl.pca(adata, svd_solver='arpack')  # Ensure PCA is computed
sc.pp.neighbors(adata, n_neighbors=params.get("n_neighbors", 10), n_pcs=params.get("n_pcs", 40))
print("✅ PCA & neighbors computed!")

print("🔄 Running UMAP...")
sc.tl.umap(adata)
print("✅ UMAP computation complete!")

# Save processed data
os.makedirs("data", exist_ok=True)
adata.write("data/processed_data.h5ad")
print("✅ Processed data saved to 'data/processed_data.h5ad'")

