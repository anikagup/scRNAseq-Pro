import json
import scanpy as sc
import os
from preprocessing import preprocess_data, load_data

# Load configuration
with open("config.json", "r") as config_file:
    config = json.load(config_file)

input_file = config.get("input_file", "data/sample_data.h5")
file_type = config.get("file_type", "auto")
params = config.get("preprocessing_params", {})

# Load and preprocess data
adata = load_data(input_file, file_type)
if adata is not None:
    adata = preprocess_data(adata, params)

    # Run Leiden clustering for visualization
    sc.tl.leiden(adata, resolution=config.get("differential_expression", {}).get("resolution", 0.5))
    
    # Identify default metadata & gene-based UMAP coloring options
    umap_features = ["leiden", "total_counts", "pct_counts_mt", "n_genes_by_counts"]
    umap_features += list(adata.var_names[:5])  # Add first 5 genes

    # Update config.json with available UMAP colors
    config["visualization"]["umap_colors"] = umap_features
    with open("config.json", "w") as config_file:
        json.dump(config, config_file, indent=4)

    # Save processed dataset
    os.makedirs("data", exist_ok=True)
    adata.write("data/processed_data.h5ad")

    print("Preprocessing complete! Processed data saved to 'data/processed_data.h5ad'")
    print(f"Available UMAP color options: {umap_features}")
else:
    print("No valid data loaded.")
