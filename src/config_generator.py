import json
import os

# Default config template
default_config = {
    "input_file": "data/sample_data.h5",
    "file_type": "auto",
    "preprocessing_params": {
        "min_genes": 200,
        "min_cells": 3,
        "target_sum": 1e4,
        "n_top_genes": 2000,
        "max_value": 10,
        "n_neighbors": 10,
        "n_pcs": 40
    }
}

# Function to get user input with default values
def get_input(prompt, default):
    user_input = input(f"{prompt} (default: {default}): ")
    return type(default)(user_input) if user_input else default

# Check if config.json already exists
if os.path.exists("config.json"):
    with open("config.json", "r") as file:
        config = json.load(file)
else:
    config = default_config

# Get inputs for file settings
config["input_file"] = input("Enter the path to the input file (default: data/sample_data.h5): ") or config["input_file"]
config["file_type"] = input("Enter file type (10x, loom, h5ad, csv, txt, auto) (default: auto): ") or config["file_type"]

# Get inputs for preprocessing parameters
print("\n--- Preprocessing Parameters ---")
params = config["preprocessing_params"]
params["min_genes"] = get_input("Minimum genes per cell", params["min_genes"])
params["min_cells"] = get_input("Minimum cells per gene", params["min_cells"])
params["target_sum"] = get_input("Normalization target sum", params["target_sum"])
params["n_top_genes"] = get_input("Number of highly variable genes", params["n_top_genes"])
params["max_value"] = get_input("Maximum scale value", params["max_value"])
params["n_neighbors"] = get_input("Number of neighbors for UMAP", params["n_neighbors"])
params["n_pcs"] = get_input("Number of principal components for UMAP", params["n_pcs"])

# DEG inputs
print("\n--- Differential Expression Analysis Parameters ---")
de_config = config.get("differential_expression", {})
de_config["resolution"] = get_input("Clustering resolution", de_config.get("resolution", 0.5))
de_config["num_degs"] = get_input("Number of DEGs to display", de_config.get("num_degs", 20))
config["differential_expression"] = de_config

# Save updated config
with open("config.json", "w") as file:
    json.dump(config, file, indent=4)

print("\nConfiguration saved to config.json!")

