import os
import scanpy as sc

# Get path of current script
script_dir = os.path.dirname(os.path.abspath(__file__))
upload_path = os.path.join(script_dir, "uploads", "latest.h5ad")

# # Ensure uploads folder exists
# os.makedirs(os.path.dirname(upload_path), exist_ok=True)

file = sc.read_h5ad(upload_path)
print(file.obs)
