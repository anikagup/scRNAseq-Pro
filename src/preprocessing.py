import scanpy as sc
import anndata as ad
import pandas as pd


#This is an example for building the program 
x="/Users/anika/Downloads/umis.mtx"

def preprocess(file_path):
    print(f"Processing file: {file_path}")
    # Incorporate reading of all file types using these commands: https://scanpy.readthedocs.io/en/latest/api/reading.html

    adata=sc.read_mtx(file_path)
    adata.var_names_make_unique()
    return adata

pre=preprocess(x)
print(pre)
