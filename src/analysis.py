import scanpy as sc
import json

# Load configuration
with open("config.json", "r") as config_file:
    config = json.load(config_file)

# Load preprocessed data
adata = sc.read_h5ad("data/processed_data.h5ad")

# Perform clustering
sc.tl.leiden(adata, resolution=config.get("differential_expression", {}).get("resolution", 0.5))
sc.pl.umap(adata, color="leiden")

# Differential Gene Expression Analysis
num_degs = config.get("differential_expression", {}).get("num_degs", 20)
sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")
sc.pl.rank_genes_groups(adata, n_genes=num_degs, sharey=False)

# Save DEGs to CSV
for cluster in adata.obs["leiden"].cat.categories:
    df = sc.get.rank_genes_groups_df(adata, group=cluster)
    df.to_csv(f"data/degs_cluster_{cluster}.csv", index=False)


print("Differential expression analysis completed and results saved.")
