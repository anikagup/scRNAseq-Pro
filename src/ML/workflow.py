import pandas as pd
import numpy as np
from anndata import AnnData
from scipy.sparse import csr_matrix, hstack
from joblib import load

def align_anndata_to_genes(adata: AnnData, expected_genes: list) -> AnnData:
    """Ensure adata includes all genes in expected_genes, in correct order, with missing ones filled with zeros."""
    current_genes = adata.var_names
    present_genes = current_genes.intersection(expected_genes)
    missing_genes = [g for g in expected_genes if g not in current_genes]

    # Build sparse matrix for missing genes
    zero_matrix = csr_matrix((adata.n_obs, len(missing_genes)))

    # Combine existing and missing
    X_existing = adata[:, present_genes].X
    X_combined = hstack([X_existing, zero_matrix])

    # Reorder columns to match expected gene order
    all_genes = list(present_genes) + missing_genes
    gene_order_idx = [all_genes.index(g) for g in expected_genes]
    X_final = X_combined.tocsr()[:, gene_order_idx]

    # Return aligned AnnData
    aligned = AnnData(X=X_final, obs=adata.obs.copy())
    aligned.var_names = expected_genes
    return aligned

def predict_labels(adata, ML_data):
    # Load trained gene list (from feature importances CSV)
    HVG_df = pd.read_csv("src/ML/XGB_2kHVG_feature_importances.csv")
    
    # # Keep only genes with non-zero importance
    # HVG_list = HVG_df[HVG_df["importance"] > 0]["gene"].tolist()
    
    HVG_list = HVG_df["gene"].tolist()


    # Align ML_data to expected gene list
    ML_aligned = align_anndata_to_genes(ML_data, HVG_list)

    # Load XGBoost classifier
    classifier = load('src/ML/XGB_classifier_2kHVG.joblib')

    # Predict using aligned data
    X_input = ML_aligned.X.toarray() if hasattr(ML_aligned.X, "toarray") else ML_aligned.X
    predictions = classifier.predict(X_input)

    # Map numeric predictions to labels
    label_dict = {
        0: 'Non Epithelial',
        1: 'Non Cancer Related',
        2: 'Early (Stage I/II)',
        3: 'Advanced (Stage III/IV)'
    }
    word_labels = [label_dict[num_label] for num_label in predictions]

    # Add predictions to output AnnData object
    adata.obs["Classifier_predictions"] = word_labels

    return adata