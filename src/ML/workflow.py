import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import joblib
from joblib import load
import sklearn 
from sklearn.ensemble import RandomForestClassifier


def predict_labels(adata):
    #subset highly variable genes used by model to predict labels
    HVG_list = pd.read_csv("src/ML/HVG_model_gene_namesV2.csv")
    HVG_list = HVG_list.iloc[0].tolist()  # Convert the first row of the DataFrame to a list
    
    intersection = adata.var_names.isin(HVG_list)

    data_subet = adata[: , intersection]

    #prediction numerical label to words
    label_dict = {
    0: 'Non Epithelial',
    1: 'Non Cancer Related',
    2: 'Early (Stage I/II)',
    3: 'Advanced (Stage III/IV)'}
    
    #load in model
    classifier = load("src/ML/classifierV2.joblib")

    #apply model to data
    predictions = classifier.predict(data_subet.X)

    #convert numerical labels to word predictions
    word_labels = [label_dict[num_label] for num_label in predictions]

    #create a new column in anndata object for model predictions
    adata.obs["Classifier_predictions"] = word_labels

    return adata
