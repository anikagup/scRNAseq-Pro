import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import joblib
from joblib import load
import sklearn 
from sklearn.ensemble import RandomForestClassifier


def predict_labels(x):
    #subset highly variable genes used by model to predict labels
    HVG_list = pd.read_csv("/Users/eriknoyman/Desktop/USC/Senior/Semester 2/BME 405/scRNA Seq Automation/scRNA-seq-Automation/src/ML/HVG_model_gene_names.csv")
    HVG_list = HVG_list.iloc[0].tolist()  # Convert the first row of the DataFrame to a list
    data_subet = x[: , HVG_list]

    #prediction numerical label to words
    label_dict = {
    0: 'non epithelial',
    1: 'non cancer related',
    2: 'early',
    3: 'advanced'}
    
    #load in model
    classifier = load("/Users/eriknoyman/Desktop/USC/Senior/Semester 2/BME 405/scRNA Seq Automation/scRNA-seq-Automation/src/ML/classifierV1.joblib")

    #apply model to data
    predictions = classifier.predict(data_subet.X)

    #convert numerical labels to word predictions
    word_labels = [label_dict[num_label] for num_label in predictions]

    #create a new column in anndata object for model predictions
    x.obs["Classifier_predictions"] = word_labels



