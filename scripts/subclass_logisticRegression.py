#!/usr/bin/env python3
# coding: utf-8

# JM 08/08/24
# In this notebook, we're going to re-apply machine learning to analyze the mouse retinal data and determine proper marker genes.

import datetime
print(f'{datetime.datetime.now()} Analysis Setup')

import sklearn as sk
import anndata as ad
import scanpy as sc 
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import joblib
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
import seaborn as sns
import os

os.chdir('/project/hipaa_ycheng11lab/atlas/CAMR2024')
sc.settings.n_jobs = -1

number_of_features = 20
ncols = 2

adata_full = ad.read_h5ad('data/camr_scrublet_batch_filtered.h5ad') # 

print(f'{datetime.datetime.now()} Analysis Time')

for gene_filter in ['highly_variable', 'moderately_variable', 'complete']:

    print(f'{datetime.datetime.now()} Filter Strategy: {gene_filter}')

    if gene_filter == 'highly_variable':
        adata = adata_full[:, adata_full.var.highly_variable]
    if gene_filter == 'moderately_variable':
        adata = adata_full[:, adata_full.var['dispersions'] > min(adata_full.var.loc[adata_full.var['highly_variable'], 'dispersions'])]
    if gene_filter == 'complete':
        adata = adata_full

    for majorclass in adata.obs['majorclass'].cat.categories:
    
        print(f'{datetime.datetime.now()} Major Class: {majorclass}')

        y = adata.obs.loc[adata.obs['majorclass'] == majorclass, 'author_cell_type']
        
        if len(y.unique()) < 2:
            continue
        
        le = LabelEncoder()
        _ = le.fit_transform(y)
        
        model_filename = f'models/ovr_logisticRegression_{gene_filter}_{majorclass}_author_cell_type.pkl'
        ovr_classifier_subclass = joblib.load(model_filename)

        print(le.inverse_transform(ovr_classifier_subclass.classes_))
        
        # Get class-specific coefficients
        class_coefficients = ovr_classifier_subclass.coef_
        num_classes = len(le.classes_)
        nrows = (num_classes // ncols) + (num_classes % ncols)
        
        # Get the top 20 important features for each class for each direction
        all_top_features_df = pd.DataFrame(columns=['Cell', 'Subclass', 'Ensembl', 'Gene', 'Coefficient'])
        for idx, class_name in enumerate(le.inverse_transform(ovr_classifier_subclass.classes_)):
            
            print(f'{datetime.datetime.now()} Subclass: {class_name}')
            
            if len(ovr_classifier_subclass.classes_) == 2 and idx == 1:
                break
            
            coefficients = class_coefficients[idx]
            feature_importance = pd.DataFrame({
                'Cell': majorclass,
                'Subclass': class_name,
                'Ensembl': adata.var_names,
                'Gene': adata.var['feature_name'].astype(str),  # Ensure feature_name is string
                'Coefficient': coefficients
            })
            top_features_df = feature_importance.sort_values(by='Coefficient', ascending=False).head(number_of_features)
            bottom_features_df = feature_importance.sort_values(by='Coefficient', ascending=True).head(number_of_features)
            all_top_features_df = pd.concat([all_top_features_df, top_features_df, bottom_features_df], ignore_index=True)
        # End subclass
        all_top_features_df.to_csv(f'spreadsheets/ovr_top_20_{gene_filter}_genes_{majorclass}_author_cell_type.csv', index=False)
    # End majorclass
# End of gene_filter
