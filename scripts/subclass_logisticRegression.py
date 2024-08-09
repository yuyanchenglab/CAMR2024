#!/usr/bin/env python3
# coding: utf-8

# JM 08/08/24
# In this notebook, we're going to re-apply machine learning to analyze the mouse retinal data and determine proper marker genes.

import sklearn as sk
import anndata as ad
import scanpy as sc 
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import joblib
import datetime
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

adata_full = ad.read_h5ad('data/camr_scrublet_batch_filtered.h5ad')

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
    
        # Extract feature matrix (X) and target vector (y)
        X = adata.X[adata.obs['majorclass'] == majorclass, :]
        y = adata.obs.loc[adata.obs['majorclass'] == majorclass, 'author_cell_type']

        if len(y.unique()) < 2:
            continue
        
        le = LabelEncoder()
        y_encoded = le.fit_transform(y)
        
        # Split data into training and testing sets
        X_train, X_test, y_train_subclass, y_test_subclass = train_test_split(X, y_encoded, test_size=0.2, random_state=42)
        
        ovr_classifier_subclass = LogisticRegression(multi_class='ovr', max_iter=1000, random_state=42, n_jobs = 16)
        ovr_classifier_subclass.fit(X_train, y_train_subclass)
        
        model_filename = f'models/ovr_logisticRegression_{gene_filter}_{majorclass}_author_cell_type.pkl'
        joblib.dump(ovr_classifier_subclass, model_filename)
        
        # Validation
        y_pred_subclass_ovr = ovr_classifier_subclass.predict(X_test)
        target_names = le.inverse_transform(np.unique(y_test_subclass))
        validation_report = classification_report(y_test_subclass, y_pred_subclass_ovr, target_names=target_names)
        print(validation_report)
        with open(f"spreadsheets/{gene_filter}_{majorclass}_validation_report.txt", "w") as report_file:
            report_file.write(validation_report)

        # Generate the confusion matrix
        cm = confusion_matrix(y_test_subclass, y_pred_subclass_ovr, labels=ovr_classifier_subclass.classes_) # le.classes_
        cm_normalized = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        plt.figure(figsize=(10, 9))
        sns.heatmap(cm_normalized, annot=True, fmt='.1%', cmap='Blues', 
                    xticklabels=le.inverse_transform(ovr_classifier_subclass.classes_), 
                    yticklabels=le.inverse_transform(ovr_classifier_subclass.classes_))
        plt.xlabel('Predicted Label')
        plt.ylabel('True Label')
        plt.title(f'{gene_filter} {majorclass}: {100 * cm.diagonal().sum() / cm.sum()}% Accuracy')
        plt.savefig(f"figures/modeling/confusion_matrix_{gene_filter}_{majorclass}_author_cell_type", bbox_inches='tight')
        plt.show()
        
        # Get class-specific coefficients
        class_coefficients = ovr_classifier_subclass.coef_
        num_classes = len(le.classes_)
        nrows = (num_classes // ncols) + (num_classes % ncols)
        
        # Create a figure with subplots
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(15, 5 * nrows))
        axes = axes.flatten()  # Flatten the 2D array of axes to 1D
        
        # Get the top 20 important features for each class
        all_top_features_df = pd.DataFrame(columns=['Cell Type', 'Gene', 'Coefficient'])
        for idx, class_name in enumerate(ovr_classifier_subclass.classes_):
            
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
            bottom_features_df = feature_importance.sort_values(by='Coefficient', ascending=False).head(number_of_features)
            all_top_features_df = pd.concat([all_top_features_df, top_features_df, bottom_features_df], ignore_index=True)
        
            # Plot the top 20 important features for each class
            sns.barplot(ax=axes[idx], x='Coefficient', y='Gene', data=top_features_df)
            axes[idx].set_title(f'Top {number_of_features} {gene_filter} Important Features for {class_name} of {majorclass}')
        # End subclass
        
        # Remove any unused subplots
        for ax in axes[num_classes:]:
            ax.remove()
        
        plt.tight_layout()
        plt.savefig(f'figures/modeling_interpretation/ovr_{gene_filter}_feature_importance_by_{majorclass}_author_cell_type')
    
        top_features_df.to_csv(f'spreadsheets/ovr_top_20_{gene_filter}_genes_{majorclass}_author_cell_type.csv', index=False)
    # End majorclass
# End of gene_filter
