#!/usr/bin/env python3
# coding: utf-8

# JM 08/08/24
# In this notebook, we're going to re-apply machine learning to analyze the mouse retinal data and determine proper marker genes.
# TODO: Merge the original script here and move the original to archive so that all 3 analyses can be done in one script.
# Perhaps even leave a flag to make it so that only 1 of the three is run at the users choosing.
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

adata = ad.read_h5ad('01_QualityControl/1_camr_scrublet_batch_filtered.h5ad')

for majorclass in adata.obs['majorclass'].cat.categories:

    print(f'{datetime.datetime.now()} Major Class: {majorclass}')
    
    # Extract feature matrix (X) and target vector (y)
    X = adata.X[adata.obs['majorclass'] == majorclass, :]
    y = adata.obs.loc[adata.obs['majorclass'] == majorclass, 'author_cell_type']
    
    if len(y.unique()) < 2:
        continue
    
    le = LabelEncoder()
    _ = le.fit_transform(y)
    
    # Split data into training and testing sets
    X_train, X_test, y_train_subclass, y_test_subclass = train_test_split(X, y_encoded, test_size=0.2, random_state=42)
    
    ovr_classifier_subclass = LogisticRegression(multi_class='ovr', max_iter=1000, random_state=42, n_jobs = 16)
    ovr_classifier_subclass.fit(X_train, y_train_subclass)
    
    model_filename = f'02_Modeling/minorclass/2_ovr_LogReg_minorclass-{majorclass}.pkl'
    joblib.dump(ovr_classifier_subclass, model_filename)
    
    # Validation
    y_pred_subclass_ovr = ovr_classifier_subclass.predict(X_test)
    target_names = le.inverse_transform(np.unique(y_test_subclass))
    validation_report = classification_report(y_test_subclass, y_pred_subclass_ovr, target_names=target_names)
    print(validation_report)
    with open(f"02_Modeling/minorclass/2_minorclass-{majorclass}_validation_report.txt", "w") as report_file:
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
    plt.title(f'{majorclass}: {100 * cm.diagonal().sum() / cm.sum()}% Accuracy')
    plt.savefig(f"02_Modeling/figures/confusion_matrix_minorclass-{majorclass}.pdf", bbox_inches='tight')
    plt.show()
    
    print(f"{majorclass} Minorclass:", le.inverse_transform(ovr_classifier_subclass.classes_))
    
    # Get class-specific coefficients
    class_coefficients = ovr_classifier_subclass.coef_
    num_classes = len(le.classes_)
    nrows = (num_classes // ncols) + (num_classes % ncols)
    
    # Get the top 20 important features for each class for each direction
    all_top_features_df = pd.DataFrame(columns=['Name', 'Major_Name', 'Ensembl', 'Marker', 'Coefficient'])
    for idx, class_name in enumerate(le.inverse_transform(ovr_classifier_subclass.classes_)):
        
        print(f'{datetime.datetime.now()} Subclass: {class_name}')
        
        if len(ovr_classifier_subclass.classes_) == 2 and idx == 1:
            continue # Binary classifications use only 1 regression equation
        
        coefficients = class_coefficients[idx]
        all_feature_importance = pd.DataFrame({
            'Name': class_name,
            'Major_Name': majorclass,
            'Ensembl': adata.var_names,
            'Marker': adata.var['feature_name'].astype(str),
            'Coefficient': coefficients
        })
        top_features_df = all_feature_importance.sort_values(by='Coefficient', ascending=False).head(number_of_features)
        bottom_features_df = all_feature_importance.sort_values(by='Coefficient', ascending=True).head(number_of_features)
        all_top_features_df = pd.concat([all_top_features_df, top_features_df, bottom_features_df], ignore_index=True)
    # End subclass
    all_top_features_df.to_csv(f'02_Modeling/minorclass/2_ovr_LogReg_minorclass-{majorclass}_AbsTop{number_of_features}Markers.txt', index=False, sep ='\t')
# End majorclass
