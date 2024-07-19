#!/usr/bin/env python3
# coding: utf-8

# ## Modeling Notebook
# 
# 
# In this notebook, we're going to re-apply various machine learning approaches analyze the mouse retinal data and draw out all of the information as we can from the data.

# In[1]:


import sklearn as sk
import anndata as ad
import scanpy as sc 
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

sc.settings.n_jobs = -1


# In[2]:

import joblib

adata = ad.read_h5ad('camr_modeling_input.h5ad')
adata



from scipy.sparse import issparse
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import classification_report
import seaborn as sns


from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay

# Extract feature matrix (X) and target vector (y)
X = adata.X
y = adata.obs['author_cell_type']

le = LabelEncoder()
y_encoded = le.fit_transform(y)

# Split data into training and testing sets
X_train, X_test, y_train_subclass, y_test_subclass = train_test_split(X, y_encoded, test_size=0.2, random_state=42)

# Takes less than a minute to run
from sklearn.linear_model import LogisticRegression

from scipy.sparse import issparse
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import classification_report
import seaborn as sns


# In[ ]:


from sklearn.linear_model import LogisticRegression

ovr_classifier_subclass = LogisticRegression(multi_class='ovr', max_iter=1000, random_state=42, n_jobs = 16)
ovr_classifier_subclass.fit(X_train, y_train_subclass)

model_filename = 'models/ovr_classifier_subclass.pkl'
joblib.dump(ovr_classifier_subclass, model_filename)


# In[ ]:


y_pred_subclass_ovr = ovr_classifier_subclass.predict(X_test)
target_names = le.inverse_transform(np.unique(y_test_subclass))
print(classification_report(y_test_subclass, y_pred_subclass_ovr, target_names=target_names))


# In[ ]:


# Get class-specific coefficients
class_coefficients = ovr_classifier_subclass.coef_
num_classes = len(le.classes_)
number_of_features = 20

# Determine the grid size (e.g., 2 rows, num_classes / 2 columns)
nrows = (num_classes // 2) + (num_classes % 2)
ncols = 2

# Create a figure with subplots
fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(15, 5 * nrows))
axes = axes.flatten()  # Flatten the 2D array of axes to 1D

# Plot the top 5 important features for each class in a grid
for idx, class_name in enumerate(le.classes_):
    coefficients = class_coefficients[idx]
    feature_importance_ovr_df = pd.DataFrame({
        'ESUM': adata.var_names,
        'feature_name': adata.var['feature_name'].astype(str),  # Ensure feature_name is string
        'coefficient': coefficients
    })
    feature_importance_ovr_df = feature_importance_ovr_df.sort_values(by='coefficient', ascending=False).head(number_of_features)

    # Plot the top 20 important features for each class
    sns.barplot(ax=axes[idx], x='coefficient', y='feature_name', data=feature_importance_ovr_df)
    axes[idx].set_title(f'Top {number_of_features} Important Features for Class: {class_name}')

# Remove any unused subplots
for ax in axes[num_classes:]:
    ax.remove()

plt.tight_layout()

plt.savefig('figures/modeling_interpretation/ovr_feature_importance_by_author_cell_type')

# Save OVR feature importance to an excel file

# Create a DataFrame to store top 20 features for each cell type
top_features_df = pd.DataFrame(columns=['Cell Type', 'Gene', 'Coefficient'])

# Loop over each cell type
for idx, class_name in enumerate(le.classes_):
    coefficients = class_coefficients[idx]
    feature_importance_ovr_df = pd.DataFrame({
        'Gene': adata.var['feature_name'].astype(str), # are these the correct genes?
        'Coefficient': coefficients
    })

    top_features = feature_importance_ovr_df.sort_values(by='Coefficient', ascending=False).head(number_of_features)
    bottom_features = feature_importance_ovr_df.sort_values(by='Coefficient', ascending=True).head(number_of_features)
    
    # Combine top and bottom features
    combined_features = pd.concat([top_features, bottom_features])

    # Append the top features for this class to the DataFrame
    top_features_per_type_df = pd.DataFrame({
        'Cell Type': [class_name] * len(combined_features),
        'Gene': combined_features['Gene'],
        'Coefficient': combined_features['Coefficient']
    })
    
    top_features_df = pd.concat([top_features_df, top_features_per_type_df], ignore_index=True)


top_features_df.to_csv('spreadsheets/ovr_top_20_genes_by_sub_cell_type_reproduction.csv', index=False)
