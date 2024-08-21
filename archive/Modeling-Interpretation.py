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

# Load the RandomForest classifier
rf_classifier = joblib.load('models/rf_classifier_reproduction.pkl')

# Load the LabelEncoder
le = joblib.load('models/le_reproduction.pkl')


#adata = ad.read_h5ad('camr_scrublet_batch_filtered.h5ad')
#adata = adata[:, adata.var.highly_variable]
adata = ad.read_h5ad('camr_modeling_input.h5ad')
adata


# In[ ]:


# Takes a long time to run...
# Calculating Cluster Specific Differentially Expressed Genes
#sc.tl.rank_genes_groups(adata, groupby="leiden_res_0.02", method="wilcoxon") # , use_raw = False
#adata.write('camr_modeling_input.h5ad')


# In[ ]:


#sc.pl.rank_genes_groups_dotplot(
#    adata, groupby="leiden_res_0.02", standard_scale="var", n_genes=5 # key='leiden_res_0.02'
#)

#plt.savefig('figures/modeling_interpretation/rank_gene_groups_leiden_res_0.02.png')


# ## Cell Type Classification

# ### Random Forest Classifier

# In[5]:


from scipy.sparse import issparse
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import classification_report
import seaborn as sns


from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay

# Extract feature matrix (X) and target vector (y)
X = adata.X
y = adata.obs['majorclass']

# Convert sparse matrix to dense if necessary
if issparse(X):
    X = X.toarray() # np array

# Encode the target variable
le = LabelEncoder()
y_encoded = le.fit_transform(y)

# Split data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y_encoded, test_size=0.2, random_state=42)

# This takes a verrrrry long time to run... at least an hour with 4 cores.
# For class-specific analysis, you can inspect trees or use permutation importance
# Using permutation importance from sklearn
from sklearn.inspection import permutation_importance

result = permutation_importance(rf_classifier, X_test, y_test, n_repeats=10, random_state=42, n_jobs=16)

perm_sorted_idx = result.importances_mean.argsort()

# Plot permutation importance for the top 20 features
plt.figure(figsize=(10, 8))
plt.boxplot(result.importances[perm_sorted_idx].T, vert=False, labels=np.array(adata.var_names)[perm_sorted_idx])
plt.title('Permutation Importance (test set)')

# Save fig
plt.savefig('figures/modeling_interpretation/random_forest_permutation_feature_importance_distribution')

# ### One-vs-Rest Logistic Regression
# Want to get feature importance per class

# In[36]:


# Takes less than a minute to run
from sklearn.linear_model import LogisticRegression

from scipy.sparse import issparse
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import classification_report
import seaborn as sns


# In[43]:

# Save the model to file
import joblib

# Save the RandomForest classifier
model_filename = 'models/rf_classifier_subclass_reproduction.pkl'

# Save the LabelEncoder
le_filename = 'models/le_subclass_reproduction.pkl'


# In[ ]:



# Load the RandomForest classifier
rf_classifier_subclass = joblib.load(model_filename)

# Load the LabelEncoder
le_subclass = joblib.load(le_filename)


# In[ ]:


from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay

# Encode the target variable
le_subclass = LabelEncoder()
y_encoded_subclass = le_subclass.fit_transform(y_subclass)

# Split data into training and testing sets
X_train, X_test, y_train_subclass, y_test_subclass = train_test_split(
    X, y_encoded_subclass, test_size=0.2, random_state=42)

# Make predictions on the test set
y_pred_subclass_rf = rf_classifier_subclass.predict(X_test)

# Verify the number of unique classes
num_classes_in_y = len(np.unique(y_encoded_subclass))
num_classes_in_le = len(le_subclass.classes_)

target_names = le_subclass.inverse_transform(np.unique(y_test_subclass))

# TODO sort by cell count
print(classification_report(y_test_subclass, y_pred_subclass_rf, target_names=target_names))

# Generate the confusion matrix
cm = confusion_matrix(y_test, y_pred_subclass_rf, labels=rf_classifier_subclass.classes_)

# Make percentage
cm_normalized = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]

# Plot the confusion matrix
plt.figure(figsize=(10, 9))
sns.heatmap(cm_normalized, cmap='Blues', 
            xticklabels=le_subclass.inverse_transform(rf_classifier_subclass.classes_), 
            yticklabels=le_subclass.inverse_transform(rf_classifier_subclass.classes_))

# Create Labels
plt.xlabel('Predicted Label')
plt.ylabel('True Label')
plt.title('Confusion Matrix Subclass (author_cell_type, Random Forest)')

# Save Figure to File
plt.savefig("figures/modeling_interpretation/random_forest_subclass_confusion_matrix", bbox_inches='tight')
plt.show()


# In[ ]:


# Testing model manually

# Define the index of the cell you want to test
cell_index = 1025  # Change this to the desired cell index

# Extract the single cell's features
single_cell = X[cell_index].reshape(1, -1)

# Predict the log-probabilities for the single cell
log_proba = rf_classifier_subclass.predict_log_proba(single_cell)

# Get the top 10 classes with the highest log-probabilities
top_10_indices = np.argsort(log_proba[0])[-10:][::-1]
top_10_log_proba = log_proba[0][top_10_indices]
top_10_classes = le_subclass.inverse_transform(top_10_indices)

# Predict the class for the single cell
predicted_class_encoded = rf_classifier_subclass.predict(single_cell)

# Decode the predicted class
predicted_class = le_subclass.inverse_transform(predicted_class_encoded)

# Get the true class for the single cell
true_class = y_subclass.iloc[cell_index]

# Print the result
print(f'Predicted class for cell index {cell_index}: {predicted_class[0]}')
print(f'True class for cell index {cell_index}: {true_class}')

# Print the top 10 log-probabilities for each class
print(f'Top 10 log-probabilities for cell index {cell_index}:')
for class_name, log_prob in zip(top_10_classes, top_10_log_proba):
    print(f'Class: {class_name}, Log-probability: {log_prob:.4f}')


# In[ ]:


import numpy as np
from sklearn.tree import _tree, export_text

# Extract the single cell's features
single_cell = X[cell_index].reshape(1, -1)

# Select the most important tree (for example, the first tree)
tree_index = 0  # Change this to select a different tree
tree = rf_classifier_subclass.estimators_[tree_index]

# Get the decision path for the single cell in the selected tree
node_indicator = tree.decision_path(single_cell)
leaf_id = tree.apply(single_cell)

# Extract the feature and threshold information
feature = tree.tree_.feature
threshold = tree.tree_.threshold

# Get the node indices along the decision path
node_index = node_indicator.indices[node_indicator.indptr[0]:node_indicator.indptr[1]]

# Print the decision path
print(f'Decision path for cell index {cell_index} in tree {tree_index}:')
for node_id in node_index:
    if leaf_id[0] == node_id:
        print(f"Node {node_id}: leaf node.")
    else:
        # Get feature name (gene) and threshold
        feature_name = adata.var['feature_name'][feature[node_id]]
        threshold_value = threshold[node_id]
        if single_cell[0, feature[node_id]] <= threshold_value:
            threshold_sign = "<="
        else:
            threshold_sign = ">"
        print(f"Node {node_id}: split on {feature_name} {threshold_sign} {threshold_value:.4f}")

# Optional: Print the text representation of the tree for reference
tree_rules = export_text(tree, feature_names=adata.var['feature_name'])
print(tree_rules)


# In[ ]:


le_subclass.inverse_transform([128])


# In[ ]:


# Get feature importances
feature_importances_subclass_rf = rf_classifier_subclass.feature_importances_

# Create a DataFrame for better visualization
feature_importance_subclass_rf_df = pd.DataFrame({
    'ESUM': adata.var_names,
    'feature_name': adata.var['feature_name'].astype(str),  # Ensure feature_name is string
    'Importance': feature_importances_subclass_rf
})

feature_importance_subclass_rf_df = feature_importance_subclass_rf_df.sort_values(
    by='Importance', ascending=False).reset_index(drop=True)

# Plot the top 20 important features overall
plt.figure(figsize=(10, 8))
sns.barplot(x='Importance', y='feature_name', data=feature_importance_subclass_rf_df.head(20))
plt.title('Top 20 Important Features Random Forest Subclass')

# Save Fig
plt.savefig('figures/modeling_interpretation/random_forest_subclass_top_20_feature_importance')
plt.show()


# In[ ]:


feature_importance_subclass_rf_df


# In[ ]:


# Plot the distribution of feature importances
plt.figure(figsize=(10, 6))
sns.histplot(feature_importance_subclass_rf_df['Importance'], bins=30, kde=True)
plt.title('Distribution of Feature Importances Random Forest Subclass')
plt.xlabel('Importance')
plt.ylabel('Frequency')

#Save figure
plt.savefig('figures/modeling_interpretation/random_forest_subclass_feature_importance_distribution')

plt.show()


# In[ ]:


# Calculate the cumulative sum of importances
feature_importance_subclass_rf_df['Cumulative Importance'] = feature_importance_subclass_rf_df['Importance'].cumsum()

# Plot the cumulative sum
plt.figure(figsize=(10, 8))
sns.lineplot(data=feature_importance_subclass_rf_df, x=feature_importance_subclass_rf_df.index, y='Cumulative Importance')
plt.axhline(y=0.9, color='r', linestyle='--')  # Line at 90% cumulative importance
plt.xticks(ticks=range(0, len(feature_importance_subclass_rf_df), 500))  # Set x-ticks at increments of 100
plt.xlabel('Feature Rank (sorted by importance)')
plt.ylabel('Cumulative Importance')
plt.title('Cumulative Sum of Feature Importances Subclass')

#save fig
plt.savefig('figures/modeling_interpretation/random_forest_subclass_cumulative_feature_importance_distribution')

plt.show()


# ### One vs Rest

# In[ ]:


from sklearn.linear_model import LogisticRegression

ovr_classifier_subclass = LogisticRegression(multi_class='ovr', max_iter=1000, random_state=42)
ovr_classifier_subclass.fit(X_train, y_train_subclass)


# In[ ]:


y_pred_subclass_ovr = ovr_classifier.predict(X_test)
print(classification_report(y_test, y_pred_subclass_ovr, target_names=le.classes_))


# In[ ]:


# Get class-specific coefficients
class_coefficients = ovr_classifier.coef_
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

plt.savefig('figures/modeling_interpretation/ovr_feature_importance_by_majorclass')

plt.show()


# ## Quality of Life

# In[ ]:


import dill

# Save the entire session
filename = 'modeling_notebook_state.pkl'
with open(filename, 'wb') as f:
    dill.dump_session(f)

