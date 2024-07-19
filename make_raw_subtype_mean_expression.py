import sklearn as sk
import anndata as ad
import scanpy as sc 
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

sc.settings.n_jobs = -1

adata = ad.read_h5ad('camr_modeling_input.h5ad')

if not os.path.isfile('raw_subtype_mean_variable_genes.csv'):
  highly_variable = adata.raw.var['feature_name'].isin(adata.var['feature_name'])
  raw_feature_expression_pd = pd.DataFrame(adata.raw.X[:, highly_variable].toarray(), columns = adata.var["feature_name"].astype(str).tolist())
  raw_feature_expression_pd["author_cell_type"] = adata.obs["author_cell_type"].tolist()
  raw_feature_expression_pd_mean = raw_feature_expression_pd.groupby("author_cell_type").agg("mean")
  raw_feature_expression_pd_mean.to_csv('raw_subtype_mean_variable_genes.csv')
