import sklearn as sk
import anndata as ad
import scanpy as sc 
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

sc.settings.n_jobs = -1

adata = ad.read_h5ad('data/camr_scrublet_batch_filtered.h5ad')

fname = 'data/raw_majorclass_meanExpression.csv'
if not os.path.isfile(fname):
  raw_feature_expression_pd = pd.DataFrame(adata.raw.X.toarray(), columns = adata.var["feature_name"].astype(str).tolist())
  raw_feature_expression_pd["majorclass"] = adata.obs["majorclass"].tolist()
  raw_feature_expression_pd.groupby("majorclass").agg("mean").to_csv(fname)
