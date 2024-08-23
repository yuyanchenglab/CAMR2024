import sklearn as sk
import anndata as ad
import scanpy as sc 
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

sc.settings.n_jobs = -1

adata = ad.read_h5ad('01_QualityControl/1_camr_scrublet_batch_filtered.h5ad')
genes = adata.var["feature_name"].astype(str).tolist()

fmajorname = 'data/raw_meanExpression_majorclass.txt'
if not os.path.isfile(fname):
  raw_feature_expression_pd = pd.DataFrame(adata.raw.X.toarray(), columns = genes)
  raw_feature_expression_pd["majorclass"] = adata.obs["majorclass"].tolist()
  raw_feature_expression_pd.groupby("majorclass").agg("mean").to_csv(fmajorname, sep = '\t')

fminorname = 'data/raw_meanExpression_minorclass.txt'
if not os.path.isfile(fminorname):
  raw_feature_expression_pd = pd.DataFrame(adata.raw.X.toarray(), columns = genes)
  raw_feature_expression_pd["author_cell_type"] = adata.obs["author_cell_type"].tolist()
  raw_feature_expression_pd.groupby("author_cell_type").agg("mean").to_csv(fminorname, sep = '\t')