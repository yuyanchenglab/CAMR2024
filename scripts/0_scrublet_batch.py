import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd

adata = sc.read_h5ad("data/d0183df5-815d-48c2-bcfe-fbf9b716505c.h5ad")
# adata = sc.read_h5ad("data/0ca84fa3-eaa3-456e-a24c-e13e225d7ba6.h5ad")

sc.pp.scrublet(adata, batch_key='reference')

adata.write(filename='00_scrublet_batch/0_camr_scrublet_batch.h5ad')
