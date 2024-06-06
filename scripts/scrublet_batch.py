import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd

adata = sc.read_h5ad("d0183df5-815d-48c2-bcfe-fbf9b716505c.h5ad")

sc.pp.scrublet(adata, batch_key='reference')

adata.write(filename='camr_scrublet_batch.h5ad')