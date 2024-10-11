# JM 9/12/2024
# Make the adata varnames be the gene symbols rather than the Ensembl ID
# TODO: Merge this with 05

print("Setup")

import sklearn as sk
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

os.chdir('/project/ycheng11lab/jfmaurer/mouse_retina_atlas_chen_2024/')

adata = ad.read_h5ad('01_QualityControl/1_camr_scrublet_batch_filtered.h5ad')

print("Give the subtypes their designations")

adata.obs["author_cell_type"] = adata.obs["author_cell_type"].astype(str)
is_unassigned = adata.obs["author_cell_type"] == adata.obs["majorclass"]
is_subtype = adata.obs["author_cell_type"].isin(["AC", "BC", "Microglia", "RGC"])
unassigned_subtypes = adata.obs["author_cell_type"].loc[is_unassigned & is_subtype]
adata.obs.loc[is_unassigned & is_subtype, "author_cell_type"] = ["UNASSIGNED_" + uv for uv in unassigned_subtypes]
adata.obs["minorclass"] = adata.obs["author_cell_type"].astype(str)
adata.obs["Major_Name"] = adata.obs["majorclass"].astype(str)

print("Make sure we can use symbols in the shiny app")

adata.var["Ensembl"] = adata.var.index.astype(str)
# adata.var.index = adata.var["feature_name"].astype(str)
# The below error occurs when the name is not removed by `tolist()`
# AttributeError: 'Categorical' object has no attribute 'get_values'
adata.var_names = adata.var["feature_name"].astype(str).tolist()

print("Save!")

adata.write('10_Make_Shiny/10_Shiny_Input.h5ad')

print("Done!")
