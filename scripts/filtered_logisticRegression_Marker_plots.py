#!/usr/bin/env python3
# coding: utf-8

# JM 08/08/24
# In this notebook, we're going to re-apply machine learning to analyze the mouse retinal data and determine proper marker genes.

import datetime
print(f'{datetime.datetime.now()} Analysis Setup')

import sklearn as sk
import anndata as ad
import scanpy as sc 
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import os

os.chdir('/project/hipaa_ycheng11lab/atlas/CAMR2024')
sc.settings.n_jobs = -1

ncols = 2

adata = ad.read_h5ad('data/camr_scrublet_batch_filtered.h5ad')

# Make memory available
adata.raw = None
adata.obs = adata.obs.loc[:, ["majorclass", "author_cell_type"]]
adata.var = adata.var.loc[:, ["gene_symbols", "feature_name"]]

adata.var.index = [gene.title() for gene in adata.var["feature_name"]] # subset on genes instead of booleans
adata.var_names = [gene.title() for gene in adata.var["feature_name"]] # subset on genes instead of booleans

adata.var_names_make_unique()

import gc
import ctypes
gc.collect() # Free memory
libc = ctypes.CDLL("libc.so.6") # clearing cache 
libc.malloc_trim(0)

if not os.path.isfile('spreadsheets/filtered_merged_ByCell_markers_with_Coefficients.txt'):
    merged_filtered_markers = pd.read_csv('spreadsheets/merged_curated-queried_ByCell_markers_with_Coefficients.txt', sep = '\t')
    small_coef = np.logical_or(merged_filtered_markers["Minor_Coefficient"] <= 0.3, merged_filtered_markers["Major_Coefficient"] <= 0.3)
    small_coef_in_query = np.logical_and(merged_filtered_markers["Queried"] == "Queried", small_coef)
    merged_filtered_markers = merged_filtered_markers[~small_coef_in_query]
    merged_filtered_markers.to_csv('spreadsheets/filtered_merged_ByCell_markers_with_Coefficients.txt', sep ='\t', index = False)
else:
    merged_filtered_markers = pd.read_csv('spreadsheets/filtered_merged_ByCell_markers_with_Coefficients.txt', sep = '\t')

for majorclass in adata.obs['majorclass'].cat.categories:
    
    print(f'{datetime.datetime.now()} Major Class: {majorclass}')
    
    is_majorclass = adata.obs['majorclass'].isin([majorclass]).tolist()
    majorclass = majorclass.upper()
    
    is_major_marker = np.logical_and(merged_filtered_markers["Major_Name"] == majorclass, merged_filtered_markers["Name"] == majorclass)
    cell_markers = merged_filtered_markers.loc[is_major_marker, "Marker"].tolist()
    
    is_subtype_marker = np.logical_and(merged_filtered_markers["Major_Name"] == majorclass, merged_filtered_markers["Name"] != majorclass)
    subtype_markers = merged_filtered_markers.loc[is_subtype_marker, "Marker"].tolist()
    
    markers = cell_markers + subtype_markers
    
    marker2type = merged_filtered_markers.loc[is_major_marker, "Name"].tolist() + merged_filtered_markers.loc[is_subtype_marker, "Name"].tolist()
    
    # sc.pl.dotplot throws a fit if there are duplicates
    unique_markers = []
    for i, m in enumerate(markers):
        if m not in unique_markers:
            unique_markers += [m]
        else:
            print(f'Found duplicate marker {m} in {marker2type[i]}!')

    final_markers = []
    for m in unique_markers:
        if m in adata.var_names:
            final_markers += [m]
        else:
            print(f'Found marker {m} from curate that is not in this data!')
    
    sc.pl.dotplot(adata[is_majorclass, final_markers],
                  var_names = final_markers,
                  groupby = 'author_cell_type',
                  vmax = 12,
                  vmin = 3,
                  show = False,
                  save = f"mouseRetina_{majorclass}_filteredMergedMarkers_Coeff-0.3.pdf")
# End majorclass
