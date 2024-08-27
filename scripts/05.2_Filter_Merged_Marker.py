#!/usr/bin/env python3
# coding: utf-8

# In this notebook, we're going plot the mouse retinal data and determine proper marker genes.

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
os.makedirs('05_Filter_Merged_Markers', exist_ok = True)
sc.settings.n_jobs = -1
if True:
    adata = ad.read_h5ad('01_QualityControl/1_camr_scrublet_batch_filtered.h5ad')
else:
    adata = ad.read_h5ad('01_QualityControl/1_camr_scrublet_batch_filtered.h5ad', backed='r')

sc.plotting.DotPlot.DEFAULT_SAVE_PREFIX = "" # "05_Filter_Merged_Markers/figures/5_dotplot_" # figures/05_Filter_Merged_Markers/figures/5_dotplot_mouseRetina_minorclass-ROD_curatedMarkers_rawCounts.pdf
sc.plotting.DotPlot.DEFAULT_LARGEST_DOT = 200.0

raw = False 
data_string = "normCounts" # 11GB
max_col = 3
if raw:
    data_string = "rawCounts" # 11GB
    max_col = 4

# Clean subtypes
adata.obs = adata.obs.loc[:, ["majorclass", "author_cell_type"]]
adata.obs["author_cell_type"] = adata.obs["author_cell_type"].astype(str)
is_unassigned = adata.obs["author_cell_type"] == adata.obs["majorclass"]
is_subtype = adata.obs["author_cell_type"].isin(["AC", "BC", "Microglia", "RGC"])
unassigned_subtypes = adata.obs["author_cell_type"].loc[is_unassigned & is_subtype]
adata.obs.loc[is_unassigned & is_subtype, "author_cell_type"] = ["Unassigned_" + uv for uv in unassigned_subtypes]


# Make memory available
if raw:
    adata.X = adata.raw.X
adata.raw = None
adata.var = adata.var.loc[:, ["gene_symbols", "feature_name", "feature_length"]]
adata.var["feature_name"] = adata.var["feature_name"].astype(str).str.capitalize()
adata.var["feature_length"] = adata.var["feature_length"].astype(int)
adata.var.index = adata.var["feature_name"] # subset on genes instead of booleans
adata.var_names = adata.var["feature_name"] # subset on genes instead of booleans
# adata.var_names_make_unique()

import gc
import ctypes
gc.collect() # Free memory
libc = ctypes.CDLL("libc.so.6") # clearing cache 
libc.malloc_trim(0)


merged_filtered_markers = pd.read_csv('05_Filter_Merged_Markers/5_curated_markers_annotatedKeep_lengthExpressionMissingFiltered.txt', sep ='\t')
merged_filtered_markers = merged_filtered_markers.sort_values(['Major_Name', 'Queried_Name']) # For now
# merged_filtered_markers.to_csv('05_Filter_Merged_Markers/5_merged_curated-queried_markers_coefficientLengthExpressionFiltered.txt', sep = '\t')


for majorclass in adata.obs['majorclass'].cat.categories:
    
    print(f'{datetime.datetime.now()} Major Class: {majorclass}')
    
    majorclass_original = majorclass
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
    for m in unique_markers: # Necessary to avoid "KeyError: "Values ['Cd39'], from ['Calb2', ...], are not valid obs/ var names or indices.
        if m in adata.var_names:
            final_markers += [m]
        else:
            print(f'Marker {m} from curate is not in this data!')

    print(f'Final Markers: {final_markers}')

    if final_markers == []: # Necessary to avoid "ValueError: left cannot be >= right"
        print(f'No markers available for {majorclass}!')
        continue
    
    sc.pl.dotplot(adata[adata.obs['majorclass'] == majorclass_original, final_markers],
                  var_names = final_markers,
                  groupby = 'author_cell_type',
                  categories_order = adata.obs.loc[adata.obs['majorclass'] == majorclass_original, "author_cell_type"].astype(str).sort_values().drop_duplicates(keep='first').tolist(), # Only celltypes that have a marker should be present
                  vmax = max_col,
                  vmin = 0,
                  show = False,
                  save = f"/project/hipaa_ycheng11lab/atlas/CAMR2024/05_Filter_Merged_Markers/figures/5_dotplot_mouseRetina_minorclass-{majorclass}_curatedMarkers_{data_string}.pdf")
# End majorclass
