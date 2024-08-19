#!/usr/bin/env python3
# coding: utf-8

# JM 08/18/24
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
sc.settings.n_jobs = -1

import matplotlib
matplotlib.rcParams.update({'font.size': 18})

adata = ad.read_h5ad('data/camr_scrublet_batch_filtered.h5ad')

raw = True
which_set = "CuratedONLY" #"Curated", "Coeff-0.3"
data_string = "NORMALIZED_COUNTS"
max_col = 3
if raw:
    data_string = "RAW_COUNTS"
    max_col = 4

# Clean subtypes
adata.obs = adata.obs.loc[:, ["majorclass", "author_cell_type"]]
adata.obs["author_cell_type"] = adata.obs["author_cell_type"].astype(str)
is_unassigned = adata.obs["author_cell_type"] == adata.obs["majorclass"]
is_subtype = adata.obs["author_cell_type"].isin(["AC", "BC", "Microglia", "RGC"])
unassigned_subtypes = adata.obs["author_cell_type"].loc[is_unassigned & is_subtype]
adata.obs.loc[is_unassigned & is_subtype, "author_cell_type"] = ["Unassigned_" + uv for uv in unassigned_subtypes]

majorclass_idx = adata.obs["author_cell_type"].str.upper().isin(adata.obs["majorclass"].astype(str).str.upper())
adata.obs.loc[majorclass_idx, "author_cell_type"] = adata.obs.loc[majorclass_idx, "author_cell_type"].str.upper()

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

merged_filtered_markers = pd.read_csv('spreadsheets/filtered_curated_ByCell_markers_with_length_rawcounts_Coefficients-0.3.txt', sep = '\t')
merged_filtered_markers = merged_filtered_markers.loc[merged_filtered_markers["Curated"] == "Curated"]
merged_filtered_markers1 = merged_filtered_markers.loc[merged_filtered_markers["Major_Name"] == "AC"]
merged_filtered_markers2 = merged_filtered_markers.loc[merged_filtered_markers["Major_Name"] == "BC"]
merged_filtered_markers3 = merged_filtered_markers.loc[merged_filtered_markers["Major_Name"] == "RGC"]
merged_filtered_markers4 = merged_filtered_markers.loc[merged_filtered_markers["Major_Name"] == "MG"]
merged_filtered_markers5 = merged_filtered_markers.loc[merged_filtered_markers["Major_Name"] == "RPE"]
merged_filtered_markers = pd.concat([merged_filtered_markers1,merged_filtered_markers2,merged_filtered_markers3,merged_filtered_markers4,merged_filtered_markers5], ignore_index = True)
# merged_filtered_markers.to_csv('spreadsheets/filtered_curatedONLY_ByCell_markers_with_length_rawcounts_Coefficients-0.3.txt', sep = '\t', index=False)
# unique_markers = []
# final_markers = []
# for majorclass in adata.obs['majorclass'].cat.categories:
    
#     print(f'{datetime.datetime.now()} Major Class: {majorclass}')
    
#     majorclass_original = majorclass
#     majorclass = majorclass.upper()
    
#     is_major_marker = np.logical_and(merged_filtered_markers["Major_Name"] == majorclass, merged_filtered_markers["Name"] == majorclass)
#     cell_markers = merged_filtered_markers.loc[is_major_marker, "Marker"].tolist()
    
#     is_subtype_marker = np.logical_and(merged_filtered_markers["Major_Name"] == majorclass, merged_filtered_markers["Name"] != majorclass)
#     subtype_markers = merged_filtered_markers.loc[is_subtype_marker, "Marker"].tolist()
    
#     markers = cell_markers + subtype_markers
    
#     marker2type = merged_filtered_markers.loc[is_major_marker, "Name"].tolist() + merged_filtered_markers.loc[is_subtype_marker, "Name"].tolist()
    
#     # sc.pl.dotplot throws a fit if there are duplicates
    
#     for i, m in enumerate(markers):
#         if m not in unique_markers:
#             unique_markers += [m]
#         else:
#             print(f'Found duplicate marker {m} in {marker2type[i]}!')

#     if markers == []: # Necessary to avoid "ValueError: left cannot be >= right"
#         print(f'No markers available for {majorclass}!')
#         continue

# for m in unique_markers: # Necessary to avoid "KeyError: "Values ['Cd39'], from ['Calb2', ...], are not valid obs/ var names or indices.
#     if m in adata.var_names:
#         final_markers += [m]
#     else:
#         print(f'Marker {m} from curate is not in this data!')

# print(f'Final Markers: {final_markers}')
adata.obs.loc[adata.obs["author_cell_type"].astype(str).str.contains("Microglia"), "author_cell_type"] = "MICROGLIA"
ordered_celltypes = merged_filtered_markers["Queried_Name"].drop_duplicates().astype('category').cat.remove_categories('RGC').dropna().tolist() + ["ROD","CONE","HC","MICROGLIA","ENDOTHELIAL"]
# End majorclass    
sc.pl.dotplot(adata[adata.obs['author_cell_type'].isin(ordered_celltypes), :],
              var_names = ["Ptn","Cntn6","Nxph1","Cpne4","Cbln4","Etv1","Epha3","Trhde", # AC
                        "Neto1","Erbb4","Grik1","Nnat","Cabp5","Sox6","Cck","Cpne9","Prkca","Hcn1", # BC
                          "Rbpms","Penk","Spp1","Coch","Opn4","Il1rapl2","Tbx20","Tac1","Cartpt", # RGC
                          "Aqp4","Glul","Rlbp1","Slc1a3","Vim", # MG
                          "Rpe65", "Rho", "Arr3","Lhx1","Onecut1", "Cd74", "Cldn5"], # RPE
              groupby = 'author_cell_type',
              categories_order = ordered_celltypes,
              vmax = max_col,
              vmin = 0,
              show = False,
              # s = 3,
              figsize = (12,10),
              save = f"BIG_mouseRetina_together_filteredMergedMarkers_{which_set}_{data_string}.pdf")

