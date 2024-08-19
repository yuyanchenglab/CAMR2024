#!/usr/bin/env python3
# coding: utf-8

# JM 08/17/24
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

adata = ad.read_h5ad('data/camr_scrublet_batch_filtered.h5ad')

raw = True
which_set = "Curated" #"Curated", "Coeff-0.3"
data_string = "NORMALIZED_COUNTS"
if raw:
    data_string = "RAW_COUNTS"

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


# if not os.path.isfile('spreadsheets/spreadsheets/filtered_curated_ByCell_markers_with_length_rawcounts_Coefficients-0.3.txt'):
merged_filtered_markers = pd.read_csv('spreadsheets/merged_curated-queried_ByCell_markers_with_Coefficients.txt', sep = '\t')
small_coef = np.logical_or(merged_filtered_markers["Minor_Coefficient"] <= 0.3, merged_filtered_markers["Major_Coefficient"] <= 0.3)
small_coef_in_query = np.logical_and(merged_filtered_markers["Queried"] == "Queried", small_coef)
merged_filtered_markers = merged_filtered_markers[~small_coef_in_query]
merged_filtered_markers.to_csv('spreadsheets/filtered_merged_ByCell_markers_with_Coefficients.txt', sep ='\t', index = False)

# curated entries need a queried name to work with the next section
q2n = pd.read_csv('spreadsheets/queried_to_name.txt', sep ='\t', index_col=0)
merged_filtered_markers.loc[merged_filtered_markers["Queried_Name"].isnull(), "Queried_Name"] = q2n["Queried_Name"].get(merged_filtered_markers.loc[merged_filtered_markers["Queried_Name"].isnull(), "Queried_Name"])

# raw counts
## subclass
major_sub_mean = pd.DataFrame(columns = ['Queried_Name','Marker','Raw_Mean_Expression_Subtype_Marker'])
for major_sub_class in ["AC", "BC", "Microglia", "RGC"]:
    subcell_raw_mean = pd.read_csv(f'data/raw_{major_sub_class}_subtype_meanExpression.csv')
    subcell_raw_mean.loc[subcell_raw_mean["author_cell_type"] == major_sub_class, "author_cell_type"] = "Unassigned_" + major_sub_class
    subcell_raw_mean_long = pd.melt(subcell_raw_mean, id_vars='author_cell_type', var_name='Marker', value_name='Raw_Mean_Expression_Subtype_Marker')
    subcell_raw_mean_long = subcell_raw_mean_long.rename(columns={'author_cell_type':'Queried_Name'})
    subcell_raw_mean_long["Marker"] = subcell_raw_mean_long["Marker"].str.capitalize()
    merged_filtered_markers.loc[np.logical_and(merged_filtered_markers["Queried_Name"] == major_sub_class, merged_filtered_markers["Minor_Coefficient"] > 0.3)
, "Queried_Name"] = "Unassigned_" + major_sub_class
    major_sub_mean = pd.concat([major_sub_mean, subcell_raw_mean_long], ignore_index = True)
major_sub_mean = major_sub_mean.loc[major_sub_mean['Raw_Mean_Expression_Subtype_Marker'] >= 4]
merged_filtered_markers = merged_filtered_markers.merge(major_sub_mean, on=['Queried_Name', 'Marker'], how = 'left')

# majorclass
major_raw_mean = pd.read_csv(f'data/raw_majorclass_meanExpression.csv')
major_raw_mean_long = pd.melt(major_raw_mean, id_vars='majorclass', var_name='Marker', value_name='Raw_Mean_Expression_Majorclass_Marker')
major_raw_mean_long = major_raw_mean_long.rename(columns={'majorclass':'Major_Name'})
major_raw_mean_long["Marker"] = major_raw_mean_long["Marker"].str.capitalize()
major_raw_mean_long = major_raw_mean_long.loc[major_raw_mean_long["Raw_Mean_Expression_Majorclass_Marker"] >= 4]

merged_filtered_markers = merged_filtered_markers.merge(major_raw_mean_long, on=['Major_Name', 'Marker'], how = 'left')

merged_filtered_markers.index = merged_filtered_markers["Marker"]
merged_filtered_markers = merged_filtered_markers.join(adata.var.loc[adata.var["feature_length"] >= 960, "feature_length"])
merged_filtered_markers.index = list(range(merged_filtered_markers.shape[0]))
# merged_filtered_markers["feature_length"] = merged_filtered_markers["feature_length"].to_numeric(errors = 'coerce')
# merged_filtered_markers["feature_length_div3"] = merged_filtered_markers["feature_length"].astype(int) // 3

# all majorclass
major_raw_mean = major_raw_mean.T
major_raw_mean.columns = major_raw_mean.iloc[0]
major_raw_mean["Marker"] = major_raw_mean.index.astype(str).str.capitalize()
major_raw_mean = major_raw_mean.drop(major_raw_mean.index[0])
merged_filtered_markers = merged_filtered_markers.merge(major_raw_mean, on = 'Marker', how = 'left')

merged_filtered_markers.to_csv('spreadsheets/less_filtered_curated_ByCell_markers_with_length_rawcounts_Coefficients-0.3.txt', sep ='\t', index = False)

long_enough = ~merged_filtered_markers["feature_length"].isnull()
detectable = np.logical_or(~merged_filtered_markers['Raw_Mean_Expression_Subtype_Marker'].isnull(), ~merged_filtered_markers['Raw_Mean_Expression_Majorclass_Marker'].isnull())
not_crowding = np.sum(merged_filtered_markers.loc[:, 'AC':'Rod'] > 100, axis = 1) == 0

filtered_indices = np.logical_and(long_enough, np.logical_and(detectable, not_crowding))
merged_filtered_markers = merged_filtered_markers.loc[filtered_indices]
merged_filtered_markers = merged_filtered_markers.sort_values(['Major_Name', 'Queried_Name'])
merged_filtered_markers.to_csv('spreadsheets/filtered_curated_ByCell_markers_with_length_rawcounts_Coefficients-0.3.txt', sep ='\t', index = False)

# else:
# merged_filtered_markers = pd.read_csv('spreadsheets/filtered_curated_ByCell_markers_with_length_rawcounts_Coefficients-0.3.txt', sep = '\t')


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
                  categories_order =  adata.obs.loc[adata.obs['majorclass'] == majorclass_original, "author_cell_type"].astype(str).sort_values().drop_duplicates(keep='first').tolist(),
                  vmax = 4,
                  vmin = 0,
                  show = False,
                  save = f"mouseRetina_{majorclass}_filteredMergedMarkers_{which_set}_{data_string}.pdf")
# End majorclass
