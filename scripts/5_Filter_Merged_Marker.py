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
os.makedirs('05_Filter_Merged_Markers', exist_ok = True)
sc.settings.n_jobs = -1

adata = ad.read_h5ad('01_QualityControl/1_camr_scrublet_batch_filtered.h5ad')
subcell_raw_mean = pd.read_csv('data/raw_meanExpression_minorclass.txt', sep = '\t')

sc.plotting.DotPlot.DEFAULT_SAVE_PREFIX = "05_Filter_Merged_Markers/5_dotplot_"
sc.plotting.DotPlot.DEFAULT_LARGEST_DOT = 200.0

raw = True
data_string = "normCounts"
max_col = 3
if raw:
    data_string = "rawCounts"
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


merged_filtered_markers = pd.read_csv('04_Merge_Curated_Markers/4_merged_curated-queried_markers.txt', sep = '\t')
small_coef = np.logical_or(merged_filtered_markers["Minor_Coefficient"] <= 0.3, merged_filtered_markers["Major_Coefficient"] <= 0.3)
small_coef_in_query = np.logical_and(merged_filtered_markers["Queried"] == "Queried", small_coef)
merged_filtered_markers = merged_filtered_markers[~small_coef_in_query]
merged_filtered_markers.to_csv('05_Filter_Merged_Markers/5_merged_curated-queried_markers_coefficientFiltered.txt', sep ='\t', index = False)

# Add gene length
merged_filtered_markers.index = merged_filtered_markers["Marker"]
merged_filtered_markers = merged_filtered_markers.join(adata.var.loc[adata.var["feature_length"] >= 960, "feature_length"])
merged_filtered_markers.index = list(range(merged_filtered_markers.shape[0]))

# Curated entries need a queried name to work with the next section
q2n = pd.read_csv('04_Merge_Curated_Markers/4_queried_to_name.txt', sep ='\t', index_col=0)
merged_filtered_markers.loc[merged_filtered_markers["Queried_Name"].isnull(), "Queried_Name"] = q2n["Queried_Name"].get(merged_filtered_markers.loc[merged_filtered_markers["Queried_Name"].isnull(), "Queried_Name"])

# raw counts
## minorclass
### Pre-process names
is_subtype = subcell_raw_mean.loc[subcell_raw_mean["author_cell_type"].isin(["AC", "BC", "Microglia", "RGC"])
unassigned_subtypes = subcell_raw_mean.loc[is_subtype, "author_cell_type"]
subcell_raw_mean.loc[is_subtype, "author_cell_type"] = ["Unassigned_" + uv for uv in unassigned_subtypes]

### Process table
subcell_raw_mean_long = pd.melt(subcell_raw_mean, id_vars='author_cell_type', var_name='Marker', value_name='Raw_Mean_Expression_Minorclass_Marker')
subcell_raw_mean_long = subcell_raw_mean_long.rename(columns={'author_cell_type':'Queried_Name'})
subcell_raw_mean_long["Marker"] = subcell_raw_mean_long["Marker"].str.capitalize()
major_sub_mean = subcell_raw_mean_long.loc[subcell_raw_mean_long['Raw_Mean_Expression_Minorclass_Marker'] >= 4]

### Add expression data
merged_filtered_markers = merged_filtered_markers.merge(major_sub_mean, on=['Queried_Name', 'Marker'], how = 'left')

## majorclass
major_raw_mean = pd.read_csv('data/raw_meanExpression_majorclass.csv')
major_raw_mean_long = pd.melt(major_raw_mean, id_vars='majorclass', var_name='Marker', value_name='Raw_Mean_Expression_Majorclass_Marker')
major_raw_mean_long = major_raw_mean_long.rename(columns={'majorclass':'Major_Name'})
major_raw_mean_long["Marker"] = major_raw_mean_long["Marker"].str.capitalize()
major_raw_mean_long = major_raw_mean_long.loc[major_raw_mean_long["Raw_Mean_Expression_Majorclass_Marker"] >= 4]

merged_filtered_markers = merged_filtered_markers.merge(major_raw_mean_long, on=['Major_Name', 'Marker'], how = 'left')

## all majorclass for each row
major_raw_mean = major_raw_mean.T
major_raw_mean.columns = major_raw_mean.iloc[0]
major_raw_mean["Marker"] = major_raw_mean.index.astype(str).str.capitalize()
major_raw_mean = major_raw_mean.drop(major_raw_mean.index[0])
merged_filtered_markers = merged_filtered_markers.merge(major_raw_mean, on = 'Marker', how = 'left')

# 
long_enough = ~merged_filtered_markers["feature_length"].isnull()
detectable = np.logical_or(~merged_filtered_markers['Raw_Mean_Expression_Minorclass_Marker'].isnull(), ~merged_filtered_markers['Raw_Mean_Expression_Majorclass_Marker'].isnull())
not_crowding = np.sum(merged_filtered_markers.loc[:, 'AC':'Rod'] > 100, axis = 1) == 0

filtered_indices = np.logical_and(long_enough, np.logical_and(detectable, not_crowding))
merged_filtered_markers = merged_filtered_markers.loc[filtered_indices]
merged_filtered_markers = merged_filtered_markers.sort_values(['Major_Name', 'Name']) # For future
merged_filtered_markers.to_csv('05_Filter_Merged_Markers/5_merged_curated-queried_markers_sorted_coefficientLengthExpressionFiltered.txt', sep ='\t', index = False)
merged_filtered_markers = merged_filtered_markers.sort_values(['Major_Name', 'Queried_Name']) # For now
merged_filtered_markers.to_csv('05_Filter_Merged_Markers/5_merged_curated-queried_markers_coefficientLengthExpressionFiltered.txt', sep = '\t')


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
                  save = f"mouseRetina_minorclass-{majorclass}_filteredMergedMarkers_{data_string}.pdf")
# End majorclass

majorclass_idx = adata.obs["author_cell_type"].str.upper().isin(adata.obs["majorclass"].astype(str).str.upper())
adata.obs.loc[majorclass_idx, "author_cell_type"] = adata.obs.loc[majorclass_idx, "author_cell_type"].str.upper()

merged_filtered_markers = merged_filtered_markers.loc[merged_filtered_markers["Curated"] == "Curated"]
merged_filtered_markers1 = merged_filtered_markers.loc[merged_filtered_markers["Major_Name"] == "AC"]
merged_filtered_markers2 = merged_filtered_markers.loc[merged_filtered_markers["Major_Name"] == "BC"]
merged_filtered_markers3 = merged_filtered_markers.loc[merged_filtered_markers["Major_Name"] == "RGC"]
merged_filtered_markers4 = merged_filtered_markers.loc[merged_filtered_markers["Major_Name"] == "MG"]
merged_filtered_markers5 = merged_filtered_markers.loc[merged_filtered_markers["Major_Name"] == "RPE"]
merged_filtered_markers = pd.concat([merged_filtered_markers1,merged_filtered_markers2,merged_filtered_markers3,merged_filtered_markers4,merged_filtered_markers5], ignore_index = True)

adata.obs.loc[adata.obs["author_cell_type"].astype(str).str.contains("Microglia"), "author_cell_type"] = "MICROGLIA" # No subtypes of Microglia here
ordered_celltypes = merged_filtered_markers["Queried_Name"].drop_duplicates().astype('category').cat.remove_categories('RGC').dropna().tolist() + ["ROD","CONE","HC","MICROGLIA","ENDOTHELIAL"]

final_ordered_markers = ["Ptn","Cntn6","Nxph1","Cpne4","Cbln4","Etv1","Epha3","Trhde", # AC
                        "Neto1","Erbb4","Grik1","Nnat","Cabp5","Sox6","Cck","Cpne9","Prkca","Hcn1", # BC
                        "Rbpms","Penk","Spp1","Coch","Opn4","Il1rapl2","Tbx20","Tac1","Cartpt", # RGC
                        "Aqp4","Glul","Rlbp1","Slc1a3","Vim", # MG
                        "Rpe65", # RPE
                        "Rho", "Arr3","Lhx1","Onecut1", "Cd74", "Cldn5"] # Manual

sc.pl.dotplot(adata[adata.obs['author_cell_type'].isin(ordered_celltypes), :],
              var_names = final_ordered_markers,
              groupby = 'author_cell_type',
              categories_order = ordered_celltypes,
              vmax = max_col,
              vmin = 0,
              show = False,
              figsize = (12,10),
              save = f"mouseRetina_filteredMergedMarkers_curated_{data_string}.pdf")
