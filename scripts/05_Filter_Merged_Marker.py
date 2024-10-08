#!/usr/bin/env python3
# coding: utf-8

# JM 08/17/24 # This is going to have the merged
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

os.chdir('/project/ycheng11lab/jfmaurer/mouse_retina_atlas_chen_2024/')
os.makedirs('05_Filter_Merged_Markers', exist_ok = True)
sc.settings.n_jobs = -1

adata = ad.read_h5ad('01_QualityControl/1_camr_scrublet_batch_filtered.h5ad')

sc.plotting.DotPlot.DEFAULT_SAVE_PREFIX = "05_Filter_Merged_Markers/figures/5_dotplot_"
sc.plotting.DotPlot.DEFAULT_LARGEST_DOT = 200.0

coefficient_threshold = 0.3
length_threshold = 960
raw_expression_threshold = 4

plot_major_markers = False
plot_major_cells = True
plot_minor_markers = True
plot_minor_cells = True
plot_split_majorclass = True

raw = True
xenium_filtered = True

data_string = "normCounts"
max_col = 3
if raw:
    data_string = "rawCounts"
    max_col = 4
    adata.X = adata.raw.X

if xenium_filtered:
    data_string = data_string + "_xeniumFiltered"

plot_occassion = ""
data_string = data_string + f"_{plot_occassion}"
markers = ""
if plot_occassion == "august_grant": # 05.1
    pass
elif plot_occassion == "curated_xenium_filtered": # 05.2
    pass
else:
    pass

# Clean subtypes
adata.obs = adata.obs.loc[:, ["majorclass", "author_cell_type"]]
adata.obs["author_cell_type"] = adata.obs["author_cell_type"].astype(str)
is_unassigned = adata.obs["author_cell_type"] == adata.obs["majorclass"]
is_subtype = adata.obs["author_cell_type"].isin(["AC", "BC", "Microglia", "RGC"])
unassigned_subtypes = adata.obs["author_cell_type"].loc[is_unassigned & is_subtype]
adata.obs.loc[is_unassigned & is_subtype, "author_cell_type"] = ["Unassigned_" + uv for uv in unassigned_subtypes]
adata.obs["minorclass"] = adata.obs["author_cell_type"].astype(str)
adata.obs["Name"] = adata.obs["author_cell_type"].astype(str) # This should be fed through q2n
adata.obs["Major_Name"] = adata.obs["majorclass"].astype(str) # This should be fed through q2n

# Shrink adata footprint
adata.raw = None
adata.var = adata.var.loc[:, ["gene_symbols", "feature_name", "feature_length"]]
adata.var["Ensembl"] = adata.var.index.astype(str)
adata.var["feature_name"] = adata.var["feature_name"].astype(str).str.capitalize()
adata.var["feature_length"] = adata.var["feature_length"].astype(int)
adata.var_names = adata.var["feature_name"].astype(str).tolist() # No need to manually adjust index
# adata.var_names_make_unique() # Why is this like this?


# Work starts here


# Coefficient Marking
merged_filtered_markers = pd.read_csv('04_Merge_Curated_Markers/4_merged_curated-queried_markers.txt', sep = '\t').drop_duplicates()
small_coef = np.logical_or(merged_filtered_markers["Minor_Coefficient"] <= coefficient_threshold, merged_filtered_markers["Major_Coefficient"] <= coefficient_threshold)
merged_filtered_markers["small_coef_in_query"] = np.logical_and(merged_filtered_markers["Queried"] == "Queried", small_coef)

# Add gene length
merged_filtered_markers.index = merged_filtered_markers["Marker"] # Make sure the below join keeps all rows
merged_filtered_markers = merged_filtered_markers.join(adata.var)
merged_filtered_markers["missing_length"] = merged_filtered_markers["feature_length"].isnull()
merged_filtered_markers[merged_filtered_markers["missing_length"], "feature_length"] = 0
merged_filtered_markers["long_enough"] = merged_filtered_markers["feature_length"] >= length_threshold
merged_filtered_markers.index = list(range(merged_filtered_markers.shape[0]))

# CURATED entries need a queried name to work with the next section
q2n = pd.read_csv('04_Merge_Curated_Markers/4_queried_to_name.txt', sep ='\t', index_col=0)
merged_filtered_markers.loc[merged_filtered_markers["Queried_Name"].isnull(), "Queried_Name"] = q2n["Queried_Name"].get(merged_filtered_markers.loc[merged_filtered_markers["Queried_Name"].isnull(), "Queried_Name"])

# minorclass
## Pre-process expression table names
subcell_raw_mean = pd.read_csv('data/raw_meanExpression_minorclass.txt', sep = '\t')
is_subtype = subcell_raw_mean["author_cell_type"].isin(["AC", "BC", "Microglia", "RGC"])
unassigned_subtypes = subcell_raw_mean.loc[is_subtype, "author_cell_type"]
subcell_raw_mean.loc[is_subtype, "author_cell_type"] = ["Unassigned_" + uv for uv in unassigned_subtypes]

## Process table
subcell_raw_mean_long = pd.melt(subcell_raw_mean, id_vars='author_cell_type', var_name='Marker', value_name='Raw_Mean_Expression_Minorclass_Target')
subcell_raw_mean_long = subcell_raw_mean_long.rename(columns={'author_cell_type':'Queried_Name'})
subcell_raw_mean_long["Marker"] = subcell_raw_mean_long["Marker"].str.capitalize()
subcell_raw_mean_long["Minor_Detectable_Expression"] = subcell_raw_mean_long['Raw_Mean_Expression_Minorclass_Target'] >= raw_expression_threshold
merged_filtered_markers = merged_filtered_markers.merge(subcell_raw_mean_long, on=['Queried_Name', 'Marker'], how = 'left')

# majorclass expression
major_raw_mean = pd.read_csv('data/raw_meanExpression_majorclass.txt', sep = '\t')
major_raw_mean_long = pd.melt(major_raw_mean, id_vars='majorclass', var_name='Marker', value_name='Raw_Mean_Expression_Majorclass_Target')
major_raw_mean_long = major_raw_mean_long.rename(columns={'majorclass':'Major_Name'}) # Major_Queried name?
major_raw_mean_long["Marker"] = major_raw_mean_long["Marker"].str.capitalize()
major_raw_mean_long["Major_Detectable_Expression"] = major_raw_mean_long["Raw_Mean_Expression_Majorclass_Target"] >= raw_expression_threshold
merged_filtered_markers = merged_filtered_markers.merge(major_raw_mean_long, on=['Major_Name', 'Marker'], how = 'left')

## all majorclass expression
major_raw_mean = major_raw_mean.T
major_raw_mean.columns = major_raw_mean.iloc[0]
major_raw_mean["Marker"] = major_raw_mean.index.astype(str).str.capitalize()
major_raw_mean = major_raw_mean.drop(major_raw_mean.index[0])
merged_filtered_markers = merged_filtered_markers.merge(major_raw_mean, on = 'Marker', how = 'left')

## Expression Filtering
merged_filtered_markers["Detectable_Expression"] = merged_filtered_markers["Major_Detectable_Expression"] | merged_filtered_markers["Minor_Detectable_Expression"]
merged_filtered_markers["not_crowding"] = np.sum(merged_filtered_markers.loc[:, 'AC':'Rod'] > 100, axis = 1) == 0

# Filter
merged_filtered_markers["xenium_filter"] = (merged_filtered_markers["long_enough"] & 
                                            merged_filtered_markers["Detectable_Expression"] &
                                            merged_filtered_markers["not_crowding"] &
                                            ~merged_filtered_markers["small_coef_in_query"])
merged_filtered_markers = merged_filtered_markers.sort_values(['Major_Name', 'Queried_Name']) # For now
merged_filtered_markers.to_csv('05_Filter_Merged_Markers/5_merged_markers_xeniumAnnotated.txt', sep = '\t')
merged_filtered_markers = merged_filtered_markers.loc[filtered_indices]
merged_filtered_markers.to_csv('05_Filter_Merged_Markers/5_merged_markers_xeniumFiltered.txt', sep = '\t')

# 05.2 # This is just looking at curated, xenium-filtered markers
# merged_filtered_markers = pd.read_csv('05_Filter_Merged_Markers/5_curated_markers_annotatedKeep_lengthExpressionMissingFiltered.txt', sep ='\t')
# merged_filtered_markers = merged_filtered_markers.sort_values(['Queried_Major_Name', 'Queried_Name']) # For now

def clean_markers(var_names, dirty_markers, verbose=True):
    final_markers = []
    for m in dirty_markers:
        if m in final_markers: # sc.pl.dotplot throws a fit if there are duplicates
            print(f'Found duplicate marker {m}!')
            continue
        if m not in var_names:
            print(f'Marker {m} from curate is not in this data!')
            continue
        final_markers += [m]
    if verbose:
        print(f'Final Markers: {final_markers}')
    return(final_markers)

# 05.2
if plot_major:
    assume_correct = merged_filtered_markers["Queried_Major_Name"] == merged_filtered_markers["Queried_Name"]
    hedge_unharmonized = merged_filtered_markers["Queried_Major_Name"] == merged_filtered_markers["Name"]
    curated_majorclass_markers = merged_filtered_markers.loc[assume_correct | hedge_unharmonized]
    curated_majorclass_markers.to_csv('05_Filter_Merged_Markers/5_curatedMarkers_majorclass_xeniumFiltered.txt', sep = '\t', index = False)
    
    final_markers = clean_markers(adata.var_names, curated_majorclass_markers["Marker"].tolist())
    
    sc.pl.dotplot(adata[:, final_markers],
                  var_names = final_markers,
                  groupby = 'majorclass',
                  categories_order = adata.obs["majorclass"].astype(str).drop_duplicates(keep='first').sort_values().tolist(), # Only celltypes that have a marker should be present
                  vmax = max_col,
                  vmin = 0,
                  show = False,
                  save = f"/project/hipaa_ycheng11lab/atlas/CAMR2024/05_Filter_Merged_Markers/figures/5_dotplot_mouseRetina_majorclass_curatedMarkers_{data_string}.pdf")


for majorclass in adata.obs['majorclass'].cat.categories:
    
    print(f'{datetime.datetime.now()} Major Class: {majorclass}')
    
    majorclass_original = majorclass
    majorclass = majorclass.upper()
    
    is_major_marker = np.logical_and(merged_filtered_markers["Major_Name"] == majorclass, merged_filtered_markers["Name"] == majorclass)
    cell_markers = merged_filtered_markers.loc[is_major_marker, "Marker"].tolist()
    
    is_subtype_marker = np.logical_and(merged_filtered_markers["Major_Name"] == majorclass, merged_filtered_markers["Name"] != majorclass)
    subtype_markers = merged_filtered_markers.loc[is_subtype_marker, "Marker"].tolist()
    
    markers = cell_markers + subtype_markers

    final_markers = clean_markers(adata.var_names, markers)

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

# Make majorclass in author_cell_type uppercase to match Name and Major_Name
majorclass_idx = adata.obs["author_cell_type"].str.upper().isin(adata.obs["majorclass"].astype(str).str.upper())
adata.obs.loc[majorclass_idx, "author_cell_type"] = adata.obs.loc[majorclass_idx, "author_cell_type"].str.upper()

merged_filtered_markers_Curated = merged_filtered_markers.loc[merged_filtered_markers["Curated"] == "Curated"]
merged_filtered_markers_AC  = merged_filtered_markers_Curated.loc[merged_filtered_markers_Curated["Major_Name"] == "AC"]
merged_filtered_markers_BC  = merged_filtered_markers_Curated.loc[merged_filtered_markers_Curated["Major_Name"] == "BC"]
merged_filtered_markers_RGC = merged_filtered_markers_Curated.loc[merged_filtered_markers_Curated["Major_Name"] == "RGC"]
merged_filtered_markers_MG  = merged_filtered_markers_Curated.loc[merged_filtered_markers_Curated["Major_Name"] == "MG"]
merged_filtered_markers_RPE = merged_filtered_markers_Curated.loc[merged_filtered_markers_Curated["Major_Name"] == "RPE"]
merged_filtered_markers_Curated = pd.concat([merged_filtered_markers_AC,
                                             merged_filtered_markers_BC,
                                             merged_filtered_markers_RGC,
                                             merged_filtered_markers_MG,
                                             merged_filtered_markers_RPE],
                                            ignore_index = True)

ordered_celltypes = merged_filtered_markers_Curated["Queried_Name"].drop_duplicates().astype('category').cat.remove_categories('RGC').dropna().tolist() + ["ROD","CONE","HC","MICROGLIA","ENDOTHELIAL"]

august_grant_markers = ["Ptn","Cntn6","Nxph1","Cpne4","Cbln4","Etv1","Epha3","Trhde", # AC
                        "Neto1","Erbb4","Grik1","Nnat","Cabp5","Sox6","Cck","Cpne9","Prkca","Hcn1", # BC
                        "Rbpms","Penk","Spp1","Coch","Opn4","Il1rapl2","Tbx20","Tac1","Cartpt", # RGC
                        "Aqp4","Glul","Rlbp1","Slc1a3","Vim", # MG
                        "Rpe65", # RPE
                        "Rho", "Arr3","Lhx1","Onecut1", "Cd74", "Cldn5"] # Manual

# Temporary change for this plot
adata.obs.loc[adata.obs["author_cell_type"].astype(str).str.contains("Microglia"), "author_cell_type"] = "MICROGLIA" # No subtypes of Microglia here

sc.pl.dotplot(adata[adata.obs['author_cell_type'].isin(ordered_celltypes), :],
              var_names = august_grant_markers,
              groupby = 'author_cell_type',
              categories_order = ordered_celltypes,
              vmax = max_col,
              vmin = 0,
              show = False,
              figsize = (12,10),
              save = f"mouseRetina_filteredMergedMarkers_curated_{data_string}.pdf")
