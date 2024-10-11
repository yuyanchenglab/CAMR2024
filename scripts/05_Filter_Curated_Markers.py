#!/usr/bin/env python3
# coding: utf-8

# JM 08/17/24
# In this notebook, we're going plot the mouse retinal data and determine proper marker genes.
# Also apply xenium filters to the curated data

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
os.makedirs('05_Filter_Curated_Markers', exist_ok = True)
sc.settings.n_jobs = -1

adata = ad.read_h5ad('01_QualityControl/1_camr_scrublet_batch_filtered.h5ad')

sc.plotting.DotPlot.DEFAULT_SAVE_PREFIX = "05_Filter_Curated_Markers/figures/5_dotplot_"
sc.plotting.DotPlot.DEFAULT_LARGEST_DOT = 200.0

coefficient_threshold = 0.3
length_threshold = 960
raw_expression_threshold = 4

plot_major_markers = False
plot_major_cells = True
plot_minor_markers = True
plot_minor_cells = True
plot_split_majorclass = True
plot_only_target_cells = False
plot_only_curated = False

raw = False

data_string = "normCounts"
max_col = 3
if raw:
    data_string = "rawCounts"
    max_col = 4
    adata.X = adata.raw.X

# plot_occassion = "" # Options: "august_grant": # 05.1, "curated_xenium_filtered": # 05.2
# data_string = data_string + f"_{plot_occassion}"

# Clean subtypes
adata.obs["author_cell_type"] = adata.obs["author_cell_type"].astype(str).str.upper() # From category
adata.obs["majorclass"] = adata.obs["majorclass"].astype(str).str.upper() # From category
is_unassigned = adata.obs["author_cell_type"] == adata.obs["majorclass"]
is_subtype = adata.obs["author_cell_type"].isin(["AC", "BC", "MICROGLIA", "RGC"])
unassigned_subtypes = adata.obs["author_cell_type"].loc[is_unassigned & is_subtype]
adata.obs.loc[is_unassigned & is_subtype, "author_cell_type"] = ["UNASSIGNED_" + uv for uv in unassigned_subtypes]
adata.obs["minorclass"] = adata.obs["author_cell_type"]
adata.obs["Name"] = adata.obs["author_cell_type"] # This should be fed through q2n
adata.obs["Major_Name"] = adata.obs["majorclass"].astype(str) # This should be fed through q2n

# Shrink adata footprint
adata.raw = None # 5GB saved?
adata.var["Ensembl"] = adata.var.index.tolist()
adata.var["feature_name"] = adata.var["feature_name"].astype(str).str.capitalize()
adata.var["feature_length"] = adata.var["feature_length"].astype(int)
adata.var_names = adata.var["feature_name"].astype(str).tolist() # No need to manually adjust index
adata.var_names_make_unique() # Unnecessary for this data


# Work starts here

if os.path.isfile('05_Filter_Curated_Markers/5_curated_markers_xeniumFiltered.txt'):
    filtered_markers = pd.read_csv('05_Filter_Curated_Markers/5_curated_markers_xeniumFiltered.txt', sep = '\t')
else:
    # Coefficient Marking # Jumping the shark, don't merge until the wet lab does it.
    curated_markers = pd.read_csv('04_Harmonize_Curated_Markers/4_harmonized_curated_markers.txt', sep = '\t').drop_duplicates()
    
    # Add gene length
    curated_markers["Marker"] = curated_markers["Marker"].str.capitalize() # Make sure the below join keeps all rows
    curated_markers.index = curated_markers["Marker"] # Make sure the below join keeps all rows
    curated_markers = curated_markers.join(adata.var[["feature_length", "Ensembl"]])
    curated_markers["missing_length"] = curated_markers["feature_length"].isnull()
    curated_markers[curated_markers["missing_length"], "feature_length"] = 0
    curated_markers["Long_Enough"] = curated_markers["feature_length"] >= length_threshold
    curated_markers.index = list(range(curated_markers.shape[0])) # Move this to after minorclass merges?
    
    # minorclass
    ## Pre-process expression table names
    subcell_raw_mean = pd.read_csv('data/raw_meanExpression_minorclass.txt', sep = '\t')
    subcell_raw_mean["author_cell_type"].str.upper()
    is_subtype = subcell_raw_mean["author_cell_type"].isin(["AC", "BC", "MICROGLIA", "RGC"])
    unassigned_subtypes = subcell_raw_mean.loc[is_subtype, "author_cell_type"]
    subcell_raw_mean.loc[is_subtype, "author_cell_type"] = ["UNASSIGNED_" + uv for uv in unassigned_subtypes]
    
    ## Process table and detection
    subcell_raw_mean_long = pd.melt(subcell_raw_mean, id_vars='author_cell_type', var_name='Marker', value_name='Raw_Mean_Expression_Minorclass_Target')
    subcell_raw_mean_long = subcell_raw_mean_long.rename(columns={'author_cell_type':'Queried_Name'})
    subcell_raw_mean_long["Marker"] = subcell_raw_mean_long["Marker"].str.capitalize()
    subcell_raw_mean_long["Detectable_Minor_Expression"] = subcell_raw_mean_long['Raw_Mean_Expression_Minorclass_Target'] >= raw_expression_threshold
    curated_markers = curated_markers.merge(subcell_raw_mean_long, on=['Queried_Name', 'Marker'], how = 'left')
    
    ## Optical crowding
    minorToMajor = adata.obs[["author_cell_type", "majorclass"]].drop_duplicates()
    subcell_raw_mean.index = subcell_raw_mean['author_cell_type']
    
    minor_crowding = np.zeros(sum(curated_markers["Marker"].isin(subcell_raw_mean.columns)))
    included_markers = curated_markers.loc[curated_markers["Marker"].isin(subcell_raw_mean.columns), "Marker"]
    for majorclass in ["AC", "BC", "MICROGLIA", "RGC"]:
        minorclasses = minorToMajor.loc[minorToMajor["majorclass"] == majorclass, "author_cell_type"] # Why is below line needed given the 'clean_genes()' from last script? How well does 'upper' fix this?
        sub_expr_mtx = subcell_raw_mean.loc[minorclasses, included_markers]
        minor_crowding = minor_not_crowding & (np.sum(sub_expr_mtx > 100, axis = 0) > 0)
    
    minor_crowding.name = "Minor_Crowding_Risk" # Is this not a numpy array?
    curated_markers.index = curated_markers["Marker"]
    curated_markers = curated_markers.join(minor_crowding, how = "left").drop_duplicates() # Drop duplicate necessary?
    curated_markers.loc[curated_markers["Minor_Crowding_Risk"].isnull(), "Minor_Crowding_Risk"] = True # Be alittle leaky
    
    # majorclass expression
    major_raw_mean = pd.read_csv('data/raw_meanExpression_majorclass.txt', sep = '\t')
    major_raw_mean_long = pd.melt(major_raw_mean, id_vars='majorclass', var_name='Marker', value_name='Raw_Mean_Expression_Majorclass_Target')
    major_raw_mean_long = major_raw_mean_long.rename(columns={'majorclass':'Major_Name'}) # Major_Queried name?
    major_raw_mean_long["Marker"] = major_raw_mean_long["Marker"].str.capitalize()
    major_raw_mean_long["Detectable_Major_Expression"] = major_raw_mean_long["Raw_Mean_Expression_Majorclass_Target"] >= raw_expression_threshold
    curated_markers = curated_markers.merge(major_raw_mean_long, on=['Major_Name', 'Marker'], how = 'left')
    
    ## all majorclass expression
    major_raw_mean = major_raw_mean.T
    major_raw_mean.columns = major_raw_mean.iloc[0] # majorclass in first column becomes first row
    major_raw_mean["Marker"] = major_raw_mean.index.astype(str).str.capitalize()
    major_raw_mean = major_raw_mean.drop(major_raw_mean.index[0])
    curated_markers = curated_markers.merge(major_raw_mean, on = 'Marker', how = 'left')
    curated_markers["Major_Crowding_Risk"] = np.sum(curated_markers.loc[:, adata.obs["majorclass"].unique()] > 100, axis = 1) > 0 # If subtype is high, that should be ok?
    
    # Combine Major and Minor
    curated_markers["Detectable_Expression"] = curated_markers["Detectable_Major_Expression"] | curated_markers["Detectable_Minor_Expression"]
    curated_markers["Optical_Crowding_Risk"] = curated_markers["Major_Crowding_Risk"] | curated_markers["Minor_Crowding_Risk"]
    
    # Filter
    curated_markers["Xenium_Filter"] = (curated_markers["Long_Enough"] &
                                        curated_markers["Detectable_Expression"] &
                                        ~curated_markers["Optical_Crowding_Risk"])
    curated_markers = curated_markers.sort_values(['Major_Name', 'Queried_Name']) # For now
    curated_markers.to_csv('05_Filter_Curated_Markers/5_curated_markers_xeniumAnnotated.txt', sep = '\t')
    filtered_markers = curated_markers.loc[curated_markers["Xenium_Filter"]]
    filtered_markers.to_csv('05_Filter_Curated_Markers/5_curated_markers_xeniumFiltered.txt', sep = '\t')

if plot_only_curated:
    data_string = data_string + '_CuratedOnly'
else: # Road less travelled
    queried_markers = pd.read_csv("04_Harmonize_Curated_Markers/minor_by_major/4_queried.txt", sep = '\t')
    filtered_markers = filtered_markers.join(queried_markers)

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
if plot_major_markers and plot_major_cells:
    assume_correct = filtered_markers["Queried_Major_Name"] == filtered_markers["Queried_Name"]
    hedge_unharmonized = filtered_markers["Queried_Major_Name"] == filtered_markers["Name"]
    curated_majorclass_markers = filtered_markers.loc[assume_correct | hedge_unharmonized]
    
    final_markers = clean_markers(adata.var_names, curated_majorclass_markers["Marker"].tolist())
    
    sc.pl.dotplot(adata[:, final_markers],
                  var_names = final_markers,
                  groupby = 'majorclass',
                  categories_order = adata.obs["majorclass"].astype(str).drop_duplicates(keep='first').sort_values().tolist(), # Only celltypes that have a marker should be present
                  vmax = max_col,
                  vmin = 0,
                  show = False,
                  save = f"/project/ycheng11lab/jfmaurer/mouse_retina_atlas_chen_2024/05_Filter_Curated_Markers/figures/5_dotplot_mouseRetina_majorclass_curatedMarkers_{data_string}.pdf")
# End plot_major_markers and plot_major_cells

if plot_minor_markers and plot_minor_cells:
    all_markers = []
    all_target_subtypes = []
    for majorclass in adata.obs['majorclass'].cat.categories:
        
        print(f'{datetime.datetime.now()} Major Class: {majorclass}')
        
        majorclass_original = majorclass
        majorclass = majorclass.upper()
        
        is_subtype_marker = np.logical_and(filtered_markers["Major_Name"] == majorclass, filtered_markers["Name"] != majorclass)
        markers = filtered_markers.loc[is_subtype_marker, "Marker"].tolist()
        target_subtypes = filtered_markers.loc[is_subtype_marker, "Queried_Name"].drop_duplicates(keep='first').tolist()
        
        if plot_major_markers:
            is_major_marker = np.logical_and(filtered_markers["Major_Name"] == majorclass, filtered_markers["Name"] == majorclass)
            cell_markers = filtered_markers.loc[is_major_marker, "Marker"].tolist()
            markers = cell_markers + markers
        
        final_markers = clean_markers(adata.var_names, markers)
        all_markers += final_markers
        
        
        if not plot_split_majorclass:
            continue
        
        if final_markers == []: # Necessary to avoid "ValueError: left cannot be >= right"
            print(f'No minor markers available for {majorclass}!')
            continue
        
        if not plot_only_target_cells:
            all_subtypes = adata.obs.loc[adata.obs['majorclass'] == majorclass_original, "author_cell_type"].astype(str).sort_values().drop_duplicates(keep='first')
            target_subtypes += all_subtypes[!all_subtypes.isin(target_subtypes)].tolist()
        all_target_subtypes += target_subtypes
        
        sc.pl.dotplot(adata[adata.obs['majorclass'] == majorclass_original, final_markers],
                      var_names = final_markers,
                      groupby = 'author_cell_type',
                      categories_order = target_subtypes,
                      vmax = max_col,
                      vmin = 0,
                      show = False,
                      save = f"mouseRetina_minorclass-{majorclass}_{data_string}.pdf")
    # End majorclass loop
    
    final_markers = clean_markers(adata.var_names, all_markers)
    sc.pl.dotplot(adata[:, final_markers],
                      var_names = final_markers,
                      groupby = 'author_cell_type',
                      categories_order = all_target_subtypes,
                      vmax = max_col,
                      vmin = 0,
                      show = False,
                      save = f"mouseRetina_minorclass_{data_string}.pdf")
# End plot_minor_markers and plot_minor_cells

if plot_occassion == "august_grant": # 05.1
    
    filtered_curated_markers = merged_filtered_markers.loc[merged_filtered_markers["Curated"] == "Curated"]
    filtered_markers_AC  = filtered_curated_markers.loc[filtered_curated_markers["Major_Name"] == "AC"]
    filtered_markers_BC  = filtered_curated_markers.loc[filtered_curated_markers["Major_Name"] == "BC"]
    filtered_markers_RGC = filtered_curated_markers.loc[filtered_curated_markers["Major_Name"] == "RGC"]
    filtered_markers_MG  = filtered_curated_markers.loc[filtered_curated_markers["Major_Name"] == "MG"]
    filtered_markers_RPE = filtered_curated_markers.loc[filtered_curated_markers["Major_Name"] == "RPE"]
    filtered_curated_markers = pd.concat([filtered_markers_AC,
                                         filtered_markers_BC,
                                         filtered_markers_RGC,
                                         filtered_markers_MG,
                                         filtered_markers_RPE],
                                        ignore_index = True)
    
    ordered_celltypes = filtered_curated_markers["Queried_Name"].drop_duplicates().astype('category').cat.remove_categories('RGC').dropna().tolist() + ["ROD","CONE","HC","MICROGLIA","ENDOTHELIAL"]
    
    august_grant_markers = ["Ptn","Cntn6","Nxph1","Cpne4","Cbln4","Etv1","Epha3","Trhde", # AC
                            "Neto1","Erbb4","Grik1","Nnat","Cabp5","Sox6","Cck","Cpne9","Prkca","Hcn1", # BC
                            "Rbpms","Penk","Spp1","Coch","Opn4","Il1rapl2","Tbx20","Tac1","Cartpt", # RGC
                            "Aqp4","Glul","Rlbp1","Slc1a3","Vim", # MG
                            "Rpe65", # RPE
                            "Rho", "Arr3","Lhx1","Onecut1", "Cd74", "Cldn5"] # Manual
    
    sc.pl.dotplot(adata[adata.obs['author_cell_type'].isin(ordered_celltypes), :],
                  var_names = august_grant_markers,
                  groupby = 'author_cell_type',
                  categories_order = ordered_celltypes,
                  vmax = max_col,
                  vmin = 0,
                  show = False,
                  figsize = (12,10),
                  save = f"mouseRetina_curated_{data_string}.pdf")
# End August Grant
