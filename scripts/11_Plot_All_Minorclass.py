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
os.makedirs('11_Plot_All_Minorclass', exist_ok = True)
sc.settings.n_jobs = -1

adata = ad.read_h5ad('10_Make_Shiny/10_Shiny_Input.h5ad')

resultsPath = "09_Designer_Analysis/PanelOutputV2.txt"
markers = pd.read_csv(resultsPath, sep = '\t') # NOTE: Nrg1 & 2010007h06rik adjusted

sc.plotting.DotPlot.DEFAULT_SAVE_PREFIX = ""
sc.plotting.DotPlot.DEFAULT_LARGEST_DOT = 200.0
plot_prefix = "/project/hipaa_ycheng11lab/atlas/CAMR2024/11_Plot_All_Minorclass/"

raw = True
data_string = "normCounts"
max_col = 3

for raw in [False, True]:
  
  if raw:
    data_string = "rawCounts"
    max_col = 12
    adata.X = adata.raw.X

  all_markers = markers["Marker"].unique().tolist()
  
  sc.pl.dotplot(adata[:, all_markers],
              var_names = all_markers, gene_symbols = "feature_name",
              groupby = "minorclass",
              vmax = max_col,
              vmin = 0,
              show = False,
              save = f"{plot_prefix}All-Cell_All-Marker_{data_string}.pdf")

exit()

group_var = "majorclass"
for raw in [False, True]:
  if raw:
    data_string = "rawCounts"
    max_col = 12
    adata.X = adata.raw.X

  for cell_set in ["majorclass", "AC", "BC", "Microglia", "RGC"]:
    
    # subset adata to cells
    if cell_set == "majorclass":
      adata_set = adata
      group_var = "majorclass"
    else:
      adata_set = adata[adata.obs['Major_Name'] == cell_set]
      group_var = "minorclass"
    cat_order = adata_set.obs[group_var].astype(str).drop_duplicates(keep='first').sort_values()#.tolist(), # Only celltypes that have a marker should be present
    
    for marker_set in ["majorclass", "AC", "BC", "Microglia", "RGC"]:
    
      # subset markers
      if marker_set == "majorclass":
        marker_set_values = markers.loc[markers["Name"] == markers["Major_Name"], "Marker"]
      else:
        marker_set_values = markers.loc[marker_set == markers["Major_Name"], "Marker"]

      unique_markers = []
      for m in marker_set_values: # sc.pl.dotplot throws a fit if there are duplicates
          if m not in unique_markers:
              unique_markers += [m]
          else:
              print(f'Found duplicate marker {m}!')
      
      final_markers = []
      for m in unique_markers: # Necessary to avoid "KeyError
          if m in adata.var_names:
              final_markers += [m]
          else:
              print(f'Marker {m} from curate is not in this data!')
      
      if final_markers == []: # Necessary to avoid "ValueError: left cannot be >= right"
        print(f'No {marker_set} markers available for {cell_set}!')
        continue
      
      sc.pl.dotplot(adata_set[:, final_markers],
                    var_names = final_markers, gene_symbols = "feature_name",
                    groupby = group_var,
                    categories_order = cat_order,
                    vmax = max_col,
                    vmin = 0,
                    show = False,
                    save = f"{plot_prefix}Marker-{marker_set}_Cell-{cell_set}_{data_string}.pdf")
