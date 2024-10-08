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

os.chdir('/project/ycheng11lab/jfmaurer/mouse_retina_atlas_chen_2024/')
os.makedirs('14_Plot_V4', exist_ok = True)
sc.settings.n_jobs = -1

adata = ad.read_h5ad('10_Make_Shiny/10_Shiny_Input.h5ad')

resultsPath = "09_Designer_Analysis/PanelDesignV4.txt"
markers = pd.read_csv(resultsPath, sep = '\t').sort_values(["Major_Name", "Name"]) # NOTE: Nrg1 & 2010007h06rik adjusted
markers = pd.concat([markers.loc[markers["Name"] == markers["Major_Name"], :], markers.loc[markers["Name"] != markers["Major_Name"], :]])
all_markers = markers["Marker"].unique().tolist()
major_markers = markers.loc[markers["Name"] == markers["Major_Name"], "Marker"].unique().tolist()

sc.plotting.DotPlot.DEFAULT_SAVE_PREFIX = ""
sc.plotting.DotPlot.DEFAULT_LARGEST_DOT = 200.0
plot_prefix = "/project/ycheng11lab/jfmaurer/mouse_retina_atlas_chen_2024/14_Plot_V4/V4_"

raw = True
data_string = "normCounts"
max_col = 3

for raw in [False, True]:
  
  if raw:
    data_string = "rawCounts"
    max_col = 12
    adata.X = adata.raw.X
  
  sc.pl.dotplot(adata[:, all_markers],
              var_names = all_markers, gene_symbols = "feature_name",
              groupby = "minorclass",
              vmax = max_col,
              vmin = 0,
              show = False,
              save = f"{plot_prefix}Minor-Cell_All-Marker_{data_string}.pdf")
  
  for cellType in ['AC', 'BC','RGC','Microglia']:
      sc.pl.dotplot(adata[adata.obs["majorclass"] == cellType, all_markers],
              var_names = all_markers, gene_symbols = "feature_name",
              groupby = "minorclass",
              vmax = max_col,
              vmin = 0,
              show = False,
              save = f"{plot_prefix}{cellType}_All-Marker_{data_string}.pdf")

  sc.pl.dotplot(adata[:, all_markers],
              var_names = all_markers, gene_symbols = "feature_name",
              groupby = "majorclass",
              vmax = max_col,
              vmin = 0,
              show = False,
              save = f"{plot_prefix}Major-Cell_All-Marker_{data_string}.pdf")
              
  sc.pl.dotplot(adata[:, major_markers],
              var_names = major_markers, gene_symbols = "feature_name",
              groupby = "majorclass",
              vmax = max_col,
              vmin = 0,
              show = False,
              save = f"{plot_prefix}Major-Cell_Major-Marker_{data_string}.pdf")
