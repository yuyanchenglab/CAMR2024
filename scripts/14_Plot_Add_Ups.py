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
os.makedirs('14_Add_Ups', exist_ok = True)
sc.settings.n_jobs = -1

adata = ad.read_h5ad('10_Make_Shiny/10_Shiny_Input.h5ad')


all_markers = [
'Tmcc3',
'Anks1b',
'Asic2',
'Vsx2',
'Rho',
'Pde6a',
'Nr2e3',
'Cnga1',
'Guca1b',
'Rp1',
'Rcvrn',
# Second set below
'Pax6',
'Tfap2a',
'Gad2',
'Slc6a1',
'Slc32a1',
'Slc6a9',
'Slc17a6',
'Pou4f2',
'Rbpms',
'Onecut1',
'Onecut2']

sc.plotting.DotPlot.DEFAULT_SAVE_PREFIX = ""
sc.plotting.DotPlot.DEFAULT_LARGEST_DOT = 200.0
plot_prefix = "/project/ycheng11lab/jfmaurer/mouse_retina_atlas_chen_2024/14_Add_Ups/AddUp_"

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

  sc.pl.dotplot(adata[:, all_markers],
              var_names = all_markers, gene_symbols = "feature_name",
              groupby = "majorclass",
              vmax = max_col,
              vmin = 0,
              show = False,
              save = f"{plot_prefix}Major-Cell_All-Marker_{data_string}.pdf")
