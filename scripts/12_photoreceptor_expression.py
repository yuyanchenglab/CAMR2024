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

os.chdir('/project/ycheng11lab/jfmaurer/mouse_retina_atlas_chen_2024')
os.makedirs('12_photoreceptor_expression', exist_ok = True)
sc.settings.n_jobs = -1

adata = ad.read_h5ad('10_Make_Shiny/10_Shiny_Input.h5ad')


rod_gene = ['Rho','Pde6a', 'Gngt1', 'Optn','Nrl','Nr2e3', 'Reep6', 'Cnga1', 'Guca1b', 'Pde6a', 'Rp1',
              'Cngb1','Rcvrn', 'Pdc','Syne2','Mef2c','Fyco1','Atf4']
cone_gene = ['Opn1mw','Opn1sw', 'Ccdc136', 'Optn','Nrl','Gngt2', 'Gnat2']

sc.plotting.DotPlot.DEFAULT_SAVE_PREFIX = ""
sc.plotting.DotPlot.DEFAULT_LARGEST_DOT = 200.0
plot_prefix = "/project/ycheng11lab/jfmaurer/mouse_retina_atlas_chen_2024/12_photoreceptor_expression/"

raw = True
data_string = "normCounts"
max_col = 3

for raw in [False, True]:
  
  if raw:
    data_string = "rawCounts"
    max_col = 12
    adata.X = adata.raw.X

  # rod_gene = set(rod_gene)
  # cone_gene = set(cone_gene)

  sc.pl.dotplot(adata[:, rod_gene],
              var_names = rod_gene, gene_symbols = "feature_name",
              groupby = "majorclass",
              vmax = max_col,
              vmin = 0,
              show = False,
              save = f"{plot_prefix}Major-Cell_Rod_{data_string}.pdf")

  sc.pl.dotplot(adata[:, cone_gene],
              var_names = cone_gene, gene_symbols = "feature_name",
              groupby = "majorclass",
              vmax = max_col,
              vmin = 0,
              show = False,
              save = f"{plot_prefix}Major-Cell_Cone_{data_string}.pdf")
