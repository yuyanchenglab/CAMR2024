#!/usr/bin/env python
# coding: utf-8

# Marker Determination

import sklearn as sk
import anndata as ad
import scanpy as sc 
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import joblib
import datetime as dt

sc.settings.n_jobs = -1

os.chdir("/project/hipaa_ycheng11lab/atlas/CAMR2024/")

is_verbose = True

count_lowcluster = 4 # Recommended detection limit for cell markers 
count_highcluster = 100 # Recommended detection ceiling

sc.plotting.DotPlot.DEFAULT_SAVE_PREFIX = "" # "03_Filter_Model_Markers/figures/3_dotplot_"
sc.plotting.DotPlot.DEFAULT_LARGEST_DOT = 200.0

raw = True
data_string = "normCounts" # 13 GB
if raw:
    data_string = "rawCounts"

final_majorclass_candidates_ordered = pd.read_csv('03_Filter_Model_Markers/3_ovr_LogReg_majorclass_xeniumFiltered.txt', sep = '\t') # No need to recalculate for majorclass, variable genes should be good enough
subtype_to_type = pd.read_csv('02_Modeling/2_minorToMajorClass.txt', sep = '\t')

adata = ad.read_h5ad('01_QualityControl/1_camr_scrublet_batch_filtered.h5ad')
adata.var.index = adata.var["feature_name"] # subset on genes instead of booleans for dotplots
if raw:
    adata.X = adata.raw.X
raw_mean_expression_minorclass = pd.read_csv(f'data/raw_meanExpression_minorclass.txt', sep = '\t', index_col=0)

## Subtype Markers

def get_top_coefficient_genes(majorclass: str):
    top_features_log_reg_sub = pd.read_csv(f'02_Modeling/minorclass/2_ovr_LogReg_minorclass-{majorclass}_AbsTop20Markers.txt', sep ='\t')
    if majorclass != 'Microglia':
        top_features_log_reg_sub = top_features_log_reg_sub[top_features_log_reg_sub['Coefficient'] > 0]
    # top_features_log_reg_sub['Coefficient'] = abs(top_features_log_reg_sub['Coefficient'])
    top_features_log_reg_sub.index = top_features_log_reg_sub.Marker
    
    return top_features_log_reg_sub


# Filter based on innate features of the gene itself
def filter_gene_by_genomics(adata, top_coefficient_genes, length_threshold = 960, verbose = False):
    in_regression = adata.var["feature_name"].astype(str).isin(top_coefficient_genes["Marker"])
    long_enough = adata.var["feature_length"].astype(int) >= length_threshold # It's a conservative filter
    
    is_genomics_candidate = long_enough & in_regression
    genomics_candidates = adata.var["feature_name"].astype(str)[is_genomics_candidate.tolist()].tolist()
    
    if verbose:
        print(len(genomics_candidates), genomics_candidates)
    
    return genomics_candidates


# Filter based on the filtering criteria even though adata.obs.library_platform.unique() is a mix of 4 chemistries...
def filter_gene_by_expression(adata, raw_mean_expression, count_lowcluster = 4, count_highcluster = 100, verbose = False):
    detectable_genes = (raw_mean_expression >= count_lowcluster).sum(axis=0) >= 1 # Leaky, expression may be highest in off-target cell
    optical_crowding_genes = (raw_mean_expression > count_highcluster).sum(axis=0) > 0
    
    is_expression_candidate = detectable_genes & (~optical_crowding_genes)
    expression_candidates = adata.raw.var["feature_name"].astype(str)[is_expression_candidate.tolist()].tolist()
    
    if verbose:
        print(len(expression_candidates), expression_candidates)
    
    return expression_candidates


def merge_major_minor_markers(majorclass_candidates, minorclass_candidates, majorclass, subtype_to_type):
    cell_markers = majorclass_candidates.loc[majorclass_candidates["Name"] == majorclass, "Marker"]
    subtypes = subtype_to_type.loc[subtype_to_type.majorclass == majorclass, "minorclass"].tolist()
    subtype_markers = minorclass_candidates[minorclass_candidates.isin(subtypes)].index
    
    markers = cell_markers.tolist() + subtype_markers.tolist()
    
    unique_markers = [] # sc.pl.dotplot throws a fit if there are duplicates
    for m in markers:
        if m not in unique_markers:
            unique_markers += [m]
    
    return unique_markers

        
for majorclass in ["AC", "BC", "Microglia", "RGC"]: # adata.obs['majorclass'].cat.categories, but only ones with subtypes to check

    print(f'{dt.datetime.now()} Major Class: {majorclass}')

    top_coefficient_genes = get_top_coefficient_genes(majorclass)
    top_coefficient_genes = top_coefficient_genes.loc[top_coefficient_genes["Marker"].value_counts()==1] # Unique filter
    raw_mean_expression = raw_mean_expression_minorclass.loc[subtype_to_type.loc[subtype_to_type["majorclass"].astype(str) == majorclass, "minorclass"]]

    genomics_candidates = filter_gene_by_genomics(adata, top_coefficient_genes, verbose = is_verbose)
    expression_candidates = filter_gene_by_expression(adata, raw_mean_expression, verbose = is_verbose)
    final_candidates = np.intersect1d(genomics_candidates, expression_candidates)

    if is_verbose:
        print(len(final_candidates), final_candidates)

    ordered_markers = top_coefficient_genes.loc[final_candidates.tolist()].sort_values(['Major_Name', 'Name'])
    ordered_markers.to_csv(f'03_Filter_Model_Markers/3_ovr_LogReg_minorclass-{majorclass}_xeniumFiltered.txt')
    
    major_minor_markers = merge_major_minor_markers(final_majorclass_candidates_ordered, ordered_markers, majorclass, subtype_to_type)
    
    sc.pl.dotplot(
        adata = adata[adata.obs.majorclass == majorclass, major_minor_markers],
        var_names = major_minor_markers,
        gene_symbols = "feature_name",
        groupby = 'author_cell_type',
        vmax = count_lowcluster * 3,
        vmin = count_lowcluster - 1,
        show = False,
        save = f"/project/hipaa_ycheng11lab/atlas/CAMR2024/03_Filter_Model_Markers/figures/3_dotplot_mouseRetina_minorclass-{majorclass}_xeniumFiltered." +
               f"{count_lowcluster}-{count_highcluster}_{data_string}.pdf")
