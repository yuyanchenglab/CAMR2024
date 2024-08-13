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

majorclass_candidates = pd.read_csv('spreadsheets/ovr_top_filtered_genes_majorclass_coefficients_sensitive.csv', index_col = 0)
subtype_to_type = pd.read_csv('spreadsheets/major_author.csv', index_col = 0)

adata_full = ad.read_h5ad('data/camr_scrublet_batch_filtered.h5ad')
adata_full.var.index = adata_full.var["feature_name"] # subset on genes instead of booleans for dotplots

## Subtype Markers

def get_top_coefficient_genes(gene_filter: str, majorclass: str):
    top_coefficient_path = f'spreadsheets/ovr_top_20_{gene_filter}_genes_{majorclass}_author_cell_type.csv'
    top_features_log_reg_sub = pd.read_csv(top_coefficient_path)
    if majorclass != 'Microglia':
        top_features_log_reg_pos_sub = top_features_log_reg_sub[top_features_log_reg_sub['Coefficient'] > 0]
    else:
        top_features_log_reg_pos_sub = top_features_log_reg_sub
    top_features_log_reg_pos_sub.index = top_features_log_reg_pos_sub.Gene
    
    return top_features_log_reg_pos_sub


def get_raw_expression(adata, majorclass: str):
    if not os.path.isfile(f'data/raw_{majorclass}_subtype_meanExpression.csv'):
        raw_feature_expression = pd.DataFrame(adata.raw.X.toarray(), columns = adata.raw.var["feature_name"].astype(str).tolist())
        raw_feature_expression_sub = raw_feature_expression.loc[(adata.obs['majorclass'] == majorclass).tolist()]
        raw_feature_expression_sub["author_cell_type"] = adata.obs.loc[(adata.obs['majorclass'] == majorclass).tolist(), "author_cell_type"].tolist()
        raw_feature_expression_sub_mean = raw_feature_expression_sub.groupby("author_cell_type").agg("mean")
        raw_feature_expression_sub_sum = raw_feature_expression_sub.groupby("author_cell_type").agg("sum")
        del raw_feature_expression_sub
        raw_feature_expression_sub_mean.to_csv(f'data/raw_{majorclass}_subtype_meanExpression.csv')
        raw_feature_expression_sub_sum.to_csv(f'data/raw_{majorclass}_subtype_sumExpression.csv')
    else:
        raw_feature_expression_sub_mean = pd.read_csv(f'data/raw_{majorclass}_subtype_meanExpression.csv', index_col=0)
        raw_feature_expression_sub_sum = pd.read_csv(f'data/raw_{majorclass}_subtype_sumExpression.csv', index_col=0)
    
    return (raw_feature_expression_sub_mean, raw_feature_expression_sub_sum)


# Filter based on innate features of the gene itself
def filter_gene_by_genomics(adata, top_coefficient_genes, length_threshold = 960, verbose = False):
    in_regression = adata.var["feature_name"].astype(str).isin(top_coefficient_genes["Gene"])
    long_enough = adata.var["feature_length"].astype(int) >= length_threshold # It's a conservative filter
    
    is_genomics_candidate = long_enough & in_regression
    genomics_candidates = adata.var["feature_name"].astype(str)[is_genomics_candidate.tolist()].tolist()
    
    if verbose:
        print(len(genomics_candidates), genomics_candidates)
    
    return genomics_candidates


# Filter based on the filtering criteria even though adata.obs.library_platform.unique() is a mix of 4 chemistries...
def filter_gene_by_expression(adata, raw_mean_expression, count_lowcluster = 4, count_highcluster = 100, verbose = False):
    detectable_genes = (raw_mean_expression >= count_lowcluster).sum(axis=0) >= 1
    optical_crowding_genes = (raw_mean_expression > count_highcluster).sum(axis=0) > 0
    
    is_expression_candidate = detectable_genes & (~optical_crowding_genes)
    expression_candidates = adata.raw.var["feature_name"].astype(str)[is_expression_candidate.tolist()].tolist()
    
    if verbose:
        print(len(expression_candidates), expression_candidates)
    
    return expression_candidates


def merge_major_minor_markers(majorclass_candidates, minorclass_candidates, majorclass, subtype_to_type):
    cell_markers = majorclass_candidates[majorclass_candidates == majorclass].index
    subtypes = subtype_to_type.loc[subtype_to_type.majorclass == majorclass, "author_cell_type"].tolist()
    subtype_markers = minorclass_candidates[minorclass_candidates.isin(subtypes)].index
    
    markers = cell_markers.tolist() + subtype_markers.tolist()
    
    unique_markers = [] # sc.pl.dotplot throws a fit if there are duplicates
    for m in markers:
        if m not in unique_markers:
            unique_markers += [m]
    
    return unique_markers


for gene_filter in ['highly_variable', 'moderately_variable', 'complete']:

    print(f'{dt.datetime.now()} Filter Strategy: {gene_filter}')

    if gene_filter == 'highly_variable':
        adata = adata_full[:, adata_full.var.highly_variable]
    if gene_filter == 'moderately_variable':
        adata = adata_full[:, adata_full.var['dispersions'] > min(adata_full.var.loc[adata_full.var['highly_variable'], 'dispersions'])]
    if gene_filter == 'complete':
        adata = adata_full
        
    for majorclass in ["AC", "BC", "Microglia", "RGC"]: # adata.obs['majorclass'].cat.categories, but only ones with subtypes to check
    
        print(f'{dt.datetime.now()} Major Class: {majorclass}')

        top_coefficient_genes = get_top_coefficient_genes(gene_filter, majorclass)
        raw_mean_expression, raw_sum_expression = get_raw_expression(adata, majorclass)

        genomics_candidates = filter_gene_by_genomics(adata, top_coefficient_genes, verbose = is_verbose)
        expression_candidates = filter_gene_by_expression(adata, raw_mean_expression, verbose = is_verbose)
        final_candidates = np.intersect1d(genomics_candidates, expression_candidates)

        if is_verbose:
            print(len(final_candidates), final_candidates)

        ordered_markers = top_coefficient_genes.loc[final_candidates.tolist()].sort_values('Subclass')
        ordered_markers.to_csv(f'spreadsheets/ovr_top_{gene_filter}_filtered_markers_{majorclass}_coefficients_sensitive.csv')
        
        # major_minor_markers = merge_major_minor_markers(majorclass_candidates, ordered_markers, majorclass, subtype_to_type)

        # sc.pl.dotplot(
        #     adata[adata.obs.majorclass == majorclass, major_minor_markers],
        #     major_minor_markers,
        #     gene_symbols = "feature_name",
        #     groupby = 'author_cell_type',
        #     vmax = count_lowcluster * 3,
        #     vmin = count_lowcluster - 1,
        #     show = False,
        #     save = f"mouseRetina_{gene_filter}_{majorclass}_byClass_filteredCounts." +
        #            f"{count_lowcluster}-{count_highcluster}.pdf")
