#!/bin/env Rscript3

# JM 08/12/2024
# How do my statistics methods compare to a literature search when looking for proper cell markers?

# Setup ----

library(dplyr) # pipe
library(data.table) # fast table manip
library(ggplot2)
library(stringr)
library(ggvenn) # venn
library(ggVennDiagram) # upset plot
library(eulerr)
library(magrittr)

analysisName = "4_Merge_Curated_Markers"
analysisPath = "/project/hipaa_ycheng11lab/atlas/CAMR2024/"
setwd(analysisPath)
source("scripts/4_Merge_Curated_Markers.util.R")

curatedPath = "data/CuratedMouseRetinaMarkers.txt"
queriedMajorPath = "03_Marker_Expression/3_ovr_LogReg_majorclass_xeniumFiltered.txt"

outPath = paste0(analysisPath, "04_Merge_Curated_Markers/")
dir.create(outPath, showWarnings = FALSE)

verbose = TRUE

## Reproducible Results ----

now = \() Sys.time() %>% format('%F_%H-%M-%S-%Z')
sessionInfoPath = paste0("logs/", analysisName, ".", now(), ".sessionInfo")
sessionInfo() %>% capture.output(file = sessionInfoPath)

set.seed(1)

# Data ----

majorclass = c("AC", "ASTROCYTE", "BC", "CONE", "ENDOTHELIAL", "HC", "MG", "MICROGLIA", "PERICYTE", "RGC", "ROD", "RPE")

 ]majorclass_with_minor = c("AC", "BC", "Microglia", "RGC")

curated <- fread(curatedPath)
queriedMajor <- fread(queriedMajorPath)
queriedMinor <- majorclass_with_minor %>%
  sapply(\(major) { fread(paste0("03_Marker_Expression/3_ovr_LogReg_minorclass-", major, "_xeniumFiltered.txt")) },
         USE.NAMES = TRUE, simplify = FALSE) %>%
  rbindlist()

# Munge ----

# Columns: Name, Parent_Name, Marker, Marker_Function, Marker_Source, Marker_Notes

## Curated ----

clean_curated = munge_curate(curated, verbose)

## Queried ----

setnames(queriedMinor, "Coefficient", "Minor_Coefficient")
queriedMinor[Name == "Microglia" & Minor_Coefficient < 0, Name := "CYCLING_MICROGLIA"]
queriedMinor[toupper(Name) %in% majorclass, Name := paste0("UNASSIGNED_", Name)]
queriedMinor[, Minor_Coefficient := abs(Minor_Coefficient)] # Any binary regression should have negative values

setnames(queriedMajor, "Coefficient", "Major_Coefficient")

clean_queried =
  rbindlist(list(queriedMajor, queriedMinor), fill = TRUE) %>%
  mutate(Queried = "Queried",
         Queried_Marker = Marker, # Preserve original names as it was in original data
         Queried_Name = Name,
         Queried_Major_Name = Major_Name)

## Harmonize the datasets ----
## Queried will mainly act as the base
## Curated will conform to base.

cell_clobber_check() # 173, 140, 16

### Gene Harmony ----

clean_curated = markers_to_mouse(clean_curated, verbose)
clean_queried = markers_to_mouse(clean_queried, verbose)

clean_curated = clean_genes(clean_curated, verbose) # Lost two unique genes from typos

### Cell Harmony ----

clean_curated = depluralize(clean_curated, verbose)
cell_clobber_check() # 173, 140, 20

clean_queried = set_name_case(clean_queried)
clean_curated = set_name_case(clean_curated)
cell_clobber_check() # 173, 140, 20

clean_curated = set_cellname_manually(clean_curated, CURATE2QUERY, verbose)
cell_clobber_check() # 173, 140, 28 = 20 + length(CURATE2QUERY)

clean_curated = clean_curated %>% filter(!duplicated(paste0(Name, Marker, Parent_Name))) # Remove odd duplicate

# Back to Query

clean_queried = query_sub(clean_queried)
cell_clobber_check() # 173, 140, 107

clean_queried = set_cellname_manually(clean_queried, QUERY2CURATE, verbose)
cell_clobber_check() # 173, 139, 128

# Any random cell types need manual adjustment?
clean_queried$Name %>% unique() %>% # "NOVEL_10","NT-OODSGC","NOVEL_13","OODS_CCK_14","DV-OODSGC","NOVEL_24","NOVEL_30","UNASSIGNED_AC","UNASSIGNED_BC","UNASSIGNED_RGC","UNASSIGNED_MICROGLIA"
  setdiff(c(clean_curated$Parent_Name, clean_curated$Name) %>% unique()) # Yes, the above 11 do

# Stop and Review ----

fwrite(clean_curated, paste0(outPath, "/4_curated.txt"), sep = '\t')
fwrite(clean_queried, paste0(outPath, "/4_queried.txt"), sep = '\t')

## Plot Overlap ----

plot_gene_cell_venn(clean_curated, clean_queried, outPath)

## Merge ----

clean_merged =
  merge(clean_curated, clean_queried, by = c("Marker", "Name"), all = TRUE, suffixes = c(".Curated", ".Queried")) %>%
  summarize(.by = c(Name, Marker),
            Curated, Queried, Queried_Name, Queried_Major_Name,
            Minor_Coefficient, Major_Coefficient)

all(unique(clean_merged$Parent_Name) %in% majorclass) # FALSE
all(majorclass %in% unique(clean_merged$Parent_Name)) # TRUE

clean_merged = get_major_name(clean_merged, verbose) # Check!!

## Save Results ----

clean_merged = clean_merged %>%
  reframe(Name, Marker, Curated, Queried, Minor_Coefficient, Major_Name, Major_Coefficient, Queried_Name) %>%
  arrange(Major_Name, Name, Marker) %T>%
  fwrite(paste0(outPath, "4_merged_curated-queried_markers.txt"), sep = '\t')

# break full table into files by major class
cell_outpath = paste0(outPath, "4_merged_curated-queried_majorclass/")
dir.create(cell_outpath, showWarnings = FALSE)
for (cell_type in unique(clean_merged$Major_Name)) {
  clean_merged %>%
    filter(Major_Name == cell_type) %>%
    reframe(Name, Marker, Curated, Queried, Minor_Coefficient, Major_Name, Major_Coefficient, Queried_Name) %>%
    arrange(Major_Name, Name, Marker) %>%
    fwrite(paste0(cell_outpath, "/4_merged_curated-queried_", cell_type, "_markers.txt"), sep = '\t')
}

### Search space overlap with curated ----

all_genes = scan("01_QualityControl/1_genes.txt", what = character()) %>% str_to_title()
variable_genes = scan("01_QualityControl/1_variable_genes.txt", what = character()) %>% str_to_title()
curated_markers = unique(clean_curated$Marker)

data.frame(Curated = curated_markers,
           All_Genes = curated_markers %in% all_genes,
           Highly_Variable = curated_markers %in% variable_genes) %>%
  fwrite(paste0(outPath, "4_curated_markers_in_data.txt"), sep = '\t')

### Queried_Name to Name for plotting ----

clean_merged[, c("Queried_Name", "Name", "Major_Name")] %>%
  na.omit() %>% distinct() %>%
  fwrite(paste0(outPath, "4_queried_to_name.txt"), sep = '\t')

# Scratch ----

if (FALSE) {

}
