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
source("scripts/CurationVsQuery.util.R")

analysisName = "CurationVsQueryByCell"
analysisPath = "/project/hipaa_ycheng11lab/atlas/CAMR2024/"
setwd(analysisPath)

curatedPath = "spreadsheets/CuratedMouseRetinaMarkers.txt"
queriedMajorPath = "spreadsheets/ovr_top_filtered_genes_majorclass_coefficients_sensitive.csv"

plotPath = paste0("figures/", analysisName, '/')
dir.create(plotPath, showWarnings = FALSE)

verbose = TRUE

## Reproducible Results ----

now = \() Sys.time() %>% format('%F_%H-%M-%S-%Z')
sessionInfoPath = paste0("logs/", analysisName, ".", now(), ".sessionInfo")
sessionInfo() %>% capture.output(file = sessionInfoPath)

set.seed(1)

# Data ----

majorclass = c("AC", "ASTROCYTE", "BC", "CONE", "ENDOTHELIAL", "HC", "MG", "MICROGLIA", "PERICYTE", "RGC", "ROD", "RPE")

needed_major = c("AC", "BC", "Microglia", "RGC")

curated <- fread(curatedPath)
queriedMajor <- fread(queriedMajorPath,
                      select = c("Gene", "Cell Type", "Coefficient"),
                      col.names = c("Marker", "Name", "Major_Coefficient"))
queriedMinor <- sapply(needed_major,
                       \(major) { fread(paste0("spreadsheets/ovr_top_", "complete", "_filtered_markers_", major, "_coefficients_sensitive.csv"), drop = c(1,4),
                                        col.names = c("Parent_Name", "Name", "Marker", "Minor_Coefficient")) },
                       USE.NAMES = TRUE, simplify = FALSE) %>%
  rbindlist()

# Munge ----

# Columns: Name, Parent_Name, Marker, Marker_Function, Marker_Source, Marker_Notes

## Curated ----

clean_curated = munge_curate(curated, verbose)

## Queried ----

queriedMinor[Name == "Microglia" & Minor_Coefficient < 0, Name := "CYCLING_MICROGLIA"]
queriedMinor[toupper(Name) %in% majorclass, Name := paste0("UNASSIGNED_", Name)]
queriedMinor[, Minor_Coefficient := abs(Minor_Coefficient)]

queriedMajor[, Parent_Name := Name]

clean_queried =
  rbindlist(list(queriedMajor, queriedMinor), fill = TRUE) %>%
  mutate(Queried = "Queried",
         Queried_Marker = Marker,
         Queried_Name = Name,
         Queried_Parent_Name = Parent_Name)

## Harmonize the datasets ----
## Queried will mainly act as the base
## Curated will conform to base.

cell_clobber_check() # 173, 140, 16

### Gene Harmony ----

clean_curated = markers_to_mouse(clean_curated, verbose)
clean_queried = markers_to_mouse(clean_queried, verbose)
clean_curated = clean_genes(clean_curated, verbose) # Lost two genes from typos

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
cell_clobber_check() # 107 / 140

clean_queried = set_cellname_manually(clean_queried, QUERY2CURATE, verbose)
cell_clobber_check() # 173, 139, 128

# Any random cell types need manual adjustment?
clean_queried$Name %>% unique() %>% # "NOVEL_10","NT-OODSGC","NOVEL_13","OODS_CCK_14","DV-OODSGC","NOVEL_24","NOVEL_30","UNASSIGNED_AC","UNASSIGNED_BC","UNASSIGNED_RGC","UNASSIGNED_MICROGLIA"
  setdiff(c(clean_curated$Parent_Name, clean_curated$Name) %>% unique()) # Yes, the above 11 do

# Stop and Review ----

fwrite(clean_curated, "spreadsheets/curated_table.txt")
fwrite(clean_queried, "spreadsheets/queried_by_cell_table.txt")

## Plot Overlap ----

plot_gene_cell_venn(clean_curated, clean_queried)

## Merge ----

clean_merged =
  merge(clean_curated, clean_queried, by = c("Marker", "Name"), all = TRUE, suffixes = c(".Curated", ".Queried")) %>%
  summarize(.by = c(Name, Marker),
            Parent_Name = ifelse(is.na(Parent_Name.Curated), Parent_Name.Queried, Parent_Name.Curated),
            Curated, Queried, Queried_Name, Queried_Parent_Name,
            Minor_Coefficient, Major_Coefficient)

all(unique(clean_merged$Parent_Name) %in% majorclass) # FALSE
all(majorclass %in% unique(clean_merged$Parent_Name)) # TRUE

clean_merged = get_major_name(clean_merged, verbose)

### Save Results ----

clean_merged = clean_merged %>%
  reframe(Name, Marker, Curated, Queried, Minor_Coefficient, Major_Name, Major_Coefficient) %>%
  arrange(Major_Name, Name, Marker) %T>%
  fwrite("spreadsheets/merged_curated-queried_ByCell_markers_with_Coefficients.txt", sep = '\t')

# break full table into files by major class
for (cell_type in unique(clean_merged$Major_Name)) {
  clean_merged %>%
    filter(Major_Name == cell_type) %>%
    reframe(Name, Marker, Curated, Queried) %>%
    fwrite(paste0("spreadsheets/merged_curated-queried_", cell_type, "_markers_with_Coefficients.txt"), sep = '\t')
}

### Search space overlap with curated ----

all_genes = fread("data/CAMR_all_genes.csv")
all_genes_formatted = str_to_title(all_genes$feature_name)
variable_genes = fread("data/CAMR_genes.csv")
variable_genes_formatted = str_to_title(variable_genes$feature_name)
curated_markers = unique(clean_curated$Marker)

data.frame(Curated = curated_markers,
           All_Genes = curated_markers %in% all_genes_formatted,
           Highly_Variable = curated_markers %in% variable_genes_formatted) %>%
  fwrite("spreadsheets/genes_in_curated_markers.txt", sep = '\t')

# Scratch ----

if (FALSE) {

}
