#!/bin/env Rscript3

# JM 7/31/2024
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

analysisName = "04_Merge_Curated_Markers"
analysisPath = "/project/hipaa_ycheng11lab/atlas/CAMR2024/"
setwd(analysisPath)
source("scripts/04_Merge_Curated_Markers.util.R")

curatedPath = "data/CuratedMouseRetinaMarkers.txt"
queriedMajorPath = "03_Filter_Model_Markers/3_ovr_LogReg_majorclass_xeniumFiltered.txt"
queriedMinorPath = "03_Filter_Model_Markers/3_ovr_LogReg_minorclass_xeniumFiltered.txt"
major2MinorPath = "02_Modeling/2_minorToMajorClass.txt"

outPath = paste0(analysisPath, "04_Merge_Curated_Markers/")
dir.create(outPath, showWarnings = FALSE)
plotPath = paste0(outPath, "figures/")
dir.create(plotPath, showWarnings = FALSE)

verbose = TRUE

## Reproducible Results ----

now = \() Sys.time() %>% format('%F_%H-%M-%S-%Z')
sessionInfoPath = paste0("logs/", analysisName, ".", now(), ".sessionInfo")
sessionInfo() %>% capture.output(file = sessionInfoPath)

set.seed(1)

# Data ----

majorclass = c("AC", "ASTROCYTE", "BC", "CONE", "ENDOTHELIAL", "HC", "MG", "MICROGLIA", "PERICYTE", "RGC", "ROD", "RPE")
majorclass_with_minor = c("AC", "BC", "MICROGLIA", "RGC")

curated <- fread(curatedPath)
queriedMajor <- fread(queriedMajorPath)
queriedMinor <- fread(queriedMinorPath)
major2Minor <- fread(major2MinorPath, header = TRUE, select = c("author_cell_type", "majorclass"))

# Munge ----

# Columns: Name, Parent_Name, Marker, Marker_Function, Marker_Source, Marker_Notes

## Curated ----

clean_curated = munge_curate(curated, verbose)

# Next time, I think that only a few changes should be made to the original file
# to make it human understandable and machine understandable
# 1) columns: MajorClass, Class, Minor, Gene, Marker_Source, Marker_Annotation, Gene_Annotation (Optional), Notes
# 2) Be careful with copy and paste with regards to non-ASCII characters
# 3) No need for blank rows (Optional)
# 4) No headers with spaces, underscores as the only non-alphanumeric
# 5) NA instead of n/a
# extra) Columns should have 1 type of thing. It doesn't matter what it is, as long as it's the same.

## Queried ----

setnames(queriedMinor, "Coefficient", "Minor_Coefficient")
queriedMinor[Name == "Microglia" & Minor_Coefficient < 0, Name := "CYCLING_MICROGLIA"]
queriedMinor[toupper(Name) %in% majorclass_with_minor, Name := paste0("UNASSIGNED_", Name)]
# queriedMinor[toupper(Name) %in% setdiff(majorclass, majorclass_with_minor), Name := paste0("MINOR_", Name)]
queriedMinor[, Minor_Coefficient := abs(Minor_Coefficient)] # Any binary regression should have negative values

setnames(queriedMajor, "Coefficient", "Major_Coefficient")

clean_queried =
  rbindlist(list(queriedMajor, queriedMinor), fill = TRUE) %>%
  mutate(Queried = "Queried",
         Queried_Name = Name, # Preserve original names as it was in original data
         Queried_Marker = Marker,
         Queried_Major_Name = Major_Name)

major2Minor %<>%
  reframe(Name = author_cell_type,
          Major_Name = majorclass,
          Queried_Name = author_cell_type,
          Queried_Major_Name = majorclass) %>% asDT()
setDT(major2Minor)
major2Minor[toupper(Name) %in% majorclass_with_minor, Name := paste0("UNASSIGNED_", Name)]

## Harmonize the datasets ----
## Queried will mainly act as the base as it matches original data.
## Curated will conform to base.

cell_clobber_check() # 173, 48, 9

### Gene Harmony ----

clean_curated = markers_to_mouse(clean_curated, verbose) # 248, 247 bc of space
clean_queried = markers_to_mouse(clean_queried, verbose) # 217

clean_curated = clean_genes(clean_curated, verbose) # 247 -> 245, Lost two unique genes from typos

### Cell Harmony ----

clean_curated = depluralize(clean_curated, verbose) # 114, 317 parents
cell_clobber_check() # 173, 48, 13

clean_queried = set_name_case(clean_queried)
clean_curated = set_name_case(clean_curated)
major2Minor = set_name_case(major2Minor)
cell_clobber_check() # 173, 48, 13

clean_curated = set_cellname_manually(clean_curated, CURATE2QUERY, verbose)
cell_clobber_check() # 173, 48, 21 = 13 + length(CURATE2QUERY)

clean_curated = clean_curated %>% filter(!duplicated(paste0(Name, Marker, Parent_Name))) # Remove odd duplicate

# Back to Query

clean_queried = query_sub(clean_queried)
major2Minor = query_sub(major2Minor)
cell_clobber_check() # 173, 48, 37

clean_queried = set_cellname_manually(clean_queried, QUERY2CURATE, verbose)
major2Minor = set_cellname_manually(major2Minor, QUERY2CURATE, verbose)
major2Minor[, Parent_Name := NULL]
cell_clobber_check() # 173, 48, 45

# Any random cell types need manual adjustment?
q_names = clean_queried$Name %>% unique()# "NOVEL_13" "UNASSIGNED_RGC" "UNASSIGNED_MICROGLIA"
c_names = c(clean_curated$Parent_Name, clean_curated$Name) %>% unique()
setdiff(q_names, c_names) # Yes, the above 3 do
setdiff(c_names, q_names) # Yes, 132 types are missing

## Stop and Review ----

fwrite(clean_queried, paste0(outPath, "4_harmonized_queried_markers.txt"), sep = '\t')

## Plot Overlap ----

plot_gene_cell_venn(clean_curated, clean_queried, outPath)

## Harmonize Major_Name ----

all(unique(clean_curated$Parent_Name) %in% majorclass) # FALSE
all(majorclass %in% unique(clean_curated$Parent_Name)) # FALSE

clean_curated = get_major_name(clean_curated, verbose)
all(majorclass %in% unique(clean_curated$Major_Name)) # TRUE

## Harmonize Query Name ----
major2Minor[, Major_Name := toupper(Major_Name)]
fwrite(major2Minor, paste0(outPath, "4_queried_to_name.txt"), sep = '\t')

clean_curated %<>% reframe(Name, Parent_Name, Marker, Major_Name)

clean_curated = major2Minor[clean_curated, on = c("Name", "Major_Name")]
setDT(clean_curated)
clean_curated
clean_curated[Name %in% majorclass, Queried_Name := Parent_Name]
clean_curated[Name %in% majorclass, Queried_Major_Name := Parent_Name]
clean_curated
clean_curated[Major_Name %in% majorclass, Queried_Major_Name := Major_Name]
clean_curated

## Finished ----

clean_curated = clean_curated %>%
  reframe(Name, Major_Name, Marker, Queried_Name, Queried_Major_Name) %>%
  arrange(Major_Name, Name, Marker) %T>%
  fwrite(paste0(outPath, "4_harmonized_curated_markers.txt"), sep = '\t')

### Search space overlap with curated ----

all_genes = scan("01_QualityControl/1_genes.txt", what = character()) %>% str_to_title()
variable_genes = scan("01_QualityControl/1_variable_genes.txt", what = character()) %>% str_to_title()
curated_markers = unique(clean_curated$Marker)

data.frame(Curated = curated_markers,
           All_Genes = curated_markers %in% all_genes,
           Highly_Variable = curated_markers %in% variable_genes) %>%
  fwrite(paste0(outPath, "4_curated_markers_in_data.txt"), sep = '\t')

### Queried_Name to Name for plotting ----



# Scratch ----

if (FALSE) {

  ## Inheritance ----

  # Try to cascade the gene

  # DFS
  while (nrow(sentinel)) {
    clean_curated = rbindlist(list(clean_curated, clean_curated$Marker[clean_curated$Parent_Name[sentinel]]))
  }

}
