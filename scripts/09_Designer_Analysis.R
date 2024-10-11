# JM 09/11/2024
# Look at the output of Xenium Panel Designer

# Setup ----

library(Seurat)
options(Seurat.object.assay.version = 'v5')
library(dplyr) # pipe
library(data.table) # fast table manip
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(magrittr)
library(readxl)
library(ChengLabThemes)

analysisName = "09_Designer_Analysis"
analysisPath = "/project/hipaa_ycheng11lab/atlas/CAMR2024/"
setwd(analysisPath)

outPath = paste0(analysisPath, analysisName, "/")
dir.create(outPath, showWarnings = FALSE)

## Reproducibility ----

now = \() Sys.time() %>% format('%F_%H-%M-%S-%Z')
sessionInfoPath = paste0("logs/", analysisName, ".", now(), ".sessionInfo")
sessionInfo() %>% capture.output(file = sessionInfoPath)

set.seed(1)

# Data ----

resultsPath = paste0(outPath, "PanelDesignerInputOutput.xlsx")
markers = list()
markers[["V1"]] = read_excel(resultsPath, sheet="V1")
markers[["V2"]] = read_excel(resultsPath, sheet="V2")

q2n = fread("04_Harmonize_Curated_Markers/4_queried_to_name.txt")
q2n[, Name := toupper(Name)]
q2n[, Queried_Name := toupper(Queried_Name)]

# Analysis ----

change_idx = markers[["V1"]]$Name %in% q2n$Queried_Name & markers[["V1"]]$Name != markers[["V1"]]$Major_Name
old_names = markers[["V1"]]$Name[change_idx]
markers[["V1"]]$Name[change_idx] = q2n$Name[match(old_names, q2n$Queried_Name)]

change_idx = markers[["V2"]]$Name %in% q2n$Queried_Name & markers[["V2"]]$Name != markers[["V2"]]$Major_Name
old_names = markers[["V2"]]$Name[change_idx]
markers[["V2"]]$Name[change_idx] = q2n$Name[match(old_names, q2n$Queried_Name)]

saveRDS(markers, paste0(outPath, "PanelOutput.RDS"))

qc_markers = markers[["V2"]] %>%
  filter(V2_Probecount > 0) %T>%
  fwrite("09_Designer_Analysis/PanelOutputV2.txt", sep = '\t')

major_marker_counts = qc_markers %>%
  filter(Name == Major_Name) %>%
  select(Major_Name, Marker) %>%
  distinct() %>%
  summarize(.by = Major_Name, count = n()) %T>% print() # all majornames are represented here

markers[["V1"]]$Marker[markers[["V1"]]$Major_Name == "ROD"] %in% qc_markers$Marker # 2
markers[["V1"]]$Marker[markers[["V1"]]$Major_Name == "EPITHELIAL"] %in% qc_markers$Marker # 0
markers[["V1"]]$Marker[markers[["V1"]]$Name == "DENDRITIC CELL"] %in% qc_markers$Marker # 2

p = ggplot(major_marker_counts, aes(x = reorder(Major_Name, -count), y = count)) +
  geom_col() +
  histogram_theme() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) +
  labs(y = "Major Class Markers", x = "")
ggsave(paste0(outPath, "Panel_Markers_Major.Bar.png"), p, height = 3, width = 5)

minor_marker_counts = qc_markers %>%
  filter(Name != Major_Name) %>%
  select(Major_Name, Name, Marker) %>%
  distinct() %>%
  summarize(.by = c(Major_Name, Name),
            count = n()) %T>% print()

for (major in unique(minor_marker_counts$Major_Name)) {
  p = ggplot(minor_marker_counts %>% filter(Major_Name == major),
             aes(x = reorder(Name, -count), y = count)) +
    geom_col() +
    histogram_theme() +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1)
    ) +
    labs(y = paste(major, "Minor Class Markers"), x = "")
  ggsave(paste0(outPath, "Panel_Markers_Major.", major, ".Bar.png"), p, height = 3, width = 10)
}

# Fill color and stack by Marker

# Scratch ----

if (FALSE) {
  # How many majorclasses are formatted correctly?
  # markers %>% filter(grepl("majorclass|immune", Source)) %>% reframe(a = Major_Name != Name) %>% {which(.$a)}

}
