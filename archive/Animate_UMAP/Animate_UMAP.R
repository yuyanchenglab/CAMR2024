# JM 6/25/2024
# View how the harmony software screwed up the UMAP of the data

# Setup ----

library(sceasy) # Anndata to Seurat
library(data.table)
library(ggplot2)
library(plotly)
library(magrittr) # %T>%
library(htmlwidgets) # saveWidget, need to module load pandoc before starting R on command line
# library(gapminder) # test dataset

setwd("/project/hipaa_ycheng11lab/atlas/CAMR2024/Animate_UMAP/")

# Make sure we have the version info of everything we used
now = \() Sys.time() %>% format('%F_%H-%M-%S-%Z')
sessionInfo() %>%
  capture.output(file = paste("animate_UMAP", now(), "sessionInfo", sep = "."))

# Data Munging ----

# too big for rstudio, ran on command line interactive job
reticulate::use_python("/appl/python-3.11/bin/python")
CAMR = convertFormat("../camr_scrublet_batch_filtered.h5ad",
                      from="anndata", to="seurat",
                      outFile = "camr_filtered.RDS")

# Extract the relevant information

oldUMAP = fread("CAMR_norm_UMAP.txt")
newUMAP = fread("CAMR_norm_UMAP_harmony.txt")

setnames(oldUMAP, new = c("UMAP1","UMAP2"))
setnames(newUMAP, new = c("UMAP1","UMAP2"))

oldUMAP$harmony = FALSE
newUMAP$harmony = TRUE

# Adjust UMAPs to the same orientation by flipping one of them horizontally/vertically
newUMAP[, UMAP1 := -UMAP1]

metadata = fread("CAMR_metadata.txt")
setnames(metadata, old = "V1", "cell_ID")
oldUMAP$cell_ID = newUMAP$cell_ID = metadata$cell_ID

UMAP = rbind(oldUMAP, newUMAP) # long for ggplot
UMAP = merge(UMAP, metadata, by = "cell_ID")

rm(metadata, oldUMAP, newUMAP); gc()

# Visualize ----

p <- ggplot(UMAP, aes(x = UMAP1, y = UMAP2, color = reference)) +
  geom_point(aes(frame = harmony, ids = cell_ID), size = 1)
fig <- ggplotly(p) %>%
  animation_opts(frame = 5000, easing = "cubic-in-out") %T>%
  saveWidget(file = "harmony_bastardization_project.html", selfcontained = TRUE)

# subset
set.seed(1)
some_cells = sample(UMAP$cell_ID, 1000)
p <- ggplot(UMAP %>% filter(cell_ID %in% some_cells),
            aes(x = UMAP1, y = UMAP2, color = reference)) +
  geom_point(aes(frame = harmony, ids = cell_ID), size = 1)
fig <- ggplotly(p) %>%
  animation_opts(frame = 5000, easing = "cubic-in-out") %T>%
  saveWidget(file = "harmony_bastardization_project_small.html", selfcontained = TRUE)

set.seed(1)
some_cells = sample(unique(UMAP$cell_ID), 10*1000)
p <- ggplot(UMAP %>% filter(cell_ID %in% some_cells),
            aes(x = UMAP1, y = UMAP2, color = reference)) +
  geom_point(aes(frame = harmony, ids = cell_ID), size = 1)
fig <- ggplotly(p) %>%
  animation_opts(frame = 5000, easing = "cubic-in-out") %T>%
  saveWidget(file = "harmony_bastardization_project_medium.html", selfcontained = TRUE)

set.seed(1)
some_cells = sample(unique(UMAP$cell_ID), 100*1000)
p <- ggplot(UMAP %>% filter(cell_ID %in% some_cells),
            aes(x = UMAP1, y = UMAP2, color = reference)) +
  geom_point(aes(frame = harmony, ids = cell_ID), size = 0.2)
fig <- ggplotly(p) %>%
  animation_opts(frame = 5000, easing = "cubic-in-out") %T>%
  saveWidget(file = "harmony_bastardization_project_large.html", selfcontained = TRUE)

if (interactive()) { print(fig) }

# Scratch ----

