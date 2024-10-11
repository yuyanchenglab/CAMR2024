# JM 09/12/2024
# Make a shiny app for the data I work on.

# Setup ----

library(Seurat)
library(ShinyCell)
library(reticulate)
# module load hdf5/1.12.1 is needed to install ShinyCell

analysisName = "10_Make_Shiny"
analysisPath = "/project/ycheng11lab/jfmaurer/mouse_retina_atlas_chen_2024/"
outPath = paste0(analysisPath, analysisName, "/")
dir.create(outPath, showWarnings = FALSE)
setwd(outPath)

# This must be installed...
if (!interactive()) {
  reticulate::use_python("/appl/python-3.11/bin/python")
}

# Make App ----

## Standalone ----

print("Benhar Tran Immun Alone")

seu = readRDS('/project/ycheng11lab/jfmaurer/mouse_retina_atlas_chen_2024/06_Curated_Immune/6_immune.RDS')
scConf = createConfig(seu)
makeShinyApp(seu, scConf,
             shiny.title = "Time Series Retina Crush -- Immune Subset",
             shiny.dir = "TimeSeriesRetinaCrushImmune",
             default.gene1 = "Cd74", default.gene2 = "Cd3e",
             default.multigene = c("Serpinb1b", "Serpine2", "Slc17a1", "Slc17a6", "Slc24a2", "Slc7a11", "Slc7a6", "Stxbp6", "Tac1", "Tagln2", "Tbr1", "Tbx20", "Thy1", "Tpbg"),
             gex.slot = "counts")

print("Mouse Retina Atlas Alone")

seu = '/project/ycheng11lab/jfmaurer/mouse_retina_atlas_chen_2024/10_Make_Shiny/10_Shiny_Input.h5ad'
scConf = createConfig(seu, meta.to.include = c("nCount_RNA", "nFeature_RNA", "percent.mt", "majorclass", "reference", "minorclass", "Major_Name"), maxLevels = 150)
makeShinyApp(seu, scConf,
             shiny.title = "Mouse Retina Atlas",
             shiny.dir = "MouseRetinaAtlas",
             default.gene1 = "Cd74", default.gene2 = "Cd3e",
             default.multigene = c("Tbr1", "Tbx20", "Thy1", "Tpbg"))

## Combined ----

### CAMR ----

print("CAMR")

seuAtlas = '/project/ycheng11lab/jfmaurer/mouse_retina_atlas_chen_2024/10_Make_Shiny/10_Shiny_Input.h5ad'
scConfAtlas = createConfig(seuAtlas)
makeShinyFiles(seuAtlas, scConfAtlas,
             shiny.prefix = "sc1",
             default.gene1 = "Cd74", default.gene2 = "Cd3e",
             default.multigene = c("Serpinb1b", "Serpine2", "Slc17a1", "Slc17a6", "Slc24a2", "Slc7a11", "Slc7a6", "Stxbp6", "Tac1", "Tagln2", "Tbr1", "Tbx20", "Thy1", "Tpbg"),
             shiny.dir = "MouseRetinaMarkers")

### Immune Standalone ----

print("Immune")

seuImmune = readRDS('/project/ycheng11lab/jfmaurer/mouse_retina_atlas_chen_2024/06_Curated_Immune/6_immune.RDS')
scConfImmune = createConfig(seuImmune)
makeShinyFiles(seuImmune, scConfImmune,
             shiny.prefix = "sc2",
             default.gene1 = "Cd74", default.gene2 = "Cd3e",
             default.multigene = c("Serpinb1b", "Serpine2", "Slc17a1", "Slc17a6", "Slc24a2", "Slc7a11", "Slc7a6", "Stxbp6", "Tac1", "Tagln2", "Tbr1", "Tbx20", "Thy1", "Tpbg"),
             shiny.dir = "MouseRetinaMarkers")

### Combine ----

print("Combine")

makeShinyCodesMulti(shiny.title = "Xenium Mouse Retina Marker Search",
                    shiny.prefix = c("sc1", "sc2"),
                    shiny.headers = c("Mouse Retina Atlas", "ONC Immune Subset"),
                    shiny.dir = "MouseRetinaMarkers", shiny.footnotes = "")

# Scratch ----

if (FALSE) {

}
