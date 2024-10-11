# This script is how we will prepare the dataset to be uploaded to the 10X
# Xenium Custom Panel Design software to use to optimize our gene panel for our
# future experiments. It will also be used to downsample each cluster of cells
# to get a possibly clearer picture of the genes that work best in our biology.

# Setup ----

library(Seurat)
library(dplyr)
library(DropletUtils)

analysisPath = "/project/ycheng11lab/jfmaurer/mouse_retina_atlas_chen_2024/08_xenium_panel_formatter/"
setwd(analysisPath)

subsample_rate <- 1

bundleOutputs <- function(out_dir, data, barcodes = colnames(data), cell_type = "cell_type", subset = 1:length(barcodes)) {

  if (require("data.table", quietly = TRUE)) {
    data.table::fwrite(
      data.table::data.table(
        barcode = barcodes,
        annotation = unlist(seurat_obj[[cell_type]])
      )[subset, ],
      file.path(out_dir, "annotations.csv")
    )
  } else {
    write.table(
      data.frame(
        barcode = barcodes,
        annotation = unlist(seurat_obj[[cell_type]])
      )[subset, ],
      file.path(out_dir, "annotations.csv"),
      sep = ",", row.names = FALSE
    )
  }

  bundle <- file.path(out_dir, paste0(basename(out_dir), ".zip"))

  utils::zip(
    bundle,
    list.files(out_dir, full.names = TRUE),
    zip = "zip"
  )

  if (file.info(bundle)$size / 1e6 > 500) {
    warning("The output file is more than 500 MB and will need to be subset further.")
  }
}

# Work ----

seurat_obj <- readRDS("/project/ycheng11lab/jfmaurer/mouse_retina_atlas_chen_2024/06_Curated_Immune/6_immune.RDS")

if (interactive()) {
# Check for counts
all.equal(GetAssayData(object=seurat_obj, assay="RNA", slot="counts")@x,
          as.integer(GetAssayData(object=seurat_obj, assay="RNA", slot="counts")@x))
## [1] TRUE

# Check row names
head(rownames(seurat_obj))
## [1] "ENSG00000223972" "ENSG00000227232" "ENSG00000278267" "ENSG00000243485"

# Check meta features
dplyr::glimpse(GetAssay(seurat_obj)@meta.features)
## Rows: 58,604
## Columns: 12
## $ feature_type          <fct> Gene Expression, Gene Expression, Gene Expressio…
## $ highly.variable       <lgl> FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,…
## $ mvp.mean              <dbl> 6.398244e-05, 2.274395e-03, 6.175251e-05, 1.3728…
## $ mvp.dispersion        <dbl> 0.8350443, 2.4422800, 1.2953346, 2.6563521, NaN,…
## $ ...

# Check metadata
dplyr::glimpse(seurat_obj[[]])
## Rows: 13,497
## Columns: 26
## $ ...
## $ nCounts_RNA_UMIs                         <dbl> 26642, 4053, 8439, 7967, 9405…
## $ nFeaturess_RNA                           <dbl> 3833, 1624, 3067, 1484, 3273,…
## $ cell_ontology_class                      <fct> pancreatic acinar cell, t cel…
## $ ...
}

subset <- 1:ncol(seurat_obj)
if (subsample_rate < 1) {
  subset <- sample(subset, subsample_rate * length(subset))
}

write10xCounts(
  "reference_data",
  GetAssayData(seurat_obj, assay="RNA", slot="counts")[rownames(seurat_obj), subset],
  gene.id = rownames(seurat_obj),
  gene.symbol = rownames(seurat_obj),
  barcodes = colnames(seurat_obj)[subset],
  type = "sparse",
  version = "3")

if (interactive()) {
list.files("reference_data")
#[1] "barcodes.tsv.gz" "features.tsv.gz" "matrix.mtx.gz"
}

bundleOutputs(out_dir = "/project/hipaa_ycheng11lab/atlas/CAMR2024/08_xenium_panel_formatter/reference_data/", data = seurat_obj, subset = subset, cell_type = "Major_Name")

if (interactive()) {
list.files("reference_data")
#[1] "annotations.csv"    "barcodes.tsv.gz"    "features.tsv.gz"    "matrix.mtx.gz"      "reference_data.zip"
}
