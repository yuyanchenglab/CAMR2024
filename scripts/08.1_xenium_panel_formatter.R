# Load libraries
library(Seurat)
library(dplyr)

analysisPath = "/project/hipaa_ycheng11lab/atlas/CAMR2024/08_xenium_panel_formatter/"
setwd(analysisPath)

# Read in the dataset
seurat_obj <- readRDS("/project/hipaa_ycheng11lab/atlas/CAMR2024/06_Curated_Immune/6_immune.RDS") #readRDS("/project/ycheng11lab/ychengprj01/pub-rgcs/cellChat/230720_Merged_SCTrfm.RDS")
all.equal(GetAssayData(object=seurat_obj, assay="RNA", slot="counts")@x,
          as.integer(GetAssayData(object=seurat_obj, assay="RNA", slot="counts")@x))

# Output
## [1] TRUE
head(GetAssayData(object=seurat_obj, assay="RNA", slot="counts")@x)

# Output
## [1]  1  2 18  1  2  1
# Check row names
head(rownames(seurat_obj))

# Output
## [1] "ENSG00000223972" "ENSG00000227232" "ENSG00000278267" "ENSG00000243485"
## [5] "ENSG00000284332" "ENSG00000237613"

# Check meta features
dplyr::glimpse(GetAssay(seurat_obj)@meta.features)

# Output
## Rows: 58,604
## Columns: 12
## $ feature_type          <fct> Gene Expression, Gene Expression, Gene Expressio…
## $ highly.variable       <lgl> FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,…
## $ mvp.mean              <dbl> 6.398244e-05, 2.274395e-03, 6.175251e-05, 1.3728…
## $ mvp.dispersion        <dbl> 0.8350443, 2.4422800, 1.2953346, 2.6563521, NaN,…
## $ mvp.dispersion.scaled <dbl> -0.5739472, 0.5332035, -0.2568744, 0.6806679, 0.…
## $ mean                  <dbl> 3.881082e-05, 1.079995e-03, 3.267138e-05, 4.7735…
## $ std                   <dbl> 5.573522e-03, 3.173095e-02, 5.634017e-03, 8.0407…
## $ ensembl_version       <chr> "ENSG00000223972.5", "ENSG00000227232.5", "ENSG0…
## $ feature_is_filtered   <lgl> FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,…
## $ feature_name          <fct> DDX11L1, WASH7P, MIR6859-1, MIR1302-2HG, MIR1302…
## $ feature_reference     <fct> NCBITaxon:9606, NCBITaxon:9606, NCBITaxon:9606, …
## $ feature_biotype       <fct> gene, gene, gene, gene, gene, gene, gene, gene, …
# Check metadata
dplyr::glimpse(seurat_obj[[]])

# Output
## Rows: 13,497
## Columns: 26
## $ assay_ontology_term_id                   <fct> EFO:0009922, EFO:0009922, EFO…
## $ donor_id                                 <fct> TSP9, TSP9, TSP9, TSP9, TSP9,…
## $ anatomical_information                   <fct> exocrine, exocrine, exocrine,…
## $ nCounts_RNA_UMIs                         <dbl> 26642, 4053, 8439, 7967, 9405…
## $ nFeaturess_RNA                           <dbl> 3833, 1624, 3067, 1484, 3273,…
## $ cell_ontology_class                      <fct> pancreatic acinar cell, t cel…
## $ free_annotation                          <fct> pancreatic acing cell, T cell…
## $ manually_annotated                       <lgl> TRUE, TRUE, TRUE, TRUE, TRUE,…
## $ compartment                              <fct> epithelial, immune, endotheli…
## $ sex_ontology_term_id                     <fct> PATO:0000384, PATO:0000384, P…
## $ disease_ontology_term_id                 <fct> PATO:0000461, PATO:0000461, P…
## $ is_primary_data                          <lgl> FALSE, FALSE, FALSE, FALSE, F…
## $ organism_ontology_term_id                <fct> NCBITaxon:9606, NCBITaxon:960…
## $ suspension_type                          <fct> cell, cell, cell, cell, cell,…
## $ cell_type_ontology_term_id               <fct> CL:0002064, CL:0000084, CL:00…
## $ tissue_ontology_term_id                  <fct> UBERON:0000017, UBERON:000001…
## $ development_stage_ontology_term_id       <fct> HsapDv:0000131, HsapDv:000013…
## $ self_reported_ethnicity_ontology_term_id <fct> HANCESTRO:0014, HANCESTRO:001…
## $ cell_type                                <fct> pancreatic acinar cell, T cel…
## $ assay                                    <fct> 10x 3' v3, 10x 3' v3, 10x 3' …
## $ disease                                  <fct> normal, normal, normal, norma…
## $ organism                                 <fct> Homo sapiens, Homo sapiens, H…
## $ sex                                      <fct> male, male, male, male, male,…
## $ tissue                                   <fct> exocrine pancreas, exocrine p…
## $ self_reported_ethnicity                  <fct> Hispanic or Latin American, H…
## $ development_stage                        <fct> 37-year-old human stage, 37-y…
subsample_rate <- 1

subset <- 1:ncol(seurat_obj)
if (subsample_rate < 1) {
  subset <- sample(subset, subsample_rate * length(subset))
}
# Install DropletUtils
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("DropletUtils", quietly = TRUE))
  BiocManager::install("DropletUtils")

# Load library
library(DropletUtils)
# Run function
write10xCounts(
  "reference_data",
  GetAssayData(seurat_obj, assay="RNA", slot="counts")[rownames(seurat_obj), subset],
  gene.id = rownames(seurat_obj),
  gene.symbol = rownames(seurat_obj),
  barcodes = colnames(seurat_obj)[subset],
  type = "sparse",
  version = "3"
)

# List the files
list.files("reference_data")

# Output
#[1] "barcodes.tsv.gz" "features.tsv.gz" "matrix.mtx.gz"
# Define function
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
# Run function
bundleOutputs(out_dir = "/project/hipaa_ycheng11lab/atlas/CAMR2024/08_xenium_panel_formatter/reference_data/", data = seurat_obj, subset = subset, cell_type = "Major_Name")
# List files
list.files("reference_data")

# Output
#[1] "annotations.csv"    "barcodes.tsv.gz"    "features.tsv.gz"    "matrix.mtx.gz"      "reference_data.zip"
