# This script is how we will prepare the dataset to be uploaded to the 10X
# Xenium Custom Panel Design software to use to optimize our gene panel for our
# future experiments. It will also be used to downsample each cluster of cells
# to get a possibly clearer picture of the genes that work best in our biology.

# Load the h5ad file using scanpy
import scanpy as sc
import os

os.chdir('/project/hipaa_ycheng11lab/atlas/CAMR2024')
newpath = '/project/hipaa_ycheng11lab/atlas/CAMR2024/08_xenium_panel_formatter/'
if not os.path.exists(newpath):
    os.makedirs(newpath)

# Edit the file path in this command to point to the h5ad file on your computer
ad = sc.read_h5ad("00_raw/d0183df5-815d-48c2-bcfe-fbf9b716505c.h5ad")

# rewrite with the chunked work in the get_marker_genes.ipynb, or run on big worker node.
raw_ad = ad.raw.to_adata()

# Print information for the first row
#raw_ad.var.head(1)
raw_ad.var["ensemblid"] = raw_ad.var_names

# Import libraries
import pandas as pd
import scipy
import zipfile
import tempfile
import os
import gzip
import random

random.seed(1234)

# Define function
def h5ad_to_10x(ad,
                gene_id_key="gene_id",
                gene_name_key="gene_name",
                cell_type_key="cell_type",
                output_path="matrix.zip",
                barcode_key=None,
                subsample_rate=None):
    """
    Convert an AnnData h5ad file to 10x format,
    then produce a zipfile of the final matrix.

    Arguments:
        ad: A scanpy AnnData object.
        gene_id_key: The key that the gene IDs are under.
        gene_name_key: The key that the gene names are under.
        cell_type_key: The key that the cell types are under.
        barcode_key: Optional key that the barcodes are under. If not set,
            will assume they are the index of `obs`.
        output_path: Path to write the zipfile to.
        subsample_rate: Optional argument for subsampling.
            If provided, should be between 0 and 1.
    """
    if subsample_rate:
        sc.pp.subsample(ad, subsample_rate)

    genes = ad.var.reset_index()[[gene_id_key, gene_name_key]]
    genes["feature_type"] = ["Gene Expression"] * len(genes)

    if barcode_key:
        barcodes = ad.obs[[barcode_key]]
    else:
        barcodes = pd.DataFrame(ad.obs.index)

    celltypes = ad.obs[[cell_type_key]].reset_index()
    celltypes.columns = ["barcode", "annotation"]

    with tempfile.TemporaryDirectory() as tmp_dir:

        with gzip.open(os.path.join(tmp_dir, "matrix.mtx.gz"), "w") as handle:
            scipy.io.mmwrite(handle, ad.X.T.astype(int))

        genes.to_csv(os.path.join(tmp_dir, "features.tsv.gz"), sep="\t", index=False, header=False, compression="gzip")
        barcodes.to_csv(os.path.join(tmp_dir, "barcodes.tsv.gz"), sep="\t", index=False, header=False, compression="gzip")
        celltypes.to_csv(os.path.join(tmp_dir, "celltypes.csv"), index=False)

        with zipfile.ZipFile(output_path, "w") as zip_handle:
            for file in ["matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz", "celltypes.csv"]:
                zip_handle.write(os.path.join(tmp_dir, file), arcname=file)

# Run function to save the matrix as a MEX file using the keys we identified earlier
h5ad_to_10x(raw_ad, gene_id_key="ensemblid", gene_name_key="feature_name", cell_type_key="majorclass", subsample_rate = 0.05)

os.system("unzip -l matrix.zip")

# Output
#Archive:  matrix.zip
#  Length      Date    Time    Name
#---------  ---------- -----   ----
#151244754  11-09-2023 12:14   matrix.mtx.gz
#   479701  11-09-2023 12:14   features.tsv.gz
#    72881  11-09-2023 12:14   barcodes.tsv.gz
#   941922  11-09-2023 12:14   celltypes.csv
#---------                     -------
#152739258                     4 files

# Now you are ready to upload the .zip file to the Xenium Panel Designer to build your panel.







