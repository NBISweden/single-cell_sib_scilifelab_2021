import anndata as ad
import pandas as pd
import os.path as osp
from scipy.sparse import csr_matrix

DATA_DIR = "PUT_PATH_TO_DIRECTORY_WHERE_RAW_DATA_IS_FOUND"
OUT_DIR = "PATH_TO_DIRECTORY_WHERE_DATA_SHOULD_BE_SAVED"

EXON_PTH = osp.join(DATA_DIR, "GSE115746_cells_exon_counts.csv")
META_PTH = osp.join(DATA_DIR, "GSE115746_complete_metadata_28706-cells.csv")

read_file = lambda f: pd.read_csv(f, sep=",", header=0, index_col=0)

cnt = read_file(EXON_PTH)
meta = read_file(META_PTH)
cnt = cnt.T

inter = cnt.index.intersection(meta.index)
cnt = cnt.loc[inter, :]
meta = meta.loc[inter, :]

var = pd.DataFrame(
    cnt.columns.values,
    index=cnt.columns,
    columns=["gene"],
)

adata = ad.AnnData(
    csr_matrix(cnt.values),
    var=var,
    obs=meta,
)
del cnt

adata.write_h5ad(osp.join(OUT_DIR,"VISp-sc_data.h5ad"))
