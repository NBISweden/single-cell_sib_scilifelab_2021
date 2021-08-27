import scvi
import scanpy as sc

sc.set_figure_params(figsize=(4, 4))
adata = scvi.data.heart_cell_atlas_subsampled()
sc.pp.filter_genes(adata, min_counts=3)

adata.layers["counts"] = adata.X.copy() # preserve counts
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata # freeze the state in `.raw`

sc.pp.highly_variable_genes(
  adata,
  n_top_genes=1200,
  subset=True,
  layer="counts",
  flavor="seurat_v3",
  batch_key="cell_source"
)

scvi.data.setup_anndata(
    adata,
    layer="counts",
    categorical_covariate_keys=["cell_source", "donor"],
    continuous_covariate_keys=["percent_mito", "percent_ribo"]
)



model = scvi.model.SCVI(adata)
model
model.train(max_epochs=40)




trainer.train(n_epochs=10, lr=lr)


latent = model.get_latent_representation()
adata.obsm["X_scVI"] = latent

adata_subset = adata[adata.obs.cell_type == "Fibroblast"]
latent_subset = model.get_latent_representation(adata_subset)

denoised = model.get_normalized_expression(adata_subset, library_size=1e4)
denoised.iloc[:5, :5]

adata.layers["scvi_normalized"] = model.get_normalized_expression(
    library_size=10e4
)

# use scVI latent space for UMAP generation
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.umap(adata, min_dist=0.3

sc.pl.umap(
    adata,
    color=["cell_type"],
    frameon=False,
)

sc.pl.umap(
    adata,
    color=["donor", "cell_source"],
    ncols=2,
    frameon=False,
)


