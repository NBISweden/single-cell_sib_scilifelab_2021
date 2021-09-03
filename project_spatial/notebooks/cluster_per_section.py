# Code snippet to run clustering per section


# split adata per section to a dict.
tmp = adata.raw.to_adata()
alldata = {}
for batch in sections:
    alldata[batch] = tmp[tmp.obs['library_id'] == batch,]


# function for analysis:

def analyse_section(data):
    sc.pp.highly_variable_genes(data, flavor="seurat", n_top_genes=1500, inplace=True)
    data = data[:, data.var['highly_variable']]
    sc.pp.scale(data, max_value=10)
    sc.tl.pca(data, svd_solver='arpack')
    sc.pp.neighbors(data)
    sc.tl.umap(data)
    sc.tl.leiden(data, key_added="clusters_section")
    return data

# run per section
for section in sections:
    alldata[section] = analyse_section(alldata[section])


# function for plotting
def plot_per_section(data, section):
    fig, axs = plt.subplots(2, 3, figsize=(18,10), constrained_layout=True)
    sc.pl.umap(data,  color='clusters', palette=sc.pl.palettes.default_20, show=False, ax=axs[0,0])
    sc.pl.spatial(data, color='clusters',
        img_key="hires",
        library_id=section,
        size=1.5,
        palette=[
            v
            for k, v in clusters_colors.items()
            if k in ad.obs.clusters.unique().tolist()
        ],
        legend_loc=None,
        show=False,
        ax=axs[0,1],
    )

    tmp = pd.crosstab(data.obs['clusters_section'],data.obs['clusters'], normalize='index')
    tmp.plot.bar(stacked=True, ax=axs[0,2], legend = False)

    sc.pl.umap(data,  color='clusters_section', palette=sc.pl.palettes.default_20, show=False, ax=axs[1,0])
    sc.pl.spatial(data, color='clusters_section',
        img_key="hires",
        library_id=section,
        size=1.5,
        palette=[
            v
            for k, v in clusters_colors.items()
            if k in ad.obs.clusters.unique().tolist()
        ],
        legend_loc=None,
        show=False,
        ax=axs[1,1],
    )

    tmp = pd.crosstab(data.obs['clusters'],data.obs['clusters_section'], normalize='index')
    tmp.plot.bar(stacked=True, ax=axs[1,2], legend = False)


# plot

for section in sections:
    plot_per_section(alldata[section],section)    
