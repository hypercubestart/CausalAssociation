import pickle,numpy,pandas
import matplotlib,matplotlib.pyplot
import scipy,scipy.stats
import scanpy
scanpy.settings.verbosity=5

# 0. user defined variables
adata=scanpy.read_10x_mtx('/Volumes/omics4tb2/alomana/projects/i18/results/deconvolution/cell_ranger/both/aggregated_both/outs/filtered_feature_bc_matrix',var_names='gene_symbols',cache=True)
jarFile='/Volumes/omics4tb2/alomana/projects/i18/results/deconvolution/scanpy/species.cellIDs.run2.001.pickle'

# 1. read data files
adata.var_names_make_unique() 
print(adata)

# 2. preprocess
scanpy.pp.filter_cells(adata,min_genes=200)
scanpy.pp.filter_genes(adata,min_cells=3)

# 2.1. Retrieve mouse cells only
cellIDs=adata.obs_names.tolist()
geneNames=adata.var_names.tolist()

f=open(jarFile,'rb')
[mouseCellIDs,humanCellIDs,chimericCellIDs]=pickle.load(f)
f.close()

print(len(cellIDs),len(geneNames))
print('mouse',len(mouseCellIDs))
print('human',len(humanCellIDs))
print('chimeric',len(chimericCellIDs))

# slice in mouse cells
print('before slicing mouse cells...')
print(adata)
print()
adata=adata[mouseCellIDs,:]
print('after')
print(adata)

# slice in mouse genes
mouseGenes=[element for element in geneNames if element[:4] == 'mm10']
print()
print('before slicing mouse genes...')
print(adata)
print()
adata=adata[:,mouseGenes]
print('after')
print(adata)

# 2.2. highly abundant transcripts
scanpy.pl.highest_expr_genes(adata,n_top=20,show=False,save='.pdf')

# 2.3. QC based on mitochondrial genes and number of counts
mitoGenes=[element for element in adata.var_names if element[:8] == 'mm10_mt-']
print(mitoGenes)

adata.obs['percent_mito']=numpy.sum(adata[:,mitoGenes].X,axis=1).A1/numpy.sum(adata.X,axis=1).A1
adata.obs['n_counts'] = adata.X.sum(axis=1).A1
scanpy.pl.violin(adata,['n_genes','n_counts','percent_mito'],jitter=0.4,multi_panel=True,show=False,save='.pdf')

print(adata)
scanpy.pl.scatter(adata, x='n_counts', y='percent_mito',alpha=1/4,show=False,save='.before.pdf')
scanpy.pl.scatter(adata, x='n_counts', y='n_genes',alpha=1/4,show=False,save='.before.pdf')

print('remove high mito')
adata = adata[adata.obs['percent_mito'] < 0.05, :]
print(adata)
print()

print('remove top outliers')
adata = adata[adata.obs['n_counts'] < 18000, :]
print(adata)
print()

print('remove low n_count')
adata = adata[adata.obs['n_counts'] > 1500, :]
print(adata)
print()

print('remove low n_genes')
adata = adata[adata.obs['n_genes'] > 750, :]
print(adata)
print()

scanpy.pl.scatter(adata, x='n_counts', y='percent_mito',alpha=1/4,show=False,save='.after.pdf')
scanpy.pl.scatter(adata, x='n_counts', y='n_genes',alpha=1/4,show=False,save='.after.pdf')

# 2.4. normalization and transformation
scanpy.pp.normalize_per_cell(adata, counts_per_cell_after=18e3)
scanpy.pp.log1p(adata)
adata.raw = adata

scanpy.pp.highly_variable_genes(adata,min_mean=0.0125,max_mean=3,min_disp=0.5)
scanpy.pl.highly_variable_genes(adata,show=False,save='.pdf')

print(adata)
adata = adata[:,adata.var['highly_variable']]
print(adata)

# 2.5. Regress out effects of total counts per cell and the percentage of mitochondrial genes
scanpy.pp.regress_out(adata,['n_counts','percent_mito'])
scanpy.pp.scale(adata,max_value=10)

# 2.6. Associate treatment to cell ID
cellIDs=adata.obs_names.tolist()
print(len(cellIDs))

cellConditions=[]
for cellID in cellIDs:
    if '-1' in cellID: 
        cellConditions.append('WT')
    elif '-2' in cellID:
        cellConditions.append('MAVS')
    elif '-3' in cellID:
        cellConditions.append('NLRP3')
    else:
        raise ValueError('cellID not recognized')
    
print(cellConditions,len(cellConditions))
print('number of WT cells:',cellConditions.count('WT'))
print('number of MAVS cells:',cellConditions.count('MAVS'))
print('number of NLRP3 cells:',cellConditions.count('NLRP3'))

sys.exit()

# 3. visualization
components=50

# 3.1. PCA
scanpy.tl.pca(adata, svd_solver='arpack',n_comps=components)
scanpy.pl.pca_variance_ratio(adata, log=True,n_pcs=components,show=False,save='.pdf')

adata.obs['genotype']=cellConditions

scanpy.pl.pca(adata,color='genotype',palette=['black','red','blue'],alpha=0.5,show=False,save='.genotype.pdf')
scanpy.pl.pca(adata,color='n_counts',palette='viridis',show=False,save='.counts.pdf')
scanpy.pl.pca(adata,color='n_genes',palette='viridis',show=False,save='.genes.pdf')

# 3.2. tSNE
scanpy.tl.tsne(adata)
scanpy.pl.tsne(adata,color='genotype',palette=['black','red','blue'],alpha=1/2,show=False,save='.genotype.pdf')
scanpy.pl.tsne(adata,color='n_counts',palette='viridis',show=False,save='.counts.pdf')
scanpy.pl.tsne(adata,color='n_genes',palette='viridis',show=False,save='.genes.pdf')

# 3.3. UMAP
print('UMAP resolution...')
possibleNeighbors=numpy.arange(5,50+1,1)
louvainClusters=[]
print(possibleNeighbors)
for nei in possibleNeighbors:
    scanpy.pp.neighbors(adata,n_neighbors=nei,n_pcs=components,knn=True)

    scanpy.tl.umap(adata)
    scanpy.pl.umap(adata, color='genotype',palette=['black','red','blue'],alpha=1/2,show=False,save='.genotype.nei.{}.pdf'.format(nei))
    scanpy.pl.umap(adata, color='n_counts',palette='viridis',alpha=1/2,show=False,save='.counts.nei.{}.pdf'.format(nei))
    scanpy.pl.umap(adata, color='n_genes',palette='viridis',alpha=1/2,show=False,save='.genes.nei.{}.pdf'.format(nei))

    scanpy.tl.louvain(adata)
    scanpy.pl.umap(adata, color='louvain',palette='tab10',alpha=1/2,show=False,save='.louvain.nei.{}.pdf'.format(nei))

    uniqueClusters=adata.obs['louvain'].nunique()
    louvainClusters.append(uniqueClusters)

matplotlib.pyplot.plot(possibleNeighbors,louvainClusters,'o-',color='black')
matplotlib.pyplot.xlabel('neighbors')
matplotlib.pyplot.ylabel('clusters')
matplotlib.pyplot.savefig('louvainRanks.pdf')
