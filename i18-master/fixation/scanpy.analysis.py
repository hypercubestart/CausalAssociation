#!/usr/bin/env python
# coding: utf-8

# # Script to work with fixed macrophages

# In[1]:


### conda install seaborn scikit-learn statsmodels numba pytables
### conda install -c conda-forge python-igraph louvain
### pip install scanpy

import pickle,numpy
import matplotlib,matplotlib.pyplot
import scanpy
scanpy.settings.verbosity=5


# ## 1. Read data files

# In[2]:


### define input files
adata=scanpy.read_10x_mtx('/Volumes/omics4tb2/alomana/projects/i18/results/both_aggregated/outs/filtered_feature_bc_matrix',var_names='gene_symbols',cache=True)
adata.var_names_make_unique() 
adata


# ## 2. Preprocess

# In[3]:


scanpy.pp.filter_cells(adata,min_genes=200)
scanpy.pp.filter_genes(adata,min_cells=3)


# ### 2.1. Retrieve mouse cells only

# In[4]:


cellIDs=adata.obs_names.tolist()
geneNames=adata.var_names.tolist()

jarFile='species.cellIDs.run.006.pickle'
f=open(jarFile,'rb')
[mouseCellIDs,humanCellIDs,chimericCellIDs]=pickle.load(f)
f.close()

print(len(cellIDs),len(geneNames))
print('mouse',len(mouseCellIDs))
print('human',len(humanCellIDs))
print('chimeric',len(chimericCellIDs))


# In[5]:


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


# ### 2.2. Highly abundant genes

# In[6]:


scanpy.pl.highest_expr_genes(adata,n_top=20)


# ### 2.3. QC based on mitochondrial genes and number of counts

# In[7]:


mitoGenes=[element for element in adata.var_names if element[:8] == 'mm10_mt-']
print(mitoGenes)

adata.obs['percent_mito']=numpy.sum(adata[:,mitoGenes].X,axis=1).A1/numpy.sum(adata.X,axis=1).A1
adata.obs['n_counts'] = adata.X.sum(axis=1).A1
scanpy.pl.violin(adata,['n_genes','n_counts','percent_mito'],jitter=0.4,multi_panel=True)


# In[8]:


print(adata)
scanpy.pl.scatter(adata, x='n_counts', y='percent_mito')
scanpy.pl.scatter(adata, x='n_counts', y='n_genes')


# In[9]:


adata = adata[adata.obs['n_counts'] < 30000, :]
adata = adata[adata.obs['percent_mito'] < 0.05, :]
print(adata)
scanpy.pl.scatter(adata, x='n_counts', y='percent_mito')
scanpy.pl.scatter(adata, x='n_counts', y='n_genes')


# ### 2.4. Normalization and log transform

# In[10]:


scanpy.pp.normalize_per_cell(adata, counts_per_cell_after=30e3)


# In[11]:


scanpy.pp.log1p(adata)


# In[12]:


adata.raw = adata


# In[13]:


scanpy.pp.highly_variable_genes(adata,min_mean=0.0125,max_mean=4, min_disp=0.25)
scanpy.pl.highly_variable_genes(adata)


# In[14]:


adata = adata[:,adata.var['highly_variable']]


# ### 2.5. Regress out effects of total counts per cell and the percentage of mitochondrial genes  

# In[15]:


scanpy.pp.regress_out(adata,['n_counts', 'percent_mito'])


# In[16]:


scanpy.pp.scale(adata,max_value=10)


# ### 2.6. Associate treatment to cell ID

# In[17]:


f=open('fixedLabels.pickle','rb')
cellIDsFixed=pickle.load(f)
f.close()

f=open('freshLabels.pickle','rb')
cellIDsFresh=pickle.load(f)
f.close() 


# In[18]:


cellIDs=adata.obs_names.tolist()
print(len(cellIDs))

cellConditions=[]
for cellID in cellIDs:
    if '-1' in cellID: 
        cellConditions.append('fixed')
    elif '-2' in cellID:
        cellConditions.append('fresh')
    else:
        raise ValueError('cellID not recognized')
    
print(cellConditions,len(cellConditions))


# ## 3. Visualization

# ## 3.1. PCA

# In[20]:


scanpy.tl.pca(adata, svd_solver='arpack')


# In[21]:


adata.obs['treatment']=cellConditions
scanpy.pl.pca(adata,color='treatment',palette=['red','blue'],alpha=0.5)


# In[22]:


scanpy.pl.pca_variance_ratio(adata, log=True)


# In[23]:


adata.write('resultsFile.h5ad')


# ## 3.2. UMAP

# In[24]:


scanpy.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
scanpy.tl.umap(adata)
scanpy.pl.umap(adata, color=['mm10_Atp6v1h', 'mm10_Rb1cc1', 'mm10_St18'])


# In[25]:


scanpy.pl.umap(adata, color=['mm10_Atp6v1h', 'mm10_Rb1cc1', 'mm10_St18'], use_raw=False)


# ## 3.3 Louvain

# In[26]:


scanpy.tl.louvain(adata)
scanpy.pl.umap(adata, color=['louvain'])


# In[96]:


adata.write('resultsFile.h5ad')

