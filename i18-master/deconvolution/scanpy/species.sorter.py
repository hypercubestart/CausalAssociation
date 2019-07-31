import numpy,sys,pickle,datetime
import multiprocessing,multiprocessing.pool

import scanpy
scanpy.settings.verbosity=5

def speciesRatioComputer(cellID):
    
    '''
    This function computes the ratio between mouse and human counts.
    '''
    
    sumHsa=0; sumMmu=0
    for geneName in geneNames:
        countValue=adata[cellID,geneName].X
        if geneName[:4] == 'hg19':
            sumHsa=sumHsa+countValue
        elif geneName[:4] == 'mm10':
            sumMmu=sumMmu+countValue
        else:
            raise ValueError("new species")
    ratio=sumMmu/(sumHsa+sumMmu)
    if ratio > 0.95:
        label='mouse'
        print('\t {} \t {} \t {} \t {:.3f} \t {:0.0f} \t {:0.0f}'.format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),label,cellID,ratio,sumHsa,sumMmu))
    elif ratio < 0.05:
        label='human'
        print('\t {} \t {} \t {} \t {:.3f} \t {:0.0f} \t {:0.0f}'.format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),label,cellID,ratio,sumHsa,sumMmu))
    else:
        label='chimeric'
        print('\t {} \t {} \t {} \t {:.3f} \t {:0.0f} \t {:0.0f}'.format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),label,cellID,ratio,sumHsa,sumMmu))
    
    return label

###
### MAIN
###

# 0. user defined variables
matrixFolder='/Volumes/omics4tb2/alomana/projects/i18/results/deconvolution/cell_ranger/both/aggregated_both/outs/filtered_feature_bc_matrix'

jarFile='/Volumes/omics4tb2/alomana/projects/i18/results/deconvolution/scanpy/species.cellIDs.deconvolution.002.pickle'
numberOfCPUs=8 
numberOfCells=3000 # we currently have 2,340

# 1. Read data files
print('Reading files...')
adata=scanpy.read_10x_mtx(matrixFolder,var_names='gene_symbols',cache=True)
print('\t About to make them unique...')
adata.var_names_make_unique()
print('\t unique done.')
print(adata)
print()

# 2. Preprocess
print('Preprocessing...')

# 2.1. Simple cleaning
print('\t Simple cleaning...')
scanpy.pp.filter_cells(adata,min_genes=200)
scanpy.pp.filter_genes(adata,min_cells=3)

# 2.2. Retrieve mouse cells only
print('\t Retrieving mouse cells...')

cellIDs=adata.obs_names[:numberOfCells].tolist()
geneNames=adata.var_names.tolist()
print('working number of cells',len(cellIDs))
print('working number of genes',len(geneNames))

hydra=multiprocessing.pool.Pool(numberOfCPUs)
labels=hydra.map(speciesRatioComputer,cellIDs)
hydra.close()

print('\t Cell labels computed!')

# 2.3. Storing
print('\t Storing...')
humanCellIDs=[];mouseCellIDs=[];chimericCellIDs=[]
for i in range(len(labels)):
    if labels[i] == 'mouse':
        mouseCellIDs.append(cellIDs[i])
    if labels[i] == 'human':
        humanCellIDs.append(cellIDs[i])
    if labels[i] == 'chimeric':
        chimericCellIDs.append(cellIDs[i])

print(len(cellIDs),len(labels))
print('mouse',len(mouseCellIDs))
print('human',len(humanCellIDs))
print('chimeric',len(chimericCellIDs))

f=open(jarFile,'wb')
pickle.dump([mouseCellIDs,humanCellIDs,chimericCellIDs],f)
f.close()
    
# recover labels
print()
print('recovering cell types...')
f=open(jarFile,'rb')
[mouseCellIDs,humanCellIDs,chimericCellIDs]=pickle.load(f)
f.close()

print(len(cellIDs),len(labels))
print('mouse',len(mouseCellIDs))
print('human',len(humanCellIDs))
print('chimeric',len(chimericCellIDs))
