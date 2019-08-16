import sys
import miner2 
import miner2.preprocess
import miner2.mechanistic_inference
import miner2.miner
import pandas as pd
import os
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import dill
import pickle
import seaborn as sb
import time

tf_2_genes_path = '../miner2/miner2/data/tfbsdb_tf_to_genes.pkl' # location of tfbs_db
results_dir='./results/GSM3587977_AML707B/'

# dill.dump_session(results_dir+'info/bottle.dill')
dill.load_session(results_dir+'info/bottle.dill')

# Infer transcriptional programs
reference_df = eigengenes.copy()
programs, states_ = miner2.miner.mosaic(dfr=eigengenes.copy(),clusterList=centroidClusters,minClusterSize_x=int(np.ceil(0.01*expressionData.shape[1])),minClusterSize_y=5,allow_singletons=False,max_groups=50,saveFile=os.path.join(results_dir,"regulon_activity_heatmap.pdf"),random_state=12)
transcriptional_programs, program_regulons = miner2.miner.transcriptionalPrograms(programs,referenceDictionary)
program_list = [program_regulons[("").join(["TP",str(i)])] for i in range(len(program_regulons))]
mosaicDf = reference_df.loc[np.hstack(program_list),np.hstack(states_)]
mosaicDf.to_csv(os.path.join(results_dir,"regulons_activity_heatmap.csv"))


