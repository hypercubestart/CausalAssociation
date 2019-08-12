# program used to debug
import miner2.mechanistic_inference
import dill

results_dir='./results/GSM3587977_AML707B/'
dill.load_session(results_dir+'info/bottle.dill')

axes = miner2.mechanistic_inference.get_principal_df(revised_clusters,expression_data,subkey=None,min_number_genes=1)
mechanistic_output = miner2.mechanistic_inference.enrichment(axes,revised_clusters,expression_data,correlation_threshold=min_correlation,num_cores=num_cores)