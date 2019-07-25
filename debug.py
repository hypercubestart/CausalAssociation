import dill
from scipy import stats
import math
import collections
import datetime
import time
import numpy as np

def printt(message):
    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t {}".format(message)))
    return None

def mann_whitney_u_test_debug(gene_data, regulon_data, alpha, mutation_associations):
    start_time = time.time()
    printt("starting mann-whitney-u-test...")
    sample_count = len(gene_data)
    gene_count = len(gene_data[0])
    regulon_count = len(regulon_data[0])

    pvalues = np.zeros((gene_count, regulon_count))

    mutated_genes = sorted([i for i in mutation_associations])
    for geneIndex in mutated_genes:
        # split regulon activity by with variant and without
        with_mutation = []
        without_mutation = []

        for sample_index in range(sample_count):
            if gene_data[sample_index][geneIndex] == 1:
                with_mutation.append(regulon_data[sample_index])
            else:
                without_mutation.append(regulon_data[sample_index])

        # check that each group has more than one sample
        if len(without_mutation) <= 1 or len(with_mutation) <= 1:
            for regulonIndex in [mutation_associations[geneIndex][0]]:
                pvalues[geneIndex][regulonIndex] = 1
            continue

        for regulonIndex in [mutation_associations[geneIndex][0]]:
            with_mutation_regulon_activity = []
            without_mutation_regulon_activity = []

            for i in range(len(with_mutation)):
                with_mutation_regulon_activity.append(with_mutation[i][regulonIndex])
            for i in range(len(without_mutation)):
                without_mutation_regulon_activity.append(without_mutation[i][regulonIndex])
            try:
                mwu_test = stats.mannwhitneyu(with_mutation_regulon_activity, without_mutation_regulon_activity,
                use_continuity=True, alternative='two-sided')  # perform mann-whitney u test
            except Exception:
                pvalues[geneIndex][regulonIndex] = 1
                continue

            #if mwu_test.pvalue < alpha:
            pvalues[geneIndex][regulonIndex] = mwu_test.pvalue

    printt('finished mann_whitney_u_test in {:.2f} minutes'.format((time.time() - start_time) / 60.))
    return pvalues

if __name__ == "__main__":

    dill.load_session('./debug.dill')
    mann_whitney_u_test_debug(gene_data, regulon_data, alpha, mutation_associations)
    #t_test_debug(cell_data, t_test_incorrect)


