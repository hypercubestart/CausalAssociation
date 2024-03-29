# Use mann-whitney-u-test to test for causal associations
from scipy import stats
import numpy as np
import time
import datetime
import multiprocessing
import ctypes

def printt(message):
    """Print message with timestamp
        :param message: string
    """
    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t {}".format(message)))
    return None

def mann_whitney_u_test(gene_data, regulon_data, alpha):
    """Single threaded mann whitney u test
    :param gene_data: sample_size * gene_count matrix of gene data, 0 = non-mutated, 1 = mutated
    :param regulon_data: sample_size * regulon_count matrix of regulon activity, -1 = downregulated, 0 = normal, 1 = upregulated
    :param alpha: significence level
    :return: [gene_count * regulon_count] matrix of p-values where [i, j] is p-value of variant i and regulon j
    """
    # TODO: :param alpha currently not being used
    start_time = time.time()
    printt("starting mann-whitney-u-test...")
    sample_count = len(gene_data)
    gene_count = len(gene_data[0])
    regulon_count = len(regulon_data[0])

    pvalues = np.zeros((gene_count, regulon_count))

    for geneIndex in range(gene_count):
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
            for regulonIndex in range(regulon_count):
                pvalues[geneIndex][regulonIndex] = 1
            continue

        # perform t test for current gene index against every regulon index
        for regulonIndex in range(regulon_count):
            with_mutation_regulon_activity = []
            without_mutation_regulon_activity = []

            for i in range(len(with_mutation)):
                with_mutation_regulon_activity.append(with_mutation[i][regulonIndex])
            for i in range(len(without_mutation)):
                without_mutation_regulon_activity.append(without_mutation[i][regulonIndex])
            try:
                mwu_test = stats.mannwhitneyu(with_mutation_regulon_activity, without_mutation_regulon_activity,
                                              use_continuity=True,
                                              alternative='two-sided')  # perform mann-whitney u test
            except Exception:
                pvalues[geneIndex][regulonIndex] = 1
                continue

            # save all p-values into pvalues matrix
            pvalues[geneIndex][regulonIndex] = mwu_test.pvalue

    printt('finished mann_whitney_u_test in {:.2f} minutes'.format((time.time() - start_time) / 60.))
    return pvalues

shared_gene_data = None
shared_regulon_data = None

def init(gene_data_base, regulon_data_base, gene_data, regulon_data):
    """initiation function for multiprocessing, put gene_data and regulon_data into shared memory"""
    sample_count = len(gene_data)
    gene_count = len(gene_data[0])
    regulon_count = len(regulon_data[0])

    global shared_gene_data
    global shared_regulon_data
    shared_gene_data = np.ctypeslib.as_array(gene_data_base.get_obj())
    shared_gene_data = shared_gene_data.reshape(sample_count, gene_count)
    shared_regulon_data = np.ctypeslib.as_array(regulon_data_base.get_obj())
    shared_regulon_data = shared_regulon_data.reshape(sample_count, regulon_count)

    for i in range(sample_count):
        for j in range(gene_count):
            shared_gene_data[i][j] = gene_data[i][j]

    for i in range(sample_count):
        for j in range(regulon_count):
            shared_regulon_data[i][j] = regulon_data[i][j]


def mann_whitney_u_test_multiprocessing(gene_data, regulon_data, alpha, num_cores = 4):
    """Multi threaded mann whitney u test
    :param gene_data: sample_size * gene_count matrix of gene data, 0 = non-mutated, 1 = mutated
    :param regulon_data: sample_size * regulon_count matrix of regulon activity, -1 = downregulated, 0 = normal, 1 = upregulated
    :param alpha: significence level
    :param num_cores: number of cores
    :return: [gene_count * regulon_count] matrix of p-values where [i, j] is p-value of variant i and regulon j
    """
    start_time = time.time()
    printt("starting mann-whitney-u-test...")
    sample_count = len(gene_data)
    gene_count = len(gene_data[0])
    regulon_count = len(regulon_data[0])

    gene_data_base = multiprocessing.Array(ctypes.c_double, sample_count * gene_count)
    regulon_data_base = multiprocessing.Array(ctypes.c_double, sample_count * regulon_count)

    pool = multiprocessing.Pool(num_cores, initializer=init, initargs=(gene_data_base, regulon_data_base, gene_data, regulon_data))

    results = pool.map(calculateMWU, range(gene_count))

    printt('finished mann_whitney_u_test in {:.2f} minutes'.format((time.time() - start_time) / 60.))

    return results

def calculateMWU(geneIndex):
    """multithreaded mann whitney u test function
    :param geneIndex: calculate p-values
    :return list of p-values from comparison with every regulon
    """
    sample_count = len(shared_gene_data)
    regulon_count = len(shared_regulon_data[0])

    # split regulon activity by with variant and without
    with_mutation = []
    without_mutation = []
    pvalues = np.zeros(regulon_count)

    for sample_index in range(sample_count):
        if shared_gene_data[sample_index][geneIndex] == 1:
            with_mutation.append(shared_regulon_data[sample_index])
        else:
            without_mutation.append(shared_regulon_data[sample_index])

    # check that each group has more than one sample
    if len(without_mutation) <= 1 or len(with_mutation) <= 1:
        for regulonIndex in range(regulon_count):
            pvalues[regulonIndex] = 1
        return pvalues

    for regulonIndex in range(regulon_count):
        with_mutation_regulon_activity = []
        without_mutation_regulon_activity = []

        for i in range(len(with_mutation)):
            with_mutation_regulon_activity.append(with_mutation[i][regulonIndex])
        for i in range(len(without_mutation)):
            without_mutation_regulon_activity.append(without_mutation[i][regulonIndex])
        try:
            mwu_test = stats.mannwhitneyu(with_mutation_regulon_activity, without_mutation_regulon_activity,
                                          use_continuity=True,
                                          alternative='two-sided')  # perform mann-whitney u test
            pvalues[regulonIndex] = mwu_test.pvalue
        except Exception:
             pvalues[regulonIndex] = 1

    return pvalues
