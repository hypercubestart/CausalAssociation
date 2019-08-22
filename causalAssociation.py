# algorithm to identify causal associations between mutations and regulons
# ie. find relationships where a mutation cases a change in regulon activity
# given the variants and gene expression of patients (derived using scRNA-seq)

# Genes and variants are used interchangeably because we decided to reduce all variants on a single gene to the gene
# to reduce the sample size

# import dependencies
import numpy as np
import random
import datetime
from scipy import stats
import math
from single_cell import SingleCell
import collections
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import label_binarize
import matplotlib.pyplot as plt
import pandas as pd
import dill
import time
import multiprocessing
import os
import operator
from collections import Counter
import scikitplot as skplt

# import statistical tests
from t_test import t_test
from mann_whitney_u_test import mann_whitney_u_test, mann_whitney_u_test_multiprocessing

def printt(message):
    """Print message with timestamp
        :param message: string
    """
    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t {}".format(message)))
    return None

def calculatePercentage(actual_mapping, predicted_mapping, name, file):
    """Calculate percentage of correct associations
        :param actual_mapping: dictionary {variant_index -> [regulon_index, effect]
        :param predicted_mapping: dictionary {variant_index -> [regulon_index, effect]
        :param name: string of statistical test
        :param file: File object to write out
        :return set of variant_indices of incorrect associations
    """
    number_correct = 0
    correct_effects = 0
    number_wrong = 0
    wrong_variant_regulon_association = set()
    for key, value in predicted_mapping.items():
        if key in actual_mapping and value[0] == actual_mapping[key][0]:
            number_correct+=1
            if value[1] == actual_mapping[key][1]:
                correct_effects+=1
        else:
            number_wrong+=1
            wrong_variant_regulon_association.add(key)

    printt("{}: \n\t correct: {}"
           "\n\t\t\t correct effects: {}"
           "\n\t wrong: {}"
           "\n\t ignore: {} "
           "\n\t total: {}\n".
           format(name, number_correct/len(actual_mapping),
                  np.float64(correct_effects)/number_correct,
                  number_wrong/len(actual_mapping),
                  (len(actual_mapping) - number_wrong - number_correct)/len(actual_mapping),
                  len(actual_mapping)))

    file.write("{}: \n\t correct: {}"
           "\n\t\t\t correct effects: {}"
           "\n\t wrong: {}"
           "\n\t ignore: {} "
           "\n\t total: {}\n".
           format(name, number_correct/len(actual_mapping),
                  np.float64(correct_effects)/number_correct,
                  number_wrong/len(actual_mapping),
                  (len(actual_mapping) - number_wrong - number_correct)/len(actual_mapping),
                  len(actual_mapping)))
    return wrong_variant_regulon_association

def auc_score(actual_mapping, predicted_mapping, p_values, name, regulon_count, gene_count, file):
    """Calculate AUC (area under the curve) score
        :param actual_mapping: dictionary {variant_index -> [regulon_index, effect]}
        :param predicted_mapping: dictionary {variant_index -> [regulon_index, effect]}
        :param p_values: matrix of p-values [variants * regulons]
        :param name: string of statistical test
        :param regulon_count: integer number of regulons
        :param gene_count: integer number of variants (assume each variant specific to single gene)
        :param file: File object to write out
        :return auc score
    """
    y_true = [] # true classes
    filteredp_values = [] #probabilities

    for mutationIndex in range(gene_count):
        if (mutationIndex in predicted_mapping):
            y_true.append(actual_mapping[mutationIndex][0])
            filteredp_values.append(p_values[mutationIndex])

    y_true = label_binarize(y_true, classes=range(regulon_count))

    probabilities = calculateProbabilites(filteredp_values)
    y_true, probabilities = removeColumnsWithAllZeros(y_true, probabilities)



    auc_score = roc_auc_score(y_true, probabilities)
    printt("{}: roc auc score for classifying the most associated regulon {}".format(name, auc_score))
    plot_multi_class_roc_curve(y_true, probabilities)
    file.write('roc auc score: {}'.format(auc_score))
    return auc_score


def calculateProbabilites(filteredp_values):
    """Calculate probabilities using p-values
        :param filteredp_values: matrix of p-values [variants * regulons]
        :return DataFrame of probabilities [variants * regulons] where df[i, j] is probability variant i is associated with regulon j
    """
    filteredp_values = pd.DataFrame(filteredp_values)
    sums = len(filteredp_values.columns) - filteredp_values.sum(axis = 1)
    sums.map(lambda x: len(filteredp_values.columns) - x)
    probabilities = filteredp_values.apply(lambda x: (1 - x) / sums[x.name], axis = 1)
    return probabilities

def removeColumnsWithAllZeros(y_true, df):
    """Filter out columns in binary encoding of classes where entire column is 0
    :param y_true: list-like result of transformation for fixed set of labels into 0s and 1s
    :param df: DataFrame to remove same columns as in y_true
    :return: DataFrame, DataFrame
    """
    # cannot run roc auc with empty class
    y_true_df = pd.DataFrame(y_true)

    for i in range(len(y_true[0])):
        for j in range(len(y_true)):
            if y_true[j][i] != 0:
                break
            if (j == len(y_true) - 1):
                del y_true_df[i]
                del df[i]
                i-=1
    return y_true_df, df

def plot_roc_curve(fpr, tpr):
    """Plot Roc Curve
    :param fpr: false positive rate
    :param tpr: true positive rate
    """
    #TODO: untested and unused
    plt.plot(fpr, tpr, color='orange', label='ROC')
    plt.plot([0, 1], [0, 1], color='darkblue', linestyle='--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic (ROC) Curve')
    plt.legend()
    plt.show()

def plot_multi_class_roc_curve(y_true, y_probas):
    """Plot Roc Curve
    :param fpr: true classes
    :param tpr: probability values
    """
    #TODO: untested and unused
    skplt.metrics.plot_roc_curve(y_true, y_probas)
    plt.show()

# def generateData(mutation_rate = 0.1, variant_count=10, sample_size = 1000, regulons_count = 200, mutation_not_found_rate = 0.2, noise=0.2):
#     """Previous function used to generate sample data of variant to regulon associations
#     :param mutation_rate: rate at which mutation occurs
#     :param variant_count: number of variants to have one-to-one association with regulon
#     :param sample_size: number of patients
#     :param regulons_count: number of regulons
#     :param mutation_not_found_rate: rate that mutation not found
#     :param noise: rate of random noise within regulon activity
#     :returns: cell_data [SingleCell] and dict {variant_index -> [regulon_index, effect]}
#     """
#     # len = # of mutations and number of patients
#     cell_data = []
#
#     # create one to one mutation to regulon mapping
#     availableRegulons = list(range(regulons_count))
#     dict = {} # variant index -> [regulon_index, effect]
#
#     for x in range(variant_count):
#         random_index = random.randint(0, len(availableRegulons) - 1)
#         random_effect = random.randint(-1, 1)
#
#         dict[x] = [availableRegulons[random_index], random_effect]
#         del availableRegulons[random_index]
#
#     # randomize cell data
#     for _ in range(sample_size):
#         cell = SingleCell(variant_count, regulons_count)
#         for i in range(variant_count):
#             if random.uniform(0,1) < mutation_rate:
#                 if random.uniform(0, 1) > mutation_not_found_rate:
#                     cell.set_variant(i, 1)
#                 cell.set_regulon(dict[i][0], dict[i][1])
#         cell_data.append(cell)
#
#     # add noise to regulon data
#     for cell in cell_data:
#         for i in range(len(cell.regulon_activity)):
#             if random.uniform(0,1) < noise:
#                 cell.set_regulon(i, random.randint(-1,1))
#
#     return cell_data, dict

def generateData(sample_size, gene_count, regulon_count, genes_mutated_count, genes_random_rate, samples_mutated_rate,
                 regulons_random_rate, miss_mutation_rate, miss_regulon_rate):
    """ Generate sample data of variant to regulon associations

    :param sample_size: number of patients
    :param gene_count: number of genes (or variants)
    :param regulon_count: number of regulons
    :param genes_mutated_count: number of genes with variants
    :param genes_random_rate: probability not mutated gene is observed as mutated
    :param samples_mutated_rate: percentage of samples with mutated genes
    :param regulons_random_rate: random distribution of regulon activity among non-affected regulons
    :param miss_mutation_rate: probability of there being a mutation but missing it
    :param miss_regulon_rate: probability that activity of associated regulon is not expected
    :returns: (sample_size * gene_count matrix of gene information, sample_size * regulon_count matrix of regulon activity, {variant_index -> [regulon_index, effect]})
    :rtype: (numpy.ndarray, numpy.ndarray, dictionary)
    """
    printt('starting to generate data...')
    start_time = time.time()

    gene_data = np.zeros((sample_size, gene_count))
    regulon_data = np.zeros((sample_size, regulon_count))

    # get genes_mutated_count genes to set as mutated genes and create associations
    availableGenes = list(range(gene_count)) # keep track of genes without associations
    mutationAssociations = {} # {geneIndex -> [regulonIndex, regulonEffect]}
    availableRegulons = list(range(regulon_count)) # keep track of regulons not associated with gene mutation
    for _ in range(genes_mutated_count):
        random_gene_index = random.randint(0, len(availableGenes) - 1)
        random_regulon_index = random.randint(0, len(availableRegulons) - 1)

        # generate -1 or 1 with equal probability
        random_effect = random.randint(0, 1)
        if random_effect == 0:
            random_effect = -1

        mutationAssociations[availableGenes[random_gene_index]] = [availableRegulons[random_regulon_index], random_effect]
        del availableGenes[random_gene_index]
        del availableRegulons[random_regulon_index]

    # create random distribution of noise among non_mutated genes
    for i in range(sample_size):
        for j in range(gene_count):
            if j not in mutationAssociations and random.uniform(0, 1) < genes_random_rate:
                gene_data[i][j] = 1

    # create random distribution of noise among non_associated regulons
    availableRegulons = set(availableRegulons)
    for i in range(sample_size):
        for j in range(regulon_count):
            if j not in availableRegulons and random.uniform(0, 1) < regulons_random_rate:
                # generate -1 or 1 with equal probability
                random_effect = random.randint(0, 1)
                if random_effect == 0:
                    random_effect = -1
                regulon_data[i][j] = random_effect

    # get random samples as samples with the mutation
    samples_with_mutation = int(sample_size * samples_mutated_rate)
    # for each mutated gene, get subset of samples
    for key in mutationAssociations:
        # inject causal association into random samples (number = samples_with_mutation)
        for _ in range(samples_with_mutation):
            availableSamples = list(range(sample_size))
            random_sample_index = random.randint(0, len(availableSamples) - 1)

            # inject causal association
            random_sample = availableSamples[random_sample_index]
            val = mutationAssociations[key]
            if random.uniform(0, 1) > miss_mutation_rate: # random distribution of not found mutation
                gene_data[random_sample][key] = 1
            regulon_data[random_sample][val[0]] = val[1]
            if random.uniform(0, 1) < miss_regulon_rate:
                regulon_data[random_sample][val[0]] -= val[1] * random.randint(1, 2) # random distribution of associated regulon activity not being the expected

            del availableSamples[random_sample_index]

    printt('finished generating data in {:.2f} minutes'.format((time.time() - start_time)/60.))
    return gene_data, regulon_data, mutationAssociations

def get_predicted_mapping(p_values, gene_data, regulon_data, alpha):
    """Use p-values to determine statistically significant associations between genes and regulons
    :param p_values: matrix of p-values [variants * regulons]
    :param gene_data: sample_size * gene_count matrix of gene data, 0 = non-mutated, 1 = mutated
    :param regulon_data: sample_size * regulon_count matrix of regulon activity, -1 = downregulated, 0 = normal, 1 = upregulated
    :param alpha: significance level
    :return: dict {variant_index -> [regulon_index, effect]}
    """
    predicted_mapping = {}
    gene_count = len(p_values)
    sample_size = len(regulon_data)
    for gene_index in range(gene_count):
        row = p_values[gene_index]
        min_val = min(row)

        if min_val < alpha:
            min_index = np.where(row == min_val)[0][0]

            regulon_activity_mutated = []
            regulon_activity_normal = []
            for i in range(sample_size):
                if gene_data[i][gene_index] == 1:
                    regulon_activity_mutated.append(regulon_data[i][min_index])
                else:
                    regulon_activity_normal.append(regulon_data[i][min_index])
            counts = Counter(regulon_activity_mutated)
            sorted_counts = sorted(counts.items(), key=operator.itemgetter(1), reverse=True)
            mode_normal_regulon = stats.mode(regulon_activity_normal, axis=None).mode[0]

            index = 0
            if (sorted_counts[0][0] == mode_normal_regulon):
                index = 1
            predicted_mapping[gene_index] = [min_index, sorted_counts[index][0]]

    return predicted_mapping

def run_tests():
    """Run series of tests to test robustness of algorithm"""
    sample_size = 300#3000
    gene_count = 10#10000
    regulon_count = 100#1000
    genes_mutated_count = 5#100
    samples_mutated_rate = [0.05] # percentage of samples with mutated genes, we expect 0.05-0.15
    genes_random_rate = [0.05] # probability not mutated gene is observed as mutated 0.05
    regulons_random_rate = [0.1] # random distribution of regulon activity among non-affected regulons 0.1
    miss_mutation_rate = [0.1]#[0.7, 0.95] # probability of there being a mutation but missing it 0.1 - 0.5
    miss_regulon_rate = [0.15] # probability that activity of associated regulon is not expected 0.05 - 0.15

    for i in samples_mutated_rate:
        for j in genes_random_rate:
            for k in regulons_random_rate:
                for l in miss_mutation_rate:
                    for m in miss_regulon_rate:
                        causal_association(sample_size, gene_count, regulon_count, genes_mutated_count, i, j, k, l, m)
                        print('\n')

def causal_association(sample_size, gene_count, regulon_count, genes_mutated_count, samples_mutated_rate,
                       genes_random_rate, regulons_random_rate, miss_mutation_rate, miss_regulon_rate):
    """Run Causal Association algorithm
    :param sample_size: number of patients
    :param gene_count: number of genes (or variants)
    :param regulon_count: number of regulons
    :param genes_mutated_count: number of genes with variants
    :param genes_random_rate: probability not mutated gene is observed as mutated
    :param samples_mutated_rate: percentage of samples with mutated genes
    :param regulons_random_rate: random distribution of regulon activity among non-affected regulons
    :param miss_mutation_rate: probability of there being a mutation but missing it
    :param miss_regulon_rate: probability that activity of associated regulon is not expected
    :return: None
    """
    start_time_timer = time.time()
    start_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # generate cell data
    test = "mann whitney u test"

    output_folder = './output/'

    gene_data, regulon_data, mutation_associations = generateData(sample_size=sample_size, gene_count=gene_count,
                                             regulon_count=regulon_count, genes_mutated_count=genes_mutated_count,
                                             genes_random_rate=genes_random_rate, samples_mutated_rate=samples_mutated_rate,
                                             regulons_random_rate=regulons_random_rate, miss_mutation_rate=miss_mutation_rate,
                                                                  miss_regulon_rate=miss_regulon_rate)
    # Bonferroni correction method for statistical test using multiple comparisons
    alpha = 0.05 / (gene_count * regulon_count)
    num_cores = multiprocessing.cpu_count()
    pvalues = mann_whitney_u_test_multiprocessing(gene_data, regulon_data, alpha, num_cores=num_cores)
    predicted_mapping = get_predicted_mapping(pvalues, gene_data, regulon_data, alpha)

    end_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    file_name = (datetime.datetime.now().strftime("%Y-%m-%d-%H:%M:%S"))

    # double check that os path exists
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)
    # write output file
    with open(output_folder + file_name + '.txt', 'w') as file:
        file.write('program started: {}\n'.format(start_time))
        file.write('program ended: {}\n'.format(end_time))
        file.write('time elapsed: {:.2f} minutes\n'.format((time.time() - start_time_timer)/60))
        file.write('test: {}\n\n'.format(test))

        file.write('Parameters: \n')
        file.write('sample size: {}\n'.format(sample_size))
        file.write('gene count: {}\n'.format(gene_count))
        file.write('regulon count: {}\n'.format(regulon_count))
        file.write('mutated genes count: {} ({:.2f}%)\n'.format(genes_mutated_count, genes_mutated_count/gene_count))
        file.write('mutated samples count: {} ({:.2f}%)\n'.format(int(samples_mutated_rate * sample_size), samples_mutated_rate))
        file.write('\n')
        file.write('genes random rate: {}\n'.format(genes_random_rate))
        file.write('regulons random rate: {}\n'.format(regulons_random_rate))
        file.write('miss mutation rate: {}\n'.format(miss_mutation_rate))
        file.write('miss regulon rate: {}\n'.format(miss_regulon_rate))

        try:
            mwu_test_incorrect = calculatePercentage(mutation_associations, predicted_mapping, test, file)
            auc_score(mutation_associations, predicted_mapping, pvalues, test, regulon_count, gene_count, file)
            dill.dump_session(output_folder + file_name + '.dill')
        except:
            dill.dump_session(output_folder + file_name + '.dill')
    return None
if __name__ == "__main__":
    run_tests()
    # start_time_timer = time.time()
    # start_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # # generate and save cell data
    #
    # sample_size = 3000
    # gene_count = 10000
    # regulon_count = 1000
    # genes_mutated_count = 100
    # samples_mutated_rate = 0.2 # percentage of samples with mutated genes
    # genes_random_rate = 0.01 # probability not mutated gene is observed as mutated
    # regulons_random_rate = 0.01 # random distribution of regulon activity among non-affected regulons
    # miss_mutation_rate = 0.1 # probability of there being a mutation but missing it
    # miss_regulon_rate = 0.1 # probability that activity of associated regulon is not expected
    # test = "mann whitney u test"
    #
    # output_folder = './output/'
    #
    # gene_data, regulon_data, mutation_associations = generateData(sample_size=sample_size, gene_count=gene_count,
    #                                          regulon_count=regulon_count, genes_mutated_count=genes_mutated_count,
    #                                          genes_random_rate=genes_random_rate, samples_mutated_rate=samples_mutated_rate,
    #                                          regulons_random_rate=regulons_random_rate, miss_mutation_rate=miss_mutation_rate,
    #                                                               miss_regulon_rate=miss_regulon_rate)
    # dill.dump_session('./cell_data.dill')
    #
    # # perform mann_whitney u test
    # dill.load_session('./cell_data.dill')
    # alpha = 0.05 / (gene_count * regulon_count)
    # num_cores = multiprocessing.cpu_count()
    # pvalues = mann_whitney_u_test_multiprocessing(gene_data, regulon_data, alpha, num_cores=num_cores)
    # predicted_mapping = get_predicted_mapping(pvalues, gene_data, regulon_data, alpha)
    #
    # end_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    #
    # file_name = (datetime.datetime.now().strftime("%Y-%m-%d-%H:%M:%S"))
    #
    # # double check that os path exists
    # if not os.path.isdir(output_folder):
    #     os.mkdir(output_folder)
    # # write output file
    # with open(output_folder + file_name + '.txt', 'w') as file:
    #     file.write('program started: {}\n'.format(start_time))
    #     file.write('program ended: {}\n'.format(end_time))
    #     file.write('time elapsed: {:.2f} minutes\n'.format((time.time() - start_time_timer)/60))
    #     file.write('test: {}\n\n'.format(test))
    #
    #     file.write('Parameters: \n')
    #     file.write('sample size: {}\n'.format(sample_size))
    #     file.write('gene count: {}\n'.format(gene_count))
    #     file.write('regulon count: {}\n'.format(regulon_count))
    #     file.write('mutated genes count: {} ({:.2f}%)\n'.format(genes_mutated_count, genes_mutated_count/gene_count))
    #     file.write('mutated samples count: {} ({:.2f}%)\n'.format(int(samples_mutated_rate * sample_size), samples_mutated_rate))
    #     file.write('\n')
    #     file.write('genes random rate: {}\n'.format(genes_random_rate))
    #     file.write('regulons random rate: {}\n'.format(regulons_random_rate))
    #     file.write('miss mutation rate: {}\n'.format(miss_mutation_rate))
    #     file.write('miss regulon rate: {}\n'.format(miss_regulon_rate))
    #
    #     mwu_test_incorrect = calculatePercentage(mutation_associations, predicted_mapping, test, file)
    #     auc_score(mutation_associations, predicted_mapping, pvalues, test, regulon_count, gene_count, file)
    # # mwu_test_incorrect = calculatePercentage(mutation_associations, predicted_mapping, "mann whitney u test")
    # # auc_score(mutation_associations, predicted_mapping, p_values, "mann whitney u test", regulon_count, gene_count)
    #
    # # perform t test //TODO: old and not updated
    # # t_test_predicted, t_test_p_values = t_test(cell_data, variant_count)
    # # t_test_incorrect = calculatePercentage(actual_mapping, t_test_predicted, "t-test")
    # # auc_score(actual_mapping, t_test_predicted, t_test_p_values, "t-test", regulon_count, variant_count)
    #
    # dill.dump_session('./debug.dill')










def chisquare(cell_data):
    """Use Chi-Square to test for causal association
    :param cell_data: List of SingleCell
    """
    #TODO: frequency of each variation has to be > 5, so oftentimes chi-square doens't work
    threshold = 0.10
    contingency_table_min_freq = 1
    predicted_mapping = {}
    for mutationIndex in range(SingleCell.variant_count):
        # get indices of cells with mutation
        with_mutation = [cell for cell in cell_data if cell.get_variant(mutationIndex) == 1]
        without_mutation = [cell for cell in cell_data if cell.get_variant(mutationIndex) == 0]

        if len(with_mutation) == 0 or len(without_mutation) == 0:
            printt("Chisquare: Mutation at index {} cannot be associated because len == 0".format(mutationIndex))
            continue

        p_values = []
        for i in range(SingleCell.regulon_count):
            contingency_table = np.zeros((2,3))

            # get regulon activities for current regulon
            regulon_activity_with_mutation = [cell.get_regulon(i) for cell in with_mutation]
            regulon_activity_without_mutation = [cell.get_regulon(i) for cell in without_mutation]

            # get regulon activity counts
            mutation_count = collections.Counter(regulon_activity_with_mutation)
            no_mutation_count = collections.Counter(regulon_activity_without_mutation)

            if mutation_count[-1] < contingency_table_min_freq or mutation_count[0] < contingency_table_min_freq or \
                    mutation_count[1] < contingency_table_min_freq or no_mutation_count[-1] < contingency_table_min_freq or no_mutation_count[0] < contingency_table_min_freq or no_mutation_count[1] < contingency_table_min_freq:
                printt("Chisquare: Mutation at index {} cannot be associated because less than minimum contingency table frequency".format(mutationIndex))
                p_values.append(1)
                continue

            # TODO: normalize chi-square values (eq. percentage)
            contingency_table[0][0] = mutation_count[-1]
            contingency_table[0][1] = mutation_count[0]
            contingency_table[0][2] = mutation_count[1]

            contingency_table[1][0] = no_mutation_count[-1]
            contingency_table[1][1] = no_mutation_count[0]
            contingency_table[1][2] = no_mutation_count[1]

            p_values.append(stats.chi2_contingency(contingency_table))

        # ignore if mutation and no-mutation are not associated
        if (min(p_values) > threshold):
            continue

        # get indices with lowest p value
        min_p_value = min(p_values)
        index = p_values.index(min_p_value)
        # use most common regulon activity to determine effect of variant on regulon
        effect = collections.Counter([cell.get_regulon(index) for cell in with_mutation]).most_common()[0][0]
        predicted_mapping[mutationIndex] = [index, effect, min_p_value]
    return predicted_mapping
