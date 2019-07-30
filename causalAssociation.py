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
from t_test import t_test
import time
import multiprocessing
import os
from mann_whitney_u_test import mann_whitney_u_test, mann_whitney_u_test_multiprocessing

def printt(message):
    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t {}".format(message)))
    return None

def calculatePercentage(actual_mapping, predicted_mapping, name, file):
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

    file.write('roc auc score: {}'.format(auc_score))
    return auc_score


def calculateProbabilites(filteredp_values):
    filteredp_values = pd.DataFrame(filteredp_values)
    sums = len(filteredp_values.columns) - filteredp_values.sum(axis = 1)
    sums.map(lambda x: len(filteredp_values.columns) - x)
    probabilities = filteredp_values.apply(lambda x: (1 - x) / sums[x.name],axis = 1)
    return probabilities

def removeColumnsWithAllZeros(y_true, df):
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
    plt.plot(fpr, tpr, color='orange', label='ROC')
    plt.plot([0, 1], [0, 1], color='darkblue', linestyle='--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic (ROC) Curve')
    plt.legend()
    plt.show()

# def generateData(mutation_rate = 0.1, variant_count=10, sample_size = 1000, regulons_count = 200, mutation_not_found_rate = 0.2, noise=0.2):
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

def generateData(sample_size, gene_count, regulon_count, genes_mutated_count, genes_random_rate, samples_mutated_rate, regulons_random_rate, miss_mutation_rate, miss_regulon_rate):
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


def getPredictedMapping(p_values, regulon_data):
    predicted_mapping = {}
    gene_count = len(p_values)
    sample_size = len(regulon_data)
    for gene_index in range(gene_count):
        row = p_values[gene_index]
        min_val = min(row)

        if min_val > 0:
            min_index = np.where(row == min_val)

            regulon_activity = []
            for i in range(sample_size):
                regulon_activity.append(regulon_data[i][min_index])
            effect = stats.mode(regulon_activity).mode
            predicted_mapping[gene_index] = [min_index[0][0], effect[0][0]]

    return predicted_mapping

def get_predicted_mapping(p_values, gene_data, regulon_data, alpha):
    predicted_mapping = {}
    gene_count = len(p_values)
    sample_size = len(regulon_data)
    for gene_index in range(gene_count):
        row = p_values[gene_index]
        min_val = min(row)

        if min_val < alpha:
            min_index = np.where(row == min_val)

            regulon_activity = []
            for i in range(sample_size):
                if gene_data[i][gene_index] == 1:
                    regulon_activity.append(regulon_data[i][min_index])
            effect = stats.mode(regulon_activity).mode
            predicted_mapping[gene_index] = [min_index[0][0], effect[0][0]]

    return predicted_mapping

def compareSets(pvalues, mutation_associations, alpha):
    gene_set = set()
    regulon_set = set()

    gene_count = len(pvalues)
    regulon_count = len(pvalues[0])

    for gene_index in range(gene_count):
        for regulon_index in range(regulon_count):
            if pvalues[gene_index][regulon_index] < alpha:
                gene_set.add(gene_index)
                regulon_set.add(regulon_index)

    printt("size of predicted mutated gene set: {}".format(len(gene_set)))
    printt("size of predicted mutated regulon set: {}".format(len(regulon_set)))
    printt("size of actual mutated regulon set: {}".format(len(mutation_associations)))

    correct = 0
    for gene_index in gene_set:
        if gene_index in mutation_associations and mutation_associations[gene_index][0] in regulon_set:
            correct+=1

    printt("number correct: {}/{}".format(correct, len(mutation_associations)))
    return correct

def run_tests():
    sample_size = 3000
    gene_count = 10000
    regulon_count = 1000
    genes_mutated_count = 100
    samples_mutated_rate = [0.2] # percentage of samples with mutated genes
    genes_random_rate = [0.2] # probability not mutated gene is observed as mutated
    regulons_random_rate = [0.2, 0.5, 0.8] # random distribution of regulon activity among non-affected regulons
    miss_mutation_rate = [0.2, 0.5, 0.8] # probability of there being a mutation but missing it
    miss_regulon_rate = [0.2, 0.5, 0.8] # probability that activity of associated regulon is not expected

    for i in samples_mutated_rate:
        for j in genes_random_rate:
            for k in regulons_random_rate:
                for l in miss_mutation_rate:
                    for m in miss_regulon_rate:
                        causal_association(sample_size, gene_count, regulon_count, genes_mutated_count, i, j, k, l, m)
                        print('\n')

def causal_association(sample_size, gene_count, regulon_count, genes_mutated_count, samples_mutated_rate, genes_random_rate, regulons_random_rate, miss_mutation_rate, miss_regulon_rate):
    start_time_timer = time.time()
    start_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # generate cell data

    test = "mann whitney u test"

    output_folder = './output/v2/'

    gene_data, regulon_data, mutation_associations = generateData(sample_size=sample_size, gene_count=gene_count,
                                             regulon_count=regulon_count, genes_mutated_count=genes_mutated_count,
                                             genes_random_rate=genes_random_rate, samples_mutated_rate=samples_mutated_rate,
                                             regulons_random_rate=regulons_random_rate, miss_mutation_rate=miss_mutation_rate,
                                                                  miss_regulon_rate=miss_regulon_rate)

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
        except:
            dill.dump_session(output_folder + file_name + '.dill')

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
    # # perform t test
    # # t_test_predicted, t_test_p_values = t_test(cell_data, variant_count)
    # # t_test_incorrect = calculatePercentage(actual_mapping, t_test_predicted, "t-test")
    # # auc_score(actual_mapping, t_test_predicted, t_test_p_values, "t-test", regulon_count, variant_count)
    #
    # dill.dump_session('./debug.dill')










def chisquare(cell_data):
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
