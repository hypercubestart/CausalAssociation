import datetime
from scipy import stats
import math
from single_cell import SingleCell
import collections
import dill
import numpy

def printt(message):
    """Print message with timestamp
        :param message: string
    """
    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t {}".format(message)))
    return None

def t_test(cell_data, variant_count, incorrect_mapping = None):
    """Causal association test using t test
    :param cell_data: [SingleCell]
    :param variant_count: integer number of variants
    :param incorrect_mapping: if not None, prints out incorrect associations
    """
    #TODO: SHOULD USE MANN WHITNEY U TEST OVER T-TEST
    threshold = 0.01
    ANOVA_threshold = 0.01
    predicted_mapping = {}
    aggregate_p_values = []
    for mutationIndex in range(variant_count):

        if (incorrect_mapping is not None and mutationIndex in incorrect_mapping):
            printt('{} is wrong'.format(mutationIndex))

        # split regulons by with and without mutation
        with_mutation = []
        without_mutation = []
        for cell in cell_data:
            if cell.get_variant(mutationIndex) == 1:
                with_mutation.append(cell.copy_regulon())
            else:
                without_mutation.append(cell.copy_regulon())

        # TODO: filter to find significant difference between two groups

        # check that each list has at least one regulon
        if len(with_mutation) == 0 or len(without_mutation) == 0:
            printt("T-test: Mutation at index {} cannot be associated because len == 0".format(mutationIndex))
            continue

        # evaluate p values
        p_values = []

        for i in range(len(with_mutation[0])):
            with_mutation_regulon_activity = []
            without_mutation_regulon_activity = []

            if len(without_mutation) <= 1 or len(without_mutation) <= 1:
                # cannot determine associations because not enough data
                p_values.append(math.nan)
                continue

            for j in range(len(with_mutation)):
                with_mutation_regulon_activity.append(with_mutation[j][i])
            for j in range(len(without_mutation)):
                without_mutation_regulon_activity.append(without_mutation[j][i])

            # conduct t test between every pair of values between with_mutation_average and without_mutation_average
            t_test = stats.ttest_ind(with_mutation_regulon_activity, without_mutation_regulon_activity, equal_var=True)

            ## p-value can be nan if length <= 1 or no variance (all values are the same)
            if math.isnan(t_test.pvalue):
                # cannot determine associations
                p_values.append(1)
            else:
                p_values.append(t_test.pvalue)

        aggregate_p_values.append(p_values)

        # ignore if no statistically significant difference in mean is found
        if (min(p_values) > threshold):
            continue

        # get indices with lowest p value
        index = p_values.index(min(p_values))
        # use most common regulon activity to determine effect of variant on regulon
        counts = collections.Counter([regulon[index] for regulon in with_mutation]).most_common()[0][0]
        predicted_mapping[mutationIndex] = [index, counts]
    return predicted_mapping, aggregate_p_values

if __name__ == "__main__":
    dill.load_session('./debug.dill')
    t_test(cell_data, variant_count, t_test_incorrect)