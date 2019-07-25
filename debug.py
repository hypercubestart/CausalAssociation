import dill
from single_cell import SingleCell
from scipy import stats
import math
import collections
import datetime

def printt(message):
    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t {}".format(message)))
    return None

def t_test_debug(cell_data, incorrect):
    threshold = 0.01
    ANOVA_threshold = 0.01
    predicted_mapping = {}
    aggregate_p_values = []
    for mutationIndex in range(len(cell_data[0].variants)):

        if (mutationIndex in incorrect):
            printt("mutation index {} is incorrect".format(mutationIndex))

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

            for j in range(len(with_mutation)):
                with_mutation_regulon_activity.append(with_mutation[j][i])
            for j in range(len(without_mutation)):
                without_mutation_regulon_activity.append(without_mutation[j][i])


            if len(with_mutation_regulon_activity) <= 1 or len(without_mutation_regulon_activity) <= 1:
                p_values.append(1)
                continue

            # conduct t test between every pair of values between with_mutation_average and without_mutation_average
            t_test = stats.ttest_ind(with_mutation_regulon_activity, without_mutation_regulon_activity, equal_var=False)

            ## p-value can be nan is length <= 1 or no variance
            if math.isnan(t_test.pvalue):
                p_values.append(1)
            else:
                p_values.append(t_test.pvalue)

        # ignore if no statistically significant difference in mean is found
        if (min(p_values) > threshold):
            continue

        # get indices with lowest p value
        index = p_values.index(min(p_values))
        # use most common regulon activity to determine effect of variant on regulon
        counts = collections.Counter([regulon[index] for regulon in with_mutation]).most_common()[0][0]
        predicted_mapping[mutationIndex] = [index, counts]
        aggregate_p_values.append(p_values)
    return predicted_mapping, aggregate_p_values

if __name__ == "__main__":

    dill.load_session('./debug.dill')
    t_test_debug(cell_data, t_test_incorrect)


