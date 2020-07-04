__author__ = 'Mahmudur Rahman Hera'

from typing import List, Any, Tuple, Iterator

import numpy as np
import math
import pandas as pd
from scipy.signal import savgol_filter
import logging

threshold = 0.01


class poisson:
    def __init__(self, mean):
        self.mean = mean

    def get_probability(self, k):
        total = -self.mean + k * np.log(self.mean)
        total -= sum([math.log(i) for i in range(2, k + 1)])
        return np.exp(total)

    def update_mean(self, mean):
        self.mean = mean

    def get_mean(self):
        return self.mean


# this only exists for testing purposes
def read_histogram(filename):
    histo_data = pd.read_csv(filename, delimiter=' ', names=['index', 'value'])
    dic = {}
    for t in zip(np.array(histo_data['index']), np.array(histo_data['value'])):
        dic[t[0]] = t[1]
    return dic


def determine_points(histo_data, savgol_filter_window):
    values = histo_data.values()
    updated_values = savgol_filter(values, savgol_filter_window, 2)
    lower = higher = -1
    zipped_data = zip(histo_data.keys(), updated_values)
    for i in range(1, len(zipped_data)):
        if zipped_data[i - 1][1] > zipped_data[i][1] and zipped_data[i + 1][1] > zipped_data[i][1]:
            lower = zipped_data[i][0]
        if zipped_data[i - 1][1] < zipped_data[i][1] and zipped_data[i + 1][1] < zipped_data[i][1]:
            higher = zipped_data[i][0]
            break
    return lower, higher


def expectation_maximization(data_points_raw, init_mean, max_count=50):
    logging.info('The histogram data:')
    logging.info(data_points_raw)
    # filter data points so that too large points are not there
    data_points = {k: v for k, v in data_points_raw.items() if k <= max_count * init_mean}

    logging.info('After filtering, the histogram data:')
    logging.info(data_points)

    # initial poisson distributions, the first is intended to capture the error component
    poissons = [poisson(1.01)]
    for x in range(1, max_count):
        poissons.append(poisson(float(x * init_mean)))

    # probability table to contain probabilities
    probabilities = [{k: 0.0 for k in data_points.keys()} for i in range(max_count)]
    # responsibility table to contain degree of membership of a certain data point to the poissons
    tendencies = [{k: 0.0 for k in data_points.keys()} for i in range(max_count)]
    # priors: the mixing probabilities
    priors = [1.0 / max_count] * max_count

    # perform EM
    print('Starting EM.')
    logging.info('starting EM algorithm.')
    iteration_count = 1
    previous_lambda = init_mean
    while True:
        logging.info('Running iteration ' + str(iteration_count))
        print('Running iteration ' + str(iteration_count))

        logging.info('Calculating new probabilities...')
        for j in data_points.keys():
            for i in range(max_count):
                probabilities[i][j] = poissons[i].get_probability(j)

        logging.debug('The probabilities are:')
        logging.debug(probabilities)

        logging.info('Calculating degree of membership to belong to the poissons...')
        for i in range(max_count):
            for j in data_points.keys():
                total = 0.0
                for ii in range(max_count):
                    total = total + probabilities[ii][j] * priors[ii]
                tendencies[i][j] = 1.0 * probabilities[i][j] * priors[i] / total

        logging.debug('The degrees of membership are:')
        logging.debug(tendencies)

        logging.info('Calculating the new means...')
        for i in range(max_count):
            total1 = 0.0
            total2 = 0.0
            for j in data_points.keys():
                total1 = total1 + tendencies[i][j] * j * data_points[j]
                total2 = total2 + tendencies[i][j] * data_points[j]
            poissons[i].update_mean(1.0 * total1 / total2)

        logging.info('The new means are:')
        means = [p.get_mean() for p in poissons]
        logging.info(means)

        iteration_count = iteration_count + 1
        logging.info('Estimating the average coverage from the means...')
        lambda_weighted_sum = 0.0
        ownership_factor_sum = 0.0
        ownership_factors = []
        for i in range(1, max_count):
            lambda_ = 1.0 * poissons[i].get_mean() / float(i)
            ownership_factor_this_poisson = 0.0
            for j in data_points.keys():
                ownership_factor_this_poisson += tendencies[i][j] * data_points[j]
            ownership_factors.append(ownership_factor_this_poisson)
            lambda_weighted_sum += lambda_ * ownership_factor_this_poisson
            ownership_factor_sum += ownership_factor_this_poisson
        estimated_averaged_lambda = lambda_weighted_sum / ownership_factor_sum
        logging.info('The estimated coverage for this iteration is: ' + str(estimated_averaged_lambda))

        logging.info('Updating mean of the poissons according to the estimated average..')
        for i in range(1, max_count):
            poissons[i].update_mean(i * estimated_averaged_lambda)

        logging.info('Recalculating priors...')
        for i in range(0, max_count):
            val1 = 0.0
            val2 = 0.0
            for j in data_points.keys():
                val1 += tendencies[i][j] * data_points[j]
                val2 += data_points[j]
            priors[i] = val1 / val2

        logging.info('Updated priors are: ' + str(priors))

        logging.info('Finished iteration ' + str(iteration_count))
        logging.info('Coverage estimation at this iteration: ' + str(estimated_averaged_lambda))
        if abs(previous_lambda - estimated_averaged_lambda) < threshold:
            break
        else:
            previous_lambda = estimated_averaged_lambda

    print ('EM converged.\n')
    return estimated_averaged_lambda, priors, poissons


def determine_priors_posteriors(histo_data, max_priors, savgol_filter_window):
    """
    determines the priors, posteriors and the read-coverage
    :param histo_data: dictionary containing histogram data
    :param max_priors: number of poissons
    :param savgol_filter_window: size of the savgol filter, must be an odd positive integer
    :return: priors in an array, posteriors in an array, read-coverage, the inversion point
    """
    inv_point, init_coverage = determine_points(histo_data, savgol_filter_window)
    estimated_kmer_coverage, priors, poissons = expectation_maximization(histo_data, init_coverage, max_priors)

    logging.info('EM converged with coverage = ' + str(estimated_kmer_coverage))
    logging.debug('Priors are: ' + str(priors))

    # determine posteriors here
    print('Calculating posterior probabilities.\n')
    logging.info('Calculating the posterior probabilities...')
    max_k = int((max_priors + 1) * estimated_kmer_coverage)
    posteriors_ = [0.0] * max_k
    for k in range(1, max_k):
        for i in range(len(poissons)):
            posteriors_[k] += poissons[i].get_probability(k) * priors[i]

    logging.debug('The posterior probabilities are:')
    logging.debug(posteriors_)

    return priors, posteriors_, estimated_kmer_coverage


# this main method only exists for testing purposes
if __name__ == '__main__':
    d = read_histogram('histo_real_data')
    priors, posteriors, estimated_kmer_coverage = determine_priors_posteriors(d, 30, 5)
    print(estimated_kmer_coverage)
    print(sum(posteriors))
    # print (determine_points(d))
