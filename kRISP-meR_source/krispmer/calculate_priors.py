__author__ = 'Mahmudur Rahman Hera'

import numpy as np
import pandas as pd
from scipy.signal import savgol_filter
import logging
import time

threshold = 0.1
logsum_k = {}


class poisson:
    def __init__(self, mean):
        self.mean = mean

    def get_probability(self, k):
        total = -self.mean + k * np.log(self.mean)
        if k not in logsum_k.keys():
            logsum_k[k] = np.sum(np.log(np.arange(2, k + 1, 1.0)))
        total -= logsum_k[k]
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
    values = list(histo_data.values())
    updated_values = savgol_filter(values, savgol_filter_window, 3)
    logging.info("Updated values after smoothing:")
    for v in updated_values:
        logging.info(v)
    logging.info("")
    lower = higher = -1
    zipped_data = list(zip(histo_data.keys(), updated_values))
    for i in range(1, len(zipped_data)):
        if zipped_data[i - 1][1] > zipped_data[i][1] and zipped_data[i + 1][1] > zipped_data[i][1]:
            lower = zipped_data[i][0]
        if zipped_data[i - 1][1] < zipped_data[i][1] and zipped_data[i + 1][1] < zipped_data[i][1]:
            higher = zipped_data[i][0]
            break
    return lower, higher


def expectation_maximization(data_points_raw, init_mean, max_count, max_k):
    logging.info('The histogram data:')
    logging.info(data_points_raw)
    # filter data points so that too large points are not there

    max_key = max_k + 1

    keys = data_points_raw.keys()
    data_points = np.zeros((max_key))
    multiplied_data_points = np.zeros((max_key))
    for k in range(1, max_key):
        if k in keys:
            data_points[k] = data_points_raw[k]
            multiplied_data_points[k] = data_points_raw[k] * k

    logging.info('After filtering, the histogram data:')
    logging.info(data_points)

    # initial poisson distributions, the first is intended to capture the error component
    poissons = [poisson(1.01)]
    for x in range(1, max_count):
        poissons.append(poisson(float(x * init_mean)))

    # probability table to contain probabilities
    probabilities = np.zeros((max_count, max_key))
    # responsibility table to contain degree of membership of a certain data point to the poissons
    tendencies = np.zeros((max_count, max_key))
    # priors: the mixing probabilities
    priors = np.array([1.0 / max_count] * max_count)

    # perform EM
    print('Starting EM.')
    logging.info('starting EM algorithm.')
    iteration_count = 1
    previous_lambda = init_mean
    while True:
        logging.info('Running iteration ' + str(iteration_count))
        print('Running iteration ' + str(iteration_count))

        logging.info('Calculating new probabilities...')
        for j in range(max_key):
            for i in range(max_count):
                probabilities[i][j] = poissons[i].get_probability(j)

        logging.debug('The probabilities are:')
        logging.debug(probabilities)

        logging.info('Calculating degree of membership to belong to the poissons...')

        denominators = {}
        for j in range(max_key):
            denominators[j] = np.dot(probabilities[:, j], priors)

        for i in range(max_count):
            for j in range(max_key):
                tendencies[i][j] = 1.0 * probabilities[i][j] * priors[i] / denominators[j]

        logging.debug('The degrees of membership are:')
        logging.debug(tendencies)

        logging.info('Calculating the new means...')
        for i in range(max_count):
            total1 = np.dot(tendencies[i, :], multiplied_data_points)
            total2 = np.dot(tendencies[i, :], data_points)
            poissons[i].update_mean(1.0 * total1 / total2)

        logging.info('The new means are:')
        means = [p.get_mean() for p in poissons]
        logging.info(means)

        iteration_count = iteration_count + 1
        logging.info('Estimating the average coverage from the means...')

        lambdas = []
        ownership_factors = []
        for i in range(1, max_count):
            lambdas.append(1.0 * poissons[i].get_mean() / float(i))
            ownership_factors.append(np.dot(tendencies[i, :], data_points))

        estimated_averaged_lambda = 1.0 * np.dot(np.array(lambdas), np.array(ownership_factors)) / sum(
            ownership_factors)
        logging.info('The estimated coverage for this iteration is: ' + str(estimated_averaged_lambda))

        logging.info('Updating mean of the poissons according to the estimated average..')
        for i in range(1, max_count):
            poissons[i].update_mean(i * estimated_averaged_lambda)

        logging.info('Recalculating priors...')
        for i in range(0, max_count):
            val1 = np.dot(tendencies[i, :], data_points)
            val2 = np.sum(data_points)
            priors[i] = 1.0 * val1 / val2

        logging.info('Updated priors are: ' + str(priors))

        logging.info('Finished iteration ' + str(iteration_count))
        logging.info('Coverage estimation at this iteration: ' + str(estimated_averaged_lambda))
        if abs(previous_lambda - estimated_averaged_lambda) < threshold:
            break
        else:
            previous_lambda = estimated_averaged_lambda

    print('EM converged.\n')
    return estimated_averaged_lambda, priors, poissons, probabilities


def determine_priors_posteriors(histo_data, max_priors, savgol_filter_window):
    """
    determines the priors, posteriors and the read-coverage
    :param histo_data: dictionary containing histogram data
    :param max_priors: number of poissons
    :param savgol_filter_window: size of the savgol filter, must be an odd positive integer
    :return: priors in an array, posteriors in an array, read-coverage, the inversion point
    """
    inv_point, init_coverage = determine_points(histo_data, savgol_filter_window)
    logging.info('Initial coverage is: ' + str(init_coverage))
    max_k = int((max_priors * init_coverage))
    estimated_kmer_coverage, priors, poissons, probabilities = expectation_maximization(histo_data, init_coverage,
                                                                                         max_priors, max_k)

    logging.info('EM converged with coverage = ' + str(estimated_kmer_coverage))
    logging.debug('Priors are: ' + str(priors))

    # determine posteriors here
    print('Calculating posterior probabilities.\n')
    logging.info('Calculating the posterior probabilities...')

    posteriors_ = [0.0] * (max_k + 1)
    for k in range(1, max_k + 1):
        posteriors_[k] = np.dot(probabilities[:, k], priors)

    logging.debug('The posterior probabilities are:')
    logging.debug(posteriors_)

    return priors, posteriors_, estimated_kmer_coverage


# this main method only exists for testing purposes
if __name__ == '__main__':
    start = time.time()
    d = read_histogram('histogram-human')
    priors, posteriors, estimated_kmer_coverage = determine_priors_posteriors(d, 100, 5)
    print(estimated_kmer_coverage)
    print(sum(posteriors))
    end = time.time()
    print(end - start)
    # print (determine_points(d))
