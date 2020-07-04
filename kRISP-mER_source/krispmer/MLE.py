import pandas as pd
import numpy as np

def generate_dic_from_k_spectrum(filename):
    """
    Given the filename, generates the dictionary of key values.
    Sample input spectrum file:
        1 2
        2 6
        3 2
    This means there were two occurences where copy number = 1 was observed,
    six occurences where copy number = 2 was observed and so on.
    
    returns: the dictionary
    """
    df = pd.read_csv(filename, delimiter=' ', header=None, names = ['key', 'data'])
    return pd.Series(df.data.values,index=df.key).to_dict()

def get_target_coverage(k_spectrum_data, read_coverage):
    """
    The # of times the target appears in the genome is estimated using the spectrum data.
    The reads in genome are modeled as a mixture of poissons, with mean = coverage, 2*coverage, 3*coverage and so on.
    When the tgt appears i times in the genome, the probability of observing the k-spectrum from the poisson with mean i*coverage will be the maximum.
    
    The log likelihoods are calculated and the corresponding max count is then returned.
    """
    sum_all_values = sum(k_spectrum_data.values())
    dot = np.dot(k_spectrum_data.keys(), k_spectrum_data.values())
    max_llhood = -read_coverage*sum_all_values + np.log(read_coverage)*dot
    tgt_coverage = 1
    for i in range(2,10):
        llhood = -read_coverage*sum_all_values*i + np.log(read_coverage*i)*dot
        if llhood > max_llhood:
            max_llhood = llhood
            tgt_coverage = i
    return tgt_coverage

def get_target_coverage_after_refining(k_spectrum_data, read_coverage):
    """
    The inversion point serves the purpose of pruning the observation of the error reads.
    :param k_spectrum_data: dictionary containing k-spectrum count histogram from target string's k-mers
    :param read_coverage: estimated read coverage
    :return: # of times target appears in the genome
    """
    # seems like refining is not important. So, don't refine
    #refined_spectrum = {k:k_spectrum_data[k] for k in k_spectrum_data.keys() if k >= inversion_point}
    refined_spectrum = k_spectrum_data
    return get_target_coverage(refined_spectrum, read_coverage)

if __name__=='__main__':
    dic = generate_dic_from_k_spectrum('k-spectrum3.txt')
    tgt_cvg = get_target_coverage_after_refining(dic, 5, 2)
    print(tgt_cvg)
