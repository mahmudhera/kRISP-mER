__author__ = 'Mahmudur Rahman Hera'

import numpy as np
import argparse
from generate_personalized_target import read_target_region, detect_variant
import generate_all_candidates as candidate_generator
import pandas as pd
from calculate_priors import determine_priors_posteriors, read_histogram
import dna_jellyfish as jellyfish
from MLE import get_target_coverage
from generate_adjacent_mers import generate_adjacent_mers
from get_cfd_score import get_score
import logging
import os
import subprocess

pam = "NGG"
grna_length = 20
candidate_length = 23
jf_count_file = "krispmer_temp/jf_binary_file.jf"
probability_table = []
max_k = -1
max_limit_count = 50
savgol_filter_window = 9
hist_output = 'krispmer_temp/k_spectrum_histo_data'
read_coverage = -1
target_coverage = -1


def generate_parser():
    """
    Generates the parser for appropriate parsing, then returns the parser
    :return: the parser
    """
    parser = argparse.ArgumentParser()

    # positional arguments
    parser.add_argument("reads_file", type=str, help="provide the filename of the reads file with path")
    parser.add_argument("target_file", type=str,
                        help="provide the filename of the file where there is the target region. Specify with full "
                             "path.")
    parser.add_argument("scores_file", type=str, help="provide the name of the file where you wish to write the scores")
    parser.add_argument("max_hd", type=int,
                        help="provide the maximum hamming distance between the candidate gRNA and the genomic k-mer "
                             "that you wish for the pipeline to analyze. Make sure the hamming distance is between 0 "
                             "and 3, inclusive")

    # optional arguments
    parser.add_argument("-m", "--max_copy_number", type=int, help="enter the highest number of times you think a "
                                                                  "genome region may repeat. Default: 50.")
    parser.add_argument("-w", "--savgol_filter_window", type=int,
                        help="enter the window size of savgol filter. This must be an odd integer. Default: 9.")
    parser.add_argument("-s", "--stop", help="identify stop codons in guides and exclude those guides",
                        action="store_true")
    parser.add_argument("-v", "--detect_variant", help="uses bowtie2, samtools and pilon to detect genomic variation",
                        action="store_true")
    parser.add_argument("-n", "--consider_negative", help="scans the -ve strand for finding the guides",
                        action="store_true")
    parser.add_argument("-c", "--cutoff_score", type=float, help="enter the maximum score of the guides you wish to be")
    parser.add_argument("-a", "--altPAMs", type=str, help="Type in the NGG pams that you want to work with", nargs='+')
    parser.add_argument("-r", "--remove_temp", help="Remove temporarily created files",
                        action="store_true")
    return parser


def parse_arguments(arg_string = None):
    """
    parses the arguments
    :return: the parsed argument
    """
    print('Parsing arguments...')
    logging.info('Parsing the passed arguments...')
    parser = generate_parser()
    args_after_parsing = parser.parse_args(arg_string)

    if args_after_parsing.max_hd > 3 or args_after_parsing.max_hd < 0:
        logging.info('Out of range Hamming-distance value received!')
        raise ValueError('Out of range Hamming-distance value received!')

    if args_after_parsing.savgol_filter_window is not None:
        if (args_after_parsing.savgol_filter_window % 2) == 0 or args_after_parsing.savgol_filter_window < 0:
            logging.info('Savgol-filter window size must be odd positive integer!')
            raise ValueError('Savgol-filter window size must be odd positive integer!')

    logging.info('Finished parsing the arguments.\n')
    print('Parsing successful!\n')
    return args_after_parsing


def initial_jellyfish_run(reads_file_for_jellyfish):
    """
    k-mer counting using jellyfish. takes the reads file, counts k-mers and stores in another file named
    'jf_mer_count_file'
    :param reads_file_for_jellyfish: the reads file, fasta or fastq :return: None
    """
    if not os.path.isdir('./krispmer_temp'):
        print('Creating krispmer_temp directory as it does not exist.')
        cmd = 'mkdir krispmer_temp'
        cmd_args = cmd.split(' ')
        subprocess.call(cmd_args)
    jf_command = "jellyfish count -m " + str(
        candidate_length) + " -s 100M -o " + jf_count_file + " -t 20 -C " + reads_file_for_jellyfish
    jf_command_args = jf_command.split(" ")
    # todo: uncomment this later
    subprocess.call(jf_command_args)
    return jf_count_file


def generate_k_spectrum_histogram(jellyfish_file, histo_output_file=hist_output):
    """
    generate the histogram using jellyfish command, then store the data in a dictionary and return that
    :param jellyfish_file: the jf binary file path
    :param histo_output_file: the file where you want to write the histo data
    :return: a dictionary containing the histogram information
    """
    histo_command = 'jellyfish histo ' + jellyfish_file
    histo_command_args = histo_command.split(' ')
    res = subprocess.check_output(histo_command_args)
    with open(histo_output_file, 'w') as f:
        f.write(res)
    f.close()
    return read_histogram(histo_output_file)


def generate_k_spectrum_of_target_and_count(target_string, jellyfish_count_file, max_k_limit=200):
    """
    k-spectrum of target, then count the k-mers found within the target, then generate the histogram
    :type max_k_limit: int
    :param target_string: the target string
    :param jellyfish_count_file: jellyfish binary file name
    :param max_k_limit: max value upto which the histogram is to be generated
    :return: the histogram data in a dictionary
    """
    k = candidate_length
    target = target_string
    length = len(target)
    a = set()
    for i in range(length - k):
        a.add(target[i:i + k])
    lst = []
    qf = jellyfish.QueryMerFile(jellyfish_count_file)
    for subst in a:
        mer = jellyfish.MerDNA(subst)
        count = qf[mer]
        lst.append(count)
    dic = {}
    for i in range(max_k_limit):
        dic[i + 1] = lst.count(i + 1)
    return dic


def sort_second(val):
    return val[1]


def complement(seq):
    complement_char = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    bases = list(seq)
    bases = [complement_char[base] for base in bases]
    return ''.join(bases)


def reverse_complement(s):
    return complement(s[::-1])


def get_poisson(count, k):
    logarithm = -read_coverage * count + k * np.log(read_coverage * count)
    for i in range(2, k + 1):
        logarithm -= np.log(i)
    return np.exp(logarithm)


def get_probability(count, k):
    if k > max_k:
        return 0.0
    if probability_table[count][k] != -1:
        return probability_table[count][k]
    # total = 0.0
    for k1 in range(max_k):
        probability_table[count][k1] = get_poisson(count, k1)
        # total = total + probability_table[count][k1]
    # for k1 in range(max_k):
    #	probability_table[count][k1] = probability_table[count][k1]/total
    return probability_table[count][k]


def annotate_guides_with_score(candidates_count_dictionary, jellyfish_filename, priors, posteriors, max_hd,
                               target_string, target_coverage):
    iteration_count = 0
    list_candidates = []
    for candidate in list(candidates_count_dictionary.keys()):
        strand_type = candidates_count_dictionary[candidate]
        trie = generate_adjacent_mers(candidate, max_hd)
        value1 = value2 = 0.0
        logging.info('Processing candidate ' + candidate + '...')
        flag = True
        for mer in trie.keys():
            if strand_type == '+':
                cp = get_score(candidate, mer)
            else:
                cp = get_score(reverse_complement(candidate), reverse_complement(mer))
            qf = jellyfish.QueryMerFile(jellyfish_filename)
            merDNA = jellyfish.MerDNA(mer)
            rev_comp_merDNA = jellyfish.MerDNA(reverse_complement(mer))
            k = max(qf[merDNA], qf[rev_comp_merDNA])
            if k <= 0:
                continue
            if k >= max_k:
                flag = False
                break
            p = float(target_string.count(mer))
            accum = 0.0
            for count in range(1, max_limit_count):
                probability = get_probability(count, k)
                p_count = priors[count]
                p_k = posteriors[k]
                new_val = 1.0 * probability * count * p_count / p_k
                accum = accum + new_val
            value1 = value1 + cp * p
            value2 = value2 + cp * accum
        if value1 <= 0.0 or flag is False:
            continue
        score = 1.0 * value2 / (value1 * target_coverage)
        qf = jellyfish.QueryMerFile(jellyfish_filename)
        merDNA = jellyfish.MerDNA(candidate)
        k = max(qf[merDNA], qf[jellyfish.MerDNA(reverse_complement(candidate))])
        list_candidates.append((candidate, score, k, trie, strand_type))
        iteration_count = iteration_count + 1
        logging.info('Processed ' + str(iteration_count) + 'th gRNA: ' + candidate + ' with score= ' + str(score))
    logging.info('DONE processing all candidates! Sorting...')
    list_candidates.sort(key=sort_second)
    logging.info('Final list:')
    for annotated_candidate in list_candidates:
        logging.info(annotated_candidate)
    return list_candidates


def krispmer_main(parsed_args):
    """
    takes the arguments parsed by the arg parser, then performs the pipeline of the tool
    :param parsed_args: the arguments passed to the programs
    :return: a list of all guides annotated with off-target activity
    """
    # get all arguments
    # already have in this version

    # initial jellyfish run
    print('Doing an initial run of Jellyfish to count k-mers in reads (may take some time).')
    logging.info('Doing an initial run of Jellyfish to count k-mers in reads.')
    jellyfish_binary_file = initial_jellyfish_run(parsed_args.reads_file)
    logging.info('Completed the initial run.\n')
    print('Completed the run successfully.\n')

    # personalized gRNA as a string
    if parsed_args.detect_variant:
        logging.info('Generating personalized version of the target\n')
        print('Generating personalized version of the target\n')
        modified_target_string = detect_variant(parsed_args.target_file, parsed_args.reads_file)
        logging.info('Personalized target identified')
        logging.info('The target is: ' + modified_target_string)
        print('Personalized target identified\n')
    else:
        modified_target_string = read_target_region(parsed_args.target_file)
        logging.info('Proceeding with given target: ' + modified_target_string)
        print('Proceeding with given target: ' + modified_target_string + '\n')

    # generate k-mer spectrum histogram data
    print('Generating histogram data from Jellyfish counted database (may take some time).\n')
    logging.info('Generating histogram data from the initial Jellyfish database.')
    histogram_data_dictionary = generate_k_spectrum_histogram(jellyfish_binary_file)
    logging.info('Finished generating histogram data\n')

    # determine all candidate list
    print('Generating all potential candidates.\n')
    logging.info('Generating list of potential candidates...')
    candidates_count_dictionary = candidate_generator.get_list_of_candidates(modified_target_string, pam, grna_length,
                                                                             parsed_args.stop,
                                                                             parsed_args.consider_negative,
                                                                             parsed_args.altPAMs)
    logging.info('Finished generating list of potential candidates...\n')
    logging.info(candidates_count_dictionary)

    # determine priors, posteriors and read-coverage using EM
    global read_coverage
    global max_limit_count
    global savgol_filter_window
    if parsed_args.max_copy_number is not None:
        max_limit_count = parsed_args.max_copy_number
    if parsed_args.savgol_filter_window is not None:
        savgol_filter_window = parsed_args.savgol_filter_window
    print('Calculating priors using EM (may take some time).')
    logging.info('Calculating priors...')
    priors, posteriors, read_coverage = determine_priors_posteriors(histogram_data_dictionary,
                                                                    max_limit_count, savgol_filter_window)
    global max_k
    max_k = len(posteriors)
    logging.info('Finished calculating priors and posteriors\n')

    # initialize poisson probability table
    global probability_table
    probability_table = [[-1] * max_k for i in range(max_limit_count)]  # type: List[List[int]]

    # perform MLE to determine target coverage
    print ('Determining copy-number of target in genome.\n')
    logging.info('Starting MLE to determine copy-number of target in genome.')
    global target_coverage
    k_spectrum_data_in_target = generate_k_spectrum_of_target_and_count(modified_target_string, jellyfish_binary_file,
                                                                        max_k)
    logging.info('The k-spectrum restricted with-in the k-mers of target string:')
    logging.info(k_spectrum_data_in_target)
    target_coverage = get_target_coverage(k_spectrum_data_in_target, read_coverage)
    logging.info('The target appears ' + str(target_coverage) + ' times')

    # annotate all guides
    print('Processing total ' + str(len(list(candidates_count_dictionary.keys()))) + ' candidate gRNAs')
    logging.info('Processing total ' + str(len(list(candidates_count_dictionary.keys()))) + ' candidate gRNAs')
    list_candidates = annotate_guides_with_score(candidates_count_dictionary, jellyfish_binary_file, priors, posteriors,
                                                 parsed_args.max_hd, modified_target_string, target_coverage)

    # filter guides using cut-off score
    if parsed_args.cutoff_score is not None:
        list_candidates = [c for c in list_candidates if c[1] <= parsed_args.cutoff_score]

    print('Finished processing.')
    return list_candidates


def main_func():
    logging.basicConfig(filename='krispmer.log', filemode='w', level=logging.INFO, format='%(message)s')

    args = parse_arguments()
    logging.info(args.__dict__)

    gRNAs = krispmer_main(args)
    logging.info('Total number of grnas: ' + str(len(gRNAs)))

    print("Writing to output file.")
    logging.info('Writing to output file...')
    output_file = open(args.scores_file, 'w')
    output_file.write('tgt_in_plus,tgt_in_minus,inverse_specificity,strand\n')
    for gRNA in gRNAs:
        seq = gRNA[0]
        score = gRNA[1]
        strand = gRNA[4]
        output_file.write(seq + ',' + reverse_complement(seq) + ',' + str(score) + ',' + strand + '\n')
    output_file.close()

    if args.remove_temp:
        if os.path.isdir('./krispmer_temp'):
            logging.info('Cleaning up krispmer_temp directory.')
            print('Cleaning up krispmer_temp directory')
            cmd = 'rm -rf krispmer_temp'
            cmd_args = cmd.split(' ')
            subprocess.call(cmd_args)

    logging.info('Exiting')

# python krispmer.py -n ../run-data/read.fastq ../run-data/repeated_sequence.txt ../run-data/out 0
if __name__ == '__main__':
    main_func()