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
from collections import Counter
import time
from annotate_on_target import annotate_with_on_target_scores
import json
from multiprocessing import Process, Array

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
default_sliding_window_length = 30


def generate_parser():
    """
    Generates the parser for appropriate parsing, then returns the parser
    :return: the parser
    """
    parser = argparse.ArgumentParser()

    # positional arguments
    parser.add_argument("reads_file", type=str, help="provide the filename of the reads file")
    parser.add_argument("target_file", type=str,
                        help="provide the filename of the file where there is the target region. Specify with full "
                             "path.")
    parser.add_argument("scores_file", type=str, help="provide the name of the file where you wish to write the scores")
    parser.add_argument("max_hd", type=int,
                        help="provide the maximum hamming distance between the candidate gRNA and the genomic k-mer "
                             "that you wish for the pipeline to analyze. Make sure the hamming distance is between 0 "
                             "and 3, inclusive")

    # optional arguments
    parser.add_argument("-J", "--jf_filename", type=str,
                        help="Jellyfish binary filename (to avoid generating the file repeatedly)")
    parser.add_argument("-H", "--jf_histo_filename", type=str,
                        help="Jellyfish histogram filename (computed by 'jellyfish histo' command)")
    parser.add_argument("-eout", "--em_ouput_filename", type=str,
                        help="EM-computed values to be written in this file")
    parser.add_argument("-ein", "--em_input_filename", type=str,
                        help="EM-computed values to be read from this file")
    parser.add_argument("-m", "--max_copy_number", type=int, help="enter the highest number of times you think a "
                                                                  "genome region may repeat. Default: 50.", default=max_limit_count)
    parser.add_argument("-w", "--target_sliding_window_size", type=int, help="the size of target sliding window")
    parser.add_argument("-f", "--savgol_filter_window", type=int,
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
    parser.add_argument("-j", "--jf_threads", type=int, help="Number of jellyfish threads you want (default: 32)", default=32)
    parser.add_argument("-b", "--bt2_threads", type=int, help="Number of bowtie2 threads you want (default: 32)", default=32)
    parser.add_argument("-S", "--samtools_threads", type=int, help="Number of samtools threads you want (default: 32)", default=32)
    parser.add_argument("-B", "--sort_threads", type=int, help="Number of threads you want to sort bam file (default: "
                                                               "16)", default=16)
    parser.add_argument("-p", "--pilon_threads", type=int, help="Number of pilon threads you want (default: 32)", default=32)

    parser.add_argument("-l", "--logfilename", type=str,
                        help="Name of the log file", default="krispmer_log.log")
    parser.add_argument("-k", "--krispmer_threads", type=int, help="Number of threads you want to process grnas (default: 32)",
                        default=32)
    return parser


def parse_arguments(arg_string = None):
    """
    parses the arguments
    :return: the parsed argument
    """
    print('Parsing arguments...')
    parser = generate_parser()
    args_after_parsing = parser.parse_args(arg_string)

    # performing sanity checks on max_hd
    if args_after_parsing.max_hd > 3 or args_after_parsing.max_hd < 0:
        logging.info('Out of range Hamming-distance value received!')
        raise ValueError('Out of range Hamming-distance value received!')

    # performing sanity checks on savgol-filter window
    if args_after_parsing.savgol_filter_window is not None:
        if (args_after_parsing.savgol_filter_window % 2) == 0 or args_after_parsing.savgol_filter_window < 0:
            logging.info('Savgol-filter window size must be odd positive integer!')
            raise ValueError('Savgol-filter window size must be odd positive integer!')

    print('Parsing successful!\n')
    return args_after_parsing


def initial_jellyfish_run(reads_file_for_jellyfish, jf_threads):
    """
    k-mer counting using jellyfish. takes the reads file, counts k-mers and stores in another file named
    'jf_mer_count_file'
    :param reads_file_for_jellyfish: the reads file, fasta or fastq :return: None
    :param jf_threads: number of jellyfish threads
    """
    if not os.path.isdir('./krispmer_temp'):
        print('Creating krispmer_temp directory as it does not exist.')
        cmd = 'mkdir krispmer_temp'
        cmd_args = cmd.split(' ')
        subprocess.call(cmd_args)
    jf_command = "jellyfish count -m " + str(
        candidate_length) + " -s 100M -o " + jf_count_file + " -t " + str(jf_threads) + " -C " + reads_file_for_jellyfish
    jf_command_args = jf_command.split(" ")
    subprocess.call(jf_command_args)
    return jf_count_file


def generate_k_spectrum_histogram(jellyfish_filename, histo_output_file=hist_output):
    """
    generate the histogram using jellyfish command, then store the data in a dictionary and return that
    :param jellyfish_file: the jf binary file path
    :param histo_output_file: the file where you want to write the histo data
    :return: a dictionary containing the histogram information
    """
    histo_command = 'jellyfish histo ' + jellyfish_filename
    histo_command_args = histo_command.split(' ')
    res = subprocess.check_output(histo_command_args)
    with open(histo_output_file, 'w') as f:
        f.write(res)
    f.close()
    return read_histogram(histo_output_file)


def generate_k_spectrum_of_target_and_count(target_string, jellyfish_count_file, max_k_limit):
    """
    k-spectrum of target, then count the k-mers found within the target, then generate the histogram
    :type max_k_limit: int
    :param target_string: the target string
    :param jellyfish_count_file: jellyfish binary file (jellyfish.QueryMerFile)
    :param max_k_limit: max value upto which the histogram is to be generated
    :return: the histogram data in a dictionary as k_spectrum, and the counts of k-mers indexed as positions
    """
    # a pair is returned
    # pair.first = the k-spectrum histogram of k-mers taken only from the target region
    # pair.second = a hash-map that has keys:positions in target, values:count of a k-mer in that position
    k = candidate_length
    target = target_string
    length = len(target)
    a = set()
    counts_in_positions = {}
    k_spectrum = {}
    #qf = jellyfish.QueryMerFile(jellyfish_count_file)
    qf = jellyfish_count_file
    for i in range(length - k + 1):
        subst = target[i:i + k]
        mer = jellyfish.MerDNA(subst)
        mer.canonicalize()
        count = qf[mer]
        counts_in_positions[i] = count
        if count == 0:
            logging.info("Count = 0 for substring " + subst)
            continue
        if subst not in a:
            a.add(subst)
            if count in k_spectrum.keys():
                k_spectrum[count] += 1
            else:
                k_spectrum[count] = 1
    return k_spectrum, counts_in_positions


def sort_second(val):
    return val[1]


def complement(seq):
    complement_char = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    bases = list(seq)
    bases = [complement_char[base] for base in bases]
    return ''.join(bases)


def reverse_complement(s):
    return complement(s[::-1])


logsum = {}


def get_poisson(count, k):
    logarithm = -read_coverage * count + k * np.log(read_coverage * count)
    if k not in logsum.keys():
        logsum[k] = np.sum(np.log(np.arange(2, k + 1, 1.0)))
    logarithm -= logsum[k]
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


"""
return list should have the same number of elements as candidates_count_dictionary keys
"""
def annotate_guides_with_score_parallel(candidates_count_dictionary, jellyfish_filename, priors, posteriors,
                                        max_hd,
                                        target_string,
                                        return_list):
    index = 0
    list_candidates = []
    for candidate in list(candidates_count_dictionary.keys()):
        strand_type = candidates_count_dictionary[candidate][0]
        trie = generate_adjacent_mers(candidate, max_hd)
        value1 = value2 = 0.0
        flag = True
        for mer in trie.keys():
            if strand_type == '+':
                cp = get_score(candidate, mer)
            else:
                cp = get_score(reverse_complement(candidate), reverse_complement(mer))
            qf = jellyfish_filename
            merDNA = jellyfish.MerDNA(mer)
            merDNA.canonicalize()
            k = qf[merDNA]
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
        score = 1.0 * value2 / value1
        return_list[index] = score
        index = index + 1


def annotate_guides_with_score(candidates_count_dictionary, window_copy_numbers, jellyfish_filename, priors,
                               posteriors, max_hd, target_string, target_coverage, window_size, target_length, num_threads):
    candidates = list(candidates_count_dictionary.keys())
    candidates_per_thread = len(candidates) / num_threads
    return_lists = []
    process_list = []
    for i in range(num_threads):
        low_index = i * candidates_per_thread
        high_index = min(low_index + candidates_per_thread, len(candidates))
        candidates_count_dictionary_thread = {k:candidates_count_dictionary[k] for k in candidates[low_index:high_index]}
        return_list_this_thread = Array('d', len(list(candidates_count_dictionary_thread.keys())) * [-1.0])
        return_lists.append(return_list_this_thread)

        p = Process(target=annotate_guides_with_score_parallel, args=(candidates_count_dictionary_thread,
                                                                      jellyfish_filename, priors, posteriors, max_hd,
                                                                      target_string, return_list_this_thread))
        process_list.append(p)
        p.start()

    for process in process_list:
        process.join()

    all_scores = []
    for return_list in return_lists:
        all_scores = all_scores + list(return_list)

    annotated_candidates = []
    for (candidate, score) in list(zip(candidates, all_scores)):
        annotated_candidates.append([candidate, score, -1, None, candidates_count_dictionary[candidate][0], -1, -1])

    logging.info('DONE processing all candidates! Sorting...')
    annotated_candidates.sort(key=sort_second)

    logging.info('Final list:')
    for annotated_candidate in annotated_candidates:
        logging.info(annotated_candidate)
    return annotated_candidates


def krispmer_main(parsed_args):
    """
    takes the arguments parsed by the arg parser, then performs the pipeline of the tool
    :param parsed_args: the arguments passed to the programs
    :return: a list of all guides annotated with off-target activity
    """
    # get all arguments
    # already have in this version

    # initial jellyfish run
    if parsed_args.jf_filename is None:
        start_time = time.time()
        print('Jellyfish binary file not provided. Doing an initial run of Jellyfish to count'
              ' k-mers in reads (may take some time).')
        logging.info('No jf file given. Doing an initial run of Jellyfish to count k-mers in reads.')
        jellyfish_binary_filename = initial_jellyfish_run(parsed_args.reads_file, parsed_args.jf_threads)
        jellyfish_query_mer_file = jellyfish.QueryMerFile(jellyfish_binary_filename)
        logging.info('Completed the initial run.\n')
        print('Completed the run successfully.\n')
        end_time = time.time()
        logging.info('Time needed: ' + str(end_time-start_time) + ' seconds\n')
    else:
        print('Jellyfish binary file given as input. No need to run Jellyfish. Accessing the file.')
        logging.info('Jellyfish binary file given as input. No need to run Jellyfish. Accessing the file.\n')
        jellyfish_binary_filename = parsed_args.jf_filename
        if not os.path.isfile(jellyfish_binary_filename):
            print('The Jellyfish binary file does not exist. Exiting...')
            logging.info('The Jellyfish binary file does not exist. Exiting...')
            exit(-1)
        print('File exists.')
        jellyfish_query_mer_file = jellyfish.QueryMerFile(jellyfish_binary_filename)

    # generate k-mer spectrum histogram data
    if parsed_args.jf_histo_filename is None:
        start_time = time.time()
        print('Generating histogram data from Jellyfish counted database (may take some time).\n')
        logging.info('Generating histogram data from the initial Jellyfish database.')
        histogram_data_dictionary = generate_k_spectrum_histogram(jellyfish_binary_filename)
        logging.info('Finished generating histogram data\n')
        end_time = time.time()
        logging.info('Time needed: ' + str(end_time - start_time) + ' seconds\n')
    else:
        histo_filename = parsed_args.jf_histo_filename
        print('Jellyfish histogram file ' + histo_filename + ' is given input. Accessing file.')
        if not os.path.isfile(histo_filename):
            print('File does not exist. Exiting')
            logging.info('Histogram file does not exist. Exiting.\n')
            exit(-1)
        print('File exists. Reading the histogram.')
        logging.info('Histogram file exists. Reading the histogram.\n')
        histogram_data_dictionary = read_histogram(histo_filename)

    # determine priors, posteriors and read-coverage using EM
    # if these already available in a file, then do not do EM
    global read_coverage
    global max_limit_count
    global savgol_filter_window
    global max_k
    if parsed_args.max_copy_number is not None:
        max_limit_count = parsed_args.max_copy_number
    if parsed_args.savgol_filter_window is not None:
        savgol_filter_window = parsed_args.savgol_filter_window

    if parsed_args.em_input_filename is not None:
        print('EM-computed values in file: ' + parsed_args.em_input_filename + ', no need to run EM.')
        logging.info('EM-computed values in file: ' + parsed_args.em_input_filename)
        logging.info('Skipping EM.')

        json_file = open(parsed_args.em_input_filename, 'r')
        v = json.loads(json.load(json_file))
        priors, posteriors, read_coverage, max_k = v[0], v[1], v[2], v[3]
        json_file.close()

    else:
        print('Calculating priors using EM (may take some time).')
        logging.info('Calculating priors...')
        start_time = time.time()
        priors, posteriors, read_coverage = determine_priors_posteriors(histogram_data_dictionary,
                                                                    max_limit_count, savgol_filter_window)
        max_k = len(posteriors)
        logging.info('Finished calculating priors and posteriors\n')
        end_time = time.time()
        logging.info('Time needed: ' + str(end_time - start_time) + ' seconds\n')

    if parsed_args.em_ouput_filename is not None:
        json_file = open(parsed_args.em_ouput_filename, 'w')
        json.dump(json.dumps([priors.tolist(), posteriors, read_coverage, max_k]), json_file)
        json_file.close()

    # initialize poisson probability table
    global probability_table
    probability_table = [[-1] * max_k for i in range(max_limit_count)]  # type: List[List[int]]

    # next processing are target-specific processing. until now: general
    # personalized gRNA as a string
    start_time = time.time()
    if parsed_args.detect_variant:
        logging.info('Generating personalized version of the target\n')
        print('Generating personalized version of the target\n')
        modified_target_string = detect_variant(parsed_args.target_file, parsed_args.reads_file,
                                                parsed_args.bt2_threads,
                                                parsed_args.samtools_threads,
                                                parsed_args.sort_threads,
                                                parsed_args.pilon_threads)
        logging.info('Personalized target identified')
        logging.info('The target is: ' + modified_target_string)
        print('Personalized target identified\n')
    else:
        modified_target_string = read_target_region(parsed_args.target_file)
        logging.info('Proceeding with given target: ' + modified_target_string)
        print('Proceeding with given target: ' + modified_target_string + '\n')
    end_time = time.time()
    logging.info('Time needed: ' + str(end_time - start_time) + ' seconds\n')

    # determine all candidate list
    start_time = time.time()
    print('Generating all potential candidates.\n')
    logging.info('Generating list of potential candidates...')
    candidates_count_dictionary = candidate_generator.get_list_of_candidates(modified_target_string, pam, grna_length,
                                                                             parsed_args.stop,
                                                                             parsed_args.consider_negative,
                                                                             parsed_args.altPAMs)
    logging.info('Finished generating list of potential candidates...\n')
    end_time = time.time()
    logging.info('Time needed: ' + str(end_time - start_time) + ' seconds\n')

    # perform MLE to determine target coverage
    print ('Determining copy-number of target in genome.\n')
    logging.info('Starting MLE to determine copy-number of target in genome.')
    start_time = time.time()
    global target_coverage
    k_spectrum_data_in_target, position_count_dict = generate_k_spectrum_of_target_and_count(modified_target_string,
                                                                                             jellyfish_query_mer_file,
                                                                                             max_k)
    logging.info('The k-spectrum restricted with-in the k-mers of target string:')
    logging.info(k_spectrum_data_in_target)
    target_coverage = get_target_coverage(k_spectrum_data_in_target, read_coverage, parsed_args.max_copy_number)
    logging.info('The target is estimated_to_appear ' + str(target_coverage) + ' times')
    if int(target_coverage) > 1:
        logging.info('Warning: the target string may be a segmental-duplication region')
        print('Warning: the target string may be a segmental-duplication region')

    window_length = parsed_args.target_sliding_window_size
    if window_length is None:
        window_length = default_sliding_window_length
    window_copy_numbers = {}
    for i in range(len(modified_target_string) - window_length + 1):
        restricted_counts = [position_count_dict[key] for key in range(i,i+window_length-candidate_length+1)]
        restricted_k_spectrum_in_this_window = {k:restricted_counts.count(k) for k in restricted_counts}
        estimated_copy_number_of_window = get_target_coverage(restricted_k_spectrum_in_this_window,
                                                              read_coverage,
                                                              parsed_args.max_copy_number)
        window_copy_numbers[i] = estimated_copy_number_of_window

    logging.info(window_copy_numbers)
    end_time = time.time()
    logging.info('Time needed: ' + str(end_time - start_time) + ' seconds\n')
    # if variance too much, then throw another warning
    variance = np.var(window_copy_numbers.values())
    if variance > 1.0:
        logging.info('Warning: copy-numbers of target k-mers have high variation. The target may not be unique')
        print('Warning: copy-numbers of target k-mers have high variation. The target may not be unique')

    # annotate all guides
    print('Processing total ' + str(len(list(candidates_count_dictionary.keys()))) + ' candidate gRNAs\n')
    logging.info('Processing total ' + str(len(list(candidates_count_dictionary.keys()))) + ' candidate gRNAs')
    start_time = time.time()
    # todo: no need of window_copy_numbers now
    list_candidates = annotate_guides_with_score(candidates_count_dictionary,
                                                 window_copy_numbers,
                                                 jellyfish_query_mer_file,
                                                 priors,
                                                 posteriors,
                                                 parsed_args.max_hd,
                                                 modified_target_string,
                                                 target_coverage,
                                                 window_length,
                                                 len(modified_target_string),
                                                 parsed_args.krispmer_threads)

    # filter guides using cut-off score
    if parsed_args.cutoff_score is not None:
        list_candidates = [c for c in list_candidates if c[1] <= parsed_args.cutoff_score]

    print('Finished processing.')
    end_time = time.time()
    logging.info('Time needed: ' + str(end_time - start_time) + ' seconds\n')

    on_target_annotated_guides = annotate_with_on_target_scores(list_candidates, modified_target_string)
    for ot_annotated_guide in on_target_annotated_guides:
        logging.info(ot_annotated_guide)

    return list_candidates


def main_func():
    args = parse_arguments()
    kr_start_time = time.time()
    print(args.logfilename)
    logging.basicConfig(filename=args.logfilename, filemode='w', level=logging.INFO, format='%(message)s')
    logging.info(args.__dict__)
    logging.info('krispmer run initiated at' + str(kr_start_time))

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

    logging.info('Exiting\n')
    kr_end_time = time.time()
    logging.info('Total time needed to run: ' + str(kr_end_time-kr_start_time) + 's\n')

# krispmer -n combined.fastq staphylo-target.fasta out 0 -J krispmer_temp/jf_binary_file.jf -H krispmer_temp/k_spectrum_histo_data
if __name__ == '__main__':
    print('Hello there! Welcome to krispmerv0!')
    main_func()
