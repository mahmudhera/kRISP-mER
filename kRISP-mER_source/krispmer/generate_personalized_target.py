__author__ = 'Mahmudur Rahman Hera'

import subprocess
import logging
import os
from shutil import copyfile
from pkg_resources import resource_filename


def read_target_region(filename):
    """
    reads a fasta file and generates the target region as a string
    :param filename: fasta file (target region)
    :return: string without the line containing '>'
    """
    with open(filename) as f:
        content = f.readlines()
    content = [x.strip() for x in content][1:]
    string_ = ''.join(content)
    return string_.upper()


def detect_variant(target_filename, reads_filename):
    """
	Reads the target, aligns the seq reads on the target using bowtie2,
	uses Pilon to detect the genetic variations,
	then returns the improved target sequence.

    :param target_filename: the target region filename as a fasta
    :param reads_filename: the fastq/fasta reads file. Currently, only unpaired
    :return: returns the target region after detecting personalized variations as a string
	"""
    if not os.path.isdir('./krispmer_temp'):
        print('Creating krispmer_temp directory as it does not exist.')
        cmd = 'mkdir krispmer_temp'
        cmd_args = cmd.split(' ')
        subprocess.call(cmd_args)

    logging.info('Copying the pilon jar file into the temp directory.')
    print('Copying the pilon jar file into the temp directory.')
    pilon_source_filename = resource_filename(__name__, 'pilon-1.23.jar')
    pilon_destination_filename = 'krispmer_temp/pilon-1.23.jar'
    copyfile(pilon_source_filename, pilon_destination_filename)

    with open('krispmer_temp/commands_py.sh', 'w') as f:
        f.write('cd ./krispmer_temp\n')
        f.write('bowtie2-build -f ../' + target_filename + ' krispmer\n')
        f.write('bowtie2 --local -x krispmer -U ../' + reads_filename + ' -S sam_out.sam\n')
        f.write('samtools view -bS sam_out.sam > bam_out.bam\n')
        f.write('samtools sort bam_out.bam > bam_sorted.bam\n')
        f.write('samtools index bam_sorted.bam\n')
        f.write('java -Xmx16G -jar pilon-1.23.jar --genome ../' + target_filename + ' --unpaired bam_sorted.bam\n')
        f.close()
    return_code = subprocess.call(['bash', 'krispmer_temp/commands_py.sh'])
    if return_code != 0:
        logging.info('Something went wrong while computing personalized target!')
        logging.info('Exiting')
        exit(-1)
    variant_target = read_target_region('krispmer_temp/pilon.fasta')
    return variant_target


if __name__ == '__main__':
    variant = detect_variant('inputs/target20k.txt', 'inputs/read_staphylo.fastq')
    original_target = read_target_region('inputs/target20k.txt')
    print(len(variant))
    print(len(original_target))
