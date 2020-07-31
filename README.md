![image](https://img.shields.io/badge/%20-linux-orange)
![image](https://img.shields.io/badge/%20-python-blue)
![image](https://img.shields.io/badge/crispr-referencefree-yellowgreen)
# kRISP-mER
Reference free guide RNA designing tool for CRISPR

## What kRISP-mER is
This is a tool to generate personalized guide RNAs for CRISPR without using a reference genome. 

Most genome-wide guideRNA designer tools have to use the whole reference genome to populate their database. This limits their usage for organisms with incomplete reference genome. Instead of using a reference, kRISP-mER works using the sequenced reads, and a genomic target location (location where CRISPR cleavage is intended). Using the sequenced reads only, kRISP-mER is able to design variant-aware guideRNAs and predict those with minimized off-target activity.

This tool is designed for:
1. Linux
1. Python 2.7

## Dependencies you need to install
The following installation instructions are _only to help you out_. These installation instructions are **NOT** mandatory to follow. You can install these any way you like. However, if you are having a hard time doing so by yourself, you may find these instructions useful.
* **samtools 1.0 or higher**: You can install samtools using the following commands:
```shell script
sudo apt-get update -y
sudo apt-get install -y samtools
```
* **Biopython**: install using: `pip install biopython`
* **Python binding of Jellyfish** (gmarics project). Need SWIG 3.0 or higher. This is quite tricky. Need to assess this in detail later. Tried to do the following:
```shell script
./configure --prefix=$HOME --enable-python-binding --enable-swig
make -j 4
make install
```
Python binding means that if you try to import jellyfish from a python script (the code is: `import jellyfish`), that will work and you will be able to invoke Jellyfish program from python.
If the installation steps mentioned above does not bind with python, (although the documentation of Jellyfish does say that this should): you may try the swig binding instructions from https://github.com/gmarcais/Jellyfish/blob/master/swig/Readme.md).
```shell script
cd swig/python
python setup.py build
python setup.py install --prefix=$PREFIX
```
* **Bowtie2**: You can install this with With `Bioconda`. With `Bioconda` installed, you should be able to install Bowtie 2 with `conda install bowtie2`. Details are found here: `http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#obtaining-bowtie-2`
* **Java 1.8 or higher**: kRISP-mER uses Pilon to find a better sequence after seq reads are aligned to a base sequence. Pilon is run as a plain jar file. In order to check this requirement, you may want to test `java -version`
* **numpy**: Simple installation with pip: `pip install numpy`. You may also work with some package manager, such as anaconda, in which case numpy should come built in.
* **scipy**: Simple installation with pip: `pip install scipy`. You may also work with some package manager, such as anaconda, in which case scipy should come built in.
* **sklearn**: If the package manager you are using does not already have scikit-learn installed, you can install using `pip install scikit-learn==0.16.1` (This very specific version is important to determine on target activity scores)
* **pickle**: This python package is required to determine on-target-activity as well. This should already be installed in python 2 and 3. If not, you need to manually install this.

## Installation
Once you have the dependencies installed, installing kRISP-mER is easy. You need to:
1. **Download** the github repository
1. **Go to** the directory `kRISP-mER source`
1. **Install** by running the command: `python setup.py install`

## Inputs
kRISP-mER takes as input the following:
1. The sequenced reads as a FASTQ file
2. The target region as a FASTA file

Note that the file containing the target region should have a line `>chromosome_name` as its header.

## Running
To run kRISP-mER, enter the following command in shell after successful installation.

```shell script
krispmer READS_FILENAME TARGET_FILENAME OUTPUT_FILENAME MAX_HD
```

Here, **MAX_HD** is the maximum Hammind distance that kRISP-mER will consider when scanning target sites for a particular guideRNA. The output of the program is saved in the **OUTPUT_FILENAME** file in csv format.

## Interpreting the output
After the program exits, you will see four columns after opening the output file. The first is the guideRNA in the + strand, the second is the guideRNA in the - strand (reverse complement of the first). The third column stores the estimated inverse-specificity score assigned to a particular guideRNA. And finally, the fourth column stores the strand in which the **NGG** pam was found by kRISP-mER. 

Besides the scores output file, you will also see a directory named `krispmer_temp` and a file named `krispmer.log`. The `krispmer_temp` directory contains the temporary files created when executing the program. You can use `-r` flag to delete them automatically (see detailed usage below). `krispmer.log` file contains detailed steps of the whole pipeline.

## Detailed usage
Usage:
```shell script
krispmer [-h] [-m MAX_COPY_NUMBER] [-w SAVGOL_FILTER_WINDOW] [-s] [-v]
                [-n] [-c CUTOFF_SCORE] [-a ALTPAMS [ALTPAMS ...]] [-r]
                <reads_file> <target_file> <scores_file> <max_hd>
```
### Positional (mandatory) arguments
1. reads_file
1. target_file
1. scores_file
1. max_hd

kRISP-mER allows you to design guide RNAs with WGS shotgun reads (in a FASTQ file), and a target-region (a FASTA file). Besides these two, you also have to tell the program the number of mismatches to consider when scanning for target sites against a particular guideRNA. kRISP-mER allows upto 3 mismatches. kRISP-mER does not consider indels (like other established gRNA designing tools). You also have to tell the program the name of the output csv file, where the guideRNAs along with their inverse-specificity scores and strand information is to be stored.

### Non-positional (optional) arguments
1. `-h`: You can see help with `-h` flag
1. `-m`: You can set the maximum number of times a region may repeat in the genome using `-m` flag. Default: 50. 
1. `-w`: kRISP-mER uses savgol smoothing filter to smoothen the k-spectrum data before applying Expectation-Maximization to estimate read coverage. You can set the window size of the filter using `-w` flag
1. `-s`: You can choose to exclude the guideRNAs that contain a stop-codon using the flag `-s`
1. `-v`: You can choose to polish the target region and personalize that for the individual whose sequenced reads are being used. You can do so using the flag `-v`. A long pipeline using bowtie2, samtools and pilon will start executing.
1. `-n`: You can choose to scan the PAM sequences in the -ve strand using the flag `-n`
1. `-c`: You can set a cut-off score of the inverse-specificity using the flag `-c`. The guideRNAs with score higher than that will be dropped.
1. `-a PAM1 PAM2 ...`: You can provide kRISP-mER with a list of PAMs to consider with `-a` flag. By default, NGG PAMs are considered. 
1. `-r`: You can choose to remove the temporary files automatically using the flag `-r`
1. `-j`: You can set the number of threads you want to use to count the k-mers in the sequenced reads using Jellyfish using the flag `-j`
1. `-b`: You can set the number of bowtie2 threads with the flag `-b`
1. `-S`: You can set the number of samtools threads with the flag `-S`
1. `-B`: You can set the number of threads you want to use to sort the intermediate BAM file using the flag `-B`
1. `-p`: You can set the number of threads to be used in Pilon using the flag `-p`

### Example
```shell script
krispmer -nvr read.fastq target.fasta out 1
```
This will use the sequenced reads in the file `read.fastq` and the target sequence found in `target.fasta` and write the scores as csv in the file `out`.

This specific command will:
* Determine the target sequence for the particular individual whose reads are being used
* Scan the -ve strand (reverse complement of the target) for PAM sequences to identify guideRNAs
* Consider a Hamming-distance of 1 when scanning for target-sites for a guideRNA
* Remove all temporary files

## Issues and support
If you find any issues, please feel free to:
* either create an issue in github
* or email me at `moc.liamg@39arehdumham` with the logfile (the file named `krispmer.log`)

I will try to get back to you as soon as possible.
