from setuptools import setup, find_packages, Extension

triemodule = Extension('trie', sources = ['trie/triemodule.c', 'trie/trie.c'])

setup(
    name = 'krispmer',
    version = '0.0.7',
    author= 'Mahmudur Rahman Hera',
    author_email= 'mahmudhera93@gmail.com',
    description = 'A tool to design CRISPR guideRNA without using reference genome',
    #long_description = long_description,
    #install_requires = [
    #    'biopython>=1.66','pysam==0.8.3','pyfaidx==0.4.7.1','bx-python==0.7.3'
    #],
    packages = find_packages(),
    ext_modules = [triemodule],
    package_data = {'krispmer' : ['pilon-1.23.jar',
                                  'CFD_scoring/mismatch_score.pkl',
                                  'CFD_scoring/pam_scores.pkl',
                                  'V3_model_full.pickle',
                                  'V3_model_nopos.pickle']},
    entry_points={
        'console_scripts': [
            'krispmer = krispmer.krispmer:main_func'
        ],
    }
)
