"""
Rule Set 2 Cas9 on-target cutting efficiency score for target site from Doench et al. Nature Biotechnology 2016
http://www.nature.com/nbt/journal/v34/n2/abs/nbt.3437.html
"""

import pickle
from pkg_resources import resource_filename

mms = resource_filename(__name__,'CFD_scoring/mismatch_score.pkl')
pams = resource_filename(__name__,'CFD_scoring/pam_scores.pkl')

# this computes the reverse complement of an RNA string
def revcom(s):
    basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A'}
    letters = list(s[::-1])
    letters = [basecomp[base] for base in letters]
    return ''.join(letters)

# this calculates the CFD scores
def calc_cfd(wt, sg, pam, mm_scores, pam_scores):
    score = 1.0
    sg = sg.replace('T', 'U')
    wt = wt.replace('T', 'U')
    s_list = list(sg)
    wt_list = list(wt)
    for i, sl in enumerate(s_list):
        if wt_list[i] == sl:
            score *= 1
        else:
            try:
                key = 'r' + wt_list[i] + ':d' + revcom(sl) + ',' + str(i + 1)
                score *= mm_scores[key]
            except KeyError:
                continue
    score *= pam_scores[pam]
    return score


# get mis-match PAM scores
def get_mm_pam_scores(mms, pams):
    try:
        mm_scores = pickle.load(open(mms, 'rb'))
        pam_scores = pickle.load(open(pams, 'rb'))
        return (mm_scores, pam_scores)
    except:
        raise Exception("Could not find file with mismatch scores or PAM scores")


mm_scores, pam_scores = get_mm_pam_scores(mms, pams)


# this is the interface called from outside
def get_score(candidate, sequence):
    pam = candidate[-2:]
    sg = candidate[:-3]
    cfd_score = calc_cfd(sequence, sg, pam, mm_scores, pam_scores)
    return cfd_score
