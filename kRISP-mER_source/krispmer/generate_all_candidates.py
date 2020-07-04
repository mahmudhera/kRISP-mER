__author__ = 'Mahmudur Rahman Hera'

from Bio import trie


class NFiller:
    def __init__(self, PAM):
        self.PAM = PAM
        self.bases = ['A', 'C', 'G', 'T']

    @staticmethod
    def findall(p, s):
        """Yields all the positions of the pattern p in string s."""
        i = s.find(p)
        while i != -1:
            yield i
            i = s.find(p, i + 1)

    def recursive_N_filler(self, whole, lst, ipositions):
        if (not ipositions):
            return lst
        for c in self.bases:
            lst2 = list(lst)
            lst2[ipositions[-1]] = c
            ip = ipositions[:-1]
            whole.append(self.recursive_N_filler(whole, lst2, ip))

    def get_list(self):
        PAM = self.PAM
        ipositions = [index for index in self.findall("N", PAM)]
        if (not ipositions):
            return [PAM]
        PAMlist = list(PAM)
        whole = []
        self.recursive_N_filler(whole, PAMlist, ipositions)
        w = [x for x in whole if x is not None]
        # print w
        # print(len(w))
        s = [''.join(x) for x in w]
        return s


def findall(p, s):
    i = s.find(p)
    while i != -1:
        yield i
        i = s.find(p, i + 1)


def find_candidates(target, PAM, gRNA_len):
    """
	#~ Say, 
	#~ target: ATATATATGCATATAGCTATAGCATGCAT
	#~ pAM: 	TGC
	#~ gRNA_len: 5
	"""
    ww = [(i - gRNA_len, target[i - gRNA_len:i + len(PAM)]) for i in findall(PAM, target)]
    return [x for x in ww if x[0] >= 0]


def complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    bases = list(seq)
    bases = [complement[base] for base in bases]
    return ''.join(bases)


def reverse_complement(s):
    return complement(s[::-1])


# [(i-gRNA_len, target[i-gRNA_len:i+len(PAM)]) for i in findall(PAM, target)]
def get_list_of_candidates(target_string, PAM, gRNA_length, exclude_stop_codons, consider_negative, alt_pams):
    target = target_string
    if consider_negative:
        target_rev = reverse_complement(target)
    if alt_pams is None:
        PAMs = NFiller(PAM).get_list()
    else:
        for pam in alt_pams:
            if len(pam) != 3:
                raise ValueError('Length of one or more PAMs not set to 3')
            for character in pam:
                if character != 'A' and character != 'C' and character != 'G' and character != 'T':
                    raise ValueError('Invalid PAM has been entered')
        PAMs = alt_pams
    candidates_rev = []
    candidates = []
    for PAM in PAMs:
        candidates.extend(find_candidates(target, PAM, gRNA_length))
        if consider_negative:
            candidates_rev.extend(find_candidates(target_rev, PAM, gRNA_length))
    trie_dic = trie.trie()
    for candidate in candidates:
        key = candidate[1]
        if exclude_stop_codons and ('TAG' in key or 'TAA' in key or 'TGA' in key):
            continue
        if key not in trie_dic.keys():
            trie_dic[key] = '+'
    for candidate in candidates_rev:
        if exclude_stop_codons and ('TAG' in candidate[1] or 'TAA' in candidate[1] or 'TGA' in candidate[1]):
            continue
        key = reverse_complement(candidate[1])
        if key not in trie_dic.keys():
            trie_dic[key] = '-'
    return trie_dic