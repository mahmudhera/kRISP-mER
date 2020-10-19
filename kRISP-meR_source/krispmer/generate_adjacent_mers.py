__author__ = 'Mahmudur Rahman Hera'

import trie
from itertools import chain, combinations, product


def hamming_circle(s, n, alphabet, trie):
    """Generate strings over alphabet whose Hamming distance from s is
    exactly n.
    sorted(hamming_circle('abc', 0, 'abc'))
    ['abc']
    sorted(hamming_circle('abc', 1, 'abc'))
    ['aac', 'aba', 'abb', 'acc', 'bbc', 'cbc']
    sorted(hamming_circle('aaa', 2, 'ab'))
    ['abb', 'bab', 'bba']
    """
    for positions in combinations(range(len(s)), n):
        for replacements in product(range(len(alphabet) - 1), repeat=n):
            cousin = list(s)
            for p, r in zip(positions, replacements):
                if cousin[p] == alphabet[r]:
                    cousin[p] = alphabet[-1]
                else:
                    cousin[p] = alphabet[r]
            trie[''.join(cousin)] = 1


def hamming_ball(s, n, alphabet, trie):
    """Generate strings over alphabet whose Hamming distance from s is
	less than or equal to n.
    sorted(hamming_ball('abc', 0, 'abc'))
	['abc']
	sorted(hamming_ball('abc', 1, 'abc'))
	['aac', 'aba', 'abb', 'abc', 'acc', 'bbc', 'cbc']
	sorted(hamming_ball('aaa', 2, 'ab'))
	['aaa', 'aab', 'aba', 'abb', 'baa', 'bab', 'bba']
    """
    for i in range(n + 1):
        hamming_circle(s, i, alphabet, trie)


def generate_adjacent_mers(sequence, max_hamming_distance):
    alphabet = 'AGCT'
    t = trie.trie()
    hamming_ball(sequence, max_hamming_distance, alphabet, t)
    return t
