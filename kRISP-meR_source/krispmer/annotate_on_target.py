def complement(seq):
    complement_char = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    bases = list(seq)
    bases = [complement_char[base] for base in bases]
    return ''.join(bases)


def reverse_complement(s):
    return complement(s[::-1])

# returns with just the context of 30-mer
def annotate_with_on_target_scores(guideRNAs, target_string):
    reverse_target_string = reverse_complement(target_string)
    for annotated_guideRNA in guideRNAs:
        gRNA_string = annotated_guideRNA[0]
        strand = annotated_guideRNA[4]
        if strand == '+':
            search_string = target_string
        else:
            search_string = reverse_target_string
        pos = target_string.find(gRNA_string)
        if pos < 4:
            continue
        annotated_guideRNA.append(search_string[pos-4 : pos+26])
    return guideRNAs
