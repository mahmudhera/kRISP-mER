def complement(seq):
    complement_char = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    bases = list(seq)
    bases = [complement_char[base] for base in bases]
    return ''.join(bases)


def reverse_complement(s):
    return complement(s[::-1])


def annotate_with_on_target_scores(guideRNAs, target_string):
    """Append the 30 bp Rule Set 2 context to each guide row.

    Rule Set 2 wants 4 bp upstream + 20 bp guide + 3 bp PAM + 3 bp
    downstream, read 5'->3' on the strand the guide cuts.

    generate_all_candidates stores a guide as target[i-gRNA_len : i+3],
    so the key ALREADY INCLUDES THE PAM and is 23 bp, not 20. For
    minus-strand guides the key is reverse-complemented back to forward
    orientation, so every guide is found in target_string.

    Plus strand:  [4][----- key 23 -----][3]
    Minus strand: the key read forward is rc(guide+PAM), so in forward
    coordinates the layout is [3][----- key 23 -----][4], and the window
    must be reverse complemented after slicing.

    Before 6 Nov 2026 this searched target_string but sliced the reverse
    complement for minus-strand guides, mixing two coordinate systems.
    """
    L = len(target_string)
    for annotated_guideRNA in guideRNAs:
        key = annotated_guideRNA[0]
        strand = annotated_guideRNA[4]

        pos = target_string.find(key)
        if pos < 0:
            continue

        if strand == '+':
            start, end = pos - 4, pos + len(key) + 3
        else:
            start, end = pos - 3, pos + len(key) + 4

        if start < 0 or end > L:
            continue

        window = target_string[start:end]
        if strand == '-':
            window = reverse_complement(window)
        annotated_guideRNA.append(window)
    return guideRNAs
