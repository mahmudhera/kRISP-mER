from krispmer import reverse_complement

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
