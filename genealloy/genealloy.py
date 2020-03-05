try:
    from Levenshtein import distance
except ImportError:
    has_pkgs = False
else:
    has_pkgs = True

from .codontable import (codon_to_aa, aa_to_codon_extended, codon_extended_to_aa)

def convert_seq_to_codons(seq):
    seq_codons = [seq[i : i + 3] for i in range(0, len(seq), 3)]

    return seq_codons


def convert_codonlist_to_tuplelist(seq_codons, codon_to_codon_extended):
    codon_extended = [None] * len(seq_codons)
    for i, codon in enumerate(seq_codons):
        codon_extended[i] = codon_to_codon_extended[codon]

    return codon_extended


def find_matches(parasite_tuplelist, host_tuplelist):
    index_of_matches = dict()
    len_parasite = len(parasite_tuplelist)
    len_host = len(host_tuplelist)

    for begin_index in range(0, len_host - len_parasite + 1):
        matching_sequence = []
        parasite_index = 0
        has_match = True
        for parasite_index in range(0, len_parasite):
            common_codons = set(parasite_tuplelist[parasite_index]) & set(host_tuplelist[begin_index + parasite_index])
            if  common_codons == set():
                has_match = False
                break
            else:
                matching_sequence.append(tuple(common_codons))
                parasite_index += 1
        if has_match:
            index_of_matches[begin_index] = matching_sequence

    return index_of_matches


def report_results(index_of_matches):
    if index_of_matches == dict():
        print("The guest (parasite) sequence cannot be inserted into the host sequence.")
        return
    
    print("The number of matching substrings in host is: %d," % len(index_of_matches))
    positions = list(index_of_matches.keys())
    print("at host codon position(s): %s" % positions)


def calculate_aa_solutions(seq_codons, index_of_matches):
    parasite_aa_seq = ''.join(codon_to_aa[c] for c in seq_codons)
    solutions = dict()
    solutions['parasite'] = parasite_aa_seq
    for key, match in index_of_matches.items():
        solutions[key] = ''
        for index, codon in enumerate(seq_codons):
            aa = codon_to_aa[codon]
            if set(aa_to_codon_extended[aa]) & set(match[index]) != set():
                solutions[key] += aa
            else:
                solutions[key] += codon_extended_to_aa[match[index][0]]  # take the 1st aa

    return solutions


def calculate_levenshtein_distances(aa_solutions):
    if not has_pkgs:
        raise ImportError("Function requires the python-Levenshtein package.")

    parasite = aa_solutions['parasite']
    del aa_solutions['parasite']

    levenshtein_distances = dict()
    for key, aa_seq in aa_solutions.items():
        levenshtein_distances[key] = distance(parasite, aa_seq)

    return levenshtein_distances
