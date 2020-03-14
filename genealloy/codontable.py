"""Codon to amino acid dictionary"""
codon_to_aa = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


"""Amino acid to codon dictionary, using extended nucleotide letters"""
aa_to_codon_extended = {
    "A": ["GCX"],
    "B": ["RAY"],
    "C": ["TGY"],
    "D": ["GAY"],
    "E": ["GAR"],
    "F": ["TTY"],
    "G": ["GGX"],
    "H": ["CAY"],
    "I": ["ATH"],
    #     'J': ['-'],
    "K": ["AAR"],
    "L": ["TTR", "CTX", "YTR"],
    "M": ["ATG"],
    "N": ["AAY"],
    #     'O': ['-'],
    "P": ["CCX"],
    "Q": ["CAR"],
    "R": ["CGX", "AGR", "MGR"],
    "S": ["TCX", "AGY"],
    "T": ["ACX"],
    # "U": ["-"],
    "V": ["GTX"],
    "W": ["TGG"],
    "X": ["XXX"],
    "Y": ["TAY"],
    "Z": ["SAR"],
    ".": ["-"],  # no amino acid (deletion or gap)
    "*": ["TAR", "TRA"],  # STOP codons
}


"""Codon to amino acid dictionary, using extended nucleotide letters"""
codon_extended_to_aa = {
    "GCX": "A",
    "RAY": "B",
    "TGY": "C",
    "GAY": "D",
    "GAR": "E",
    "TTY": "F",
    "GGX": "G",
    "CAY": "H",
    "ATH": "I",
    # "-": "J",
    "AAR": "K",
    "TTR": "L",
    "CTX": "L",
    "YTR": "L",
    "ATG": "M",
    "AAY": "N",
    # "-": "O",
    "CCX": "P",
    "CAR": "Q",
    "CGX": "R",
    "AGR": "R",
    "MGR": "R",
    "TCX": "S",
    "AGY": "S",
    "ACX": "T",
    # "-": "U",
    "GTX": "V",
    "TGG": "W",
    "XXX": "X",
    "TAY": "Y",
    "SAR": "Z",
    # "-": ".",  # no amino acid (deletion or gap)
    "TAR": "*",  # STOP codon
    "TRA": "*",  # STOP codon
}


"""Extended nucleotide letter to nucleotide letter dictionary"""
ambiguity_code_to_nt_set = {
    "A": {"A"},
    "G": {"G"},
    "C": {"C"},
    "T": {"T"},
    "Y": {"C", "T"},
    "R": {"A", "G"},
    "W": {"A", "T"},
    "S": {"G", "C"},
    "K": {"T", "G"},
    "M": {"C", "A"},
    "D": {"A", "G", "T"},
    "V": {"A", "C", "G"},
    "H": {"A", "C", "T"},
    "B": {"C", "G", "T"},
    "X": {"A", "C", "G", "T"},
    "N": {"A", "C", "G", "T"},
}


"""Extended nucleotide letter to complement letter dictionary"""
complement_table = {
    "A": "T",
    "G": "C",
    "C": "G",
    "T": "A",
    "Y": "R",
    "R": "Y",
    "W": "W",
    "S": "S",
    "K": "M",
    "M": "K",
    "D": "H",
    "V": "B",
    "H": "D",
    "B": "V",
    "X": "X",
    "N": "N",
}

allowed_aa_transitions = {
    "A": ["G", "A", "V", "L", "I"],
    "B": ["D", "E", "N", "Q", "B", "Z"],
    "C": ["S", "C", "M", "T"],
    "D": ["D", "E", "N", "Q", "B", "Z"],
    "E": ["D", "E", "N", "Q", "B", "Z"],
    "F": ["F", "Y", "W"],
    "G": ["G", "A", "V", "L", "I"],
    "H": ["H", "K", "R"],
    "I": ["G", "A", "V", "L", "I"],
    #     'J': ['J'],
    "K": ["H", "K", "R"],
    "L": ["G", "A", "V", "L", "I"],
    "M": ["S", "C", "M", "T"],
    "N": ["D", "E", "N", "Q", "B", "Z"],
    #     'O': ['O'],
    "P": ["P"],
    "Q": ["D", "E", "N", "Q", "B", "Z"],
    "R": ["H", "K", "R"],
    "S": ["S", "C", "M", "T"],
    "T": ["S", "C", "M", "T"],
    # "U": ['S', 'C', 'U', 'M', 'T'],
    "V": ["G", "A", "V", "L", "I"],
    "W": ["F", "Y", "W"],
    "X": ["X"],
    "Y": ["F", "Y", "W"],
    "Z": ["D", "E", "N", "Q", "B", "Z"],
    ".": ["."],
    "*": ["*"],
}


def make_transition_dictionary(aa_to_codon_extended, allowed_aa_transitions):
    transition_dictionary = {}
    for aa, aa_list in allowed_aa_transitions.items():
        codons = []
        for aa_element in aa_list:
            triplets_to_add = aa_to_codon_extended[aa_element]
            codons = codons + triplets_to_add
        transition_dictionary[aa] = codons
    return transition_dictionary


# TODO rename aa_to_codon_extended --> transition_dictionary and references elsewhere
def generate_swaptable(codon_to_aa, aa_to_codon_extended):
    """Generate a codon to extended codon dictionary"""
    codon_to_codon_extended = dict()
    for k, v in codon_to_aa.items():
        codon_to_codon_extended[k] = tuple(aa_to_codon_extended[v])

    return codon_to_codon_extended


def compare_letters(letter1, letter2, table=ambiguity_code_to_nt_set):
    """Compare two extended nucleotide letters and return True if they match"""
    set1 = table[letter1]
    set2 = table[letter2]

    if set1 & set2 != set():
        is_match = True
    else:
        is_match = False

    return is_match
