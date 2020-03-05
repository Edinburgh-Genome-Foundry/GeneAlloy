import os
from Bio import SeqIO
import genealloy as ga

data_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "data")
records = list(SeqIO.parse(os.path.join(data_dir, "example.fa"), "fasta"))

E1 = records[0]  # Host sequence
E2 = records[1]
E3 = records[2]  # E3_Same_strand_same_frame_same_aa
E4 = records[3]  # E4_Same_strand_nonmatching
E5 = records[4]  # E5_Same_strand_multiple_matches

host = str(E1.seq)
codon_to_codon_extended = ga.generate_swaptable(ga.codon_to_aa, ga.aa_to_codon_extended)

parasite = str(E2.seq)
host_codons = ga.convert_seq_to_codons(host)
host_tuplelist = ga.convert_codonlist_to_tuplelist(host_codons, codon_to_codon_extended)

parasite_codons = ga.convert_seq_to_codons(parasite)
parasite_tuplelist = ga.convert_codonlist_to_tuplelist(
    parasite_codons, codon_to_codon_extended
)

index_of_matches = ga.find_matches(parasite_tuplelist, host_tuplelist)
aa_solutions = ga.calculate_aa_solutions(parasite_codons, index_of_matches)


def test_find_matches():
    index_of_matches = ga.find_matches(parasite_tuplelist, host_tuplelist)
    assert len(index_of_matches[55]) == 11
    assert index_of_matches[55][1] == ("GGX",)
    assert len(index_of_matches) == 2


def test_calculate_aa_solutions():
    aa_solutions = ga.calculate_aa_solutions(parasite_codons, index_of_matches)
    assert len(aa_solutions) == 3
    assert aa_solutions[55] == "RGEHFKGTTVQ"


def test_levenshtein_distances():
    levenshtein_distances = ga.calculate_levenshtein_distances(aa_solutions)
    assert levenshtein_distances[55] == 0
    assert levenshtein_distances[155] == 0
