import pytest

import genealloy as ga
from genealloy.genealloy import get_next_letter
from genealloy.genealloy import get_letter_in_next_triplet


def test_convert_seq_to_codons():
    seq_codons = ga.convert_seq_to_codons("GTACCCGCTGCG")
    assert seq_codons == ["GTA", "CCC", "GCT", "GCG"]


def test_convert_codonlist_to_tuplelist():
    codon_to_codon_extended = ga.generate_swaptable(
        ga.codon_to_aa, ga.aa_to_codon_extended
    )

    result = ga.convert_codonlist_to_tuplelist(
        ["GTA", "CCC", "GCT", "GCG"], codon_to_codon_extended
    )
    assert result == [("GTX",), ("CCX",), ("GCX",), ("GCX",)]


def test_get_next_letter():
    seq_tuplelist = [("GTX",), ("CCX",), ("GCX",), ("GCX",)]
    assert get_next_letter(seq_tuplelist, [0, 0, 2]) == ("C", [1, 0, 0])

    with pytest.raises(Exception):
        assert get_next_letter(seq_tuplelist, [3, 0, 2])

    assert get_next_letter(seq_tuplelist, [1, 0, 1]) == ("X", [1, 0, 2])


def test_get_letter_in_next_triplet():
    seq_tuplelist = [("GTX", "CCX"), ("GCX",)]
    assert get_letter_in_next_triplet(seq_tuplelist, [0, 0, 2]) == ("C", [0, 1, 0])
    with pytest.raises(Exception):
        assert get_letter_in_next_triplet(seq_tuplelist, [1, 0, 1])


def test_compare_then_get_letter_recursively():
    host = "AACCGCGGTGAGCATTTCAAGGGTACAACGGTTCAA"
    parasite = "CGCGGTGAGCATTTCAAGGGTACAACGGTTCAA"
    host_triplets = ga.convert_seq_to_codons(host)
    parasite_triplets = ga.convert_seq_to_codons(parasite)
    codon_to_codon_extended = ga.generate_swaptable(
        ga.codon_to_aa, ga.aa_to_codon_extended
    )

    host_tuplelist = ga.convert_codonlist_to_tuplelist(
        host_triplets, codon_to_codon_extended
    )
    parasite_tuplelist = ga.convert_codonlist_to_tuplelist(
        parasite_triplets, codon_to_codon_extended
    )

    host_letter = ("C", [1, 0, 0])
    parasite_letter = ("C", [0, 0, 0])

    result = ga.compare_then_get_letter_recursively(
        host_tuplelist, parasite_tuplelist, host_letter, parasite_letter
    )

    expected = (
        "Finished parasite sequence, match found! Ending host position:[11, 0, 2]"
    )

    assert result == expected

    no_result = ga.compare_then_get_letter_recursively(
        host_tuplelist, (["CAA"], ["CCC"]), host_letter, parasite_letter
    )

    assert no_result == "No match for this starting codon position"
