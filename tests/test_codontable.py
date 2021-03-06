import pytest

import genealloy as ga


def test_make_transition_dictionary():
    transition_dictionary = ga.make_transition_dictionary(
        ga.aa_to_codon_extended, ga.allowed_aa_transitions
    )
    assert transition_dictionary["A"] == [
        "GGX",
        "GCX",
        "GTX",
        "TTR",
        "CTX",
        "YTR",
        "ATH",
    ]


def test_generate_swaptable():
    codon_to_codon_extended = ga.generate_swaptable(
        ga.codon_to_aa, ga.aa_to_codon_extended
    )
    assert codon_to_codon_extended["TTA"] == ("TTR", "CTX", "YTR")


def test_compare_letters():
    assert ga.compare_letters("H", "A", ga.ambiguity_code_to_nt_set) is True
    assert ga.compare_letters("A", "Y") is False
