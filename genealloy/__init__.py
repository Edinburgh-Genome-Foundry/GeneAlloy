from .codontable import (
    codon_to_aa,
    aa_to_codon_extended,
    codon_extended_to_aa,
    ambiguity_code_to_nt_set,
    allowed_aa_transitions,
    make_transition_dictionary,
    generate_swaptable,
    compare_letters,
)
from .genealloy import (
    convert_seq_to_codons,
    convert_codonlist_to_tuplelist,
    compare_then_get_letter_recursively,
    walk_seqstep,
    compare_sequence_tuplelists,
    compare_sequence_tuplelists_in_all_frames,
    find_partial_overlaps,
    make_genealloy,
    get_complement_tuplelist,
    get_reverse_tuplelist,
    get_reverse_complement_tuplelist,
    SeqStep,
)
from .version import __version__
