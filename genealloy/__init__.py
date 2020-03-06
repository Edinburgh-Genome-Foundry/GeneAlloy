from .codontable import (
    codon_to_aa,
    aa_to_codon_extended,
    codon_extended_to_aa,
    ambiguity_code_to_nt_set,
    generate_swaptable,
    compare_letters,
)
from .genealloy import (
    convert_seq_to_codons,
    convert_codonlist_to_tuplelist,
    compare_then_get_letter_recursively,
)
# from .codon_comparison_functions import (
#     find_matches,
#     report_results,
#     calculate_aa_solutions,
#     calculate_levenshtein_distances,
# )
from .version import __version__
