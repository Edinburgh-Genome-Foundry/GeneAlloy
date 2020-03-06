from .codontable import (
    codon_to_aa,
    aa_to_codon_extended,
    codon_extended_to_aa,
    compare_letters,
)


def convert_seq_to_codons(seq):
    seq_codons = [seq[i : i + 3] for i in range(0, len(seq), 3)]

    return seq_codons


def convert_codonlist_to_tuplelist(seq_codons, codon_to_codon_extended):
    codon_extended = [None] * len(seq_codons)
    for i, codon in enumerate(seq_codons):
        codon_extended[i] = codon_to_codon_extended[codon]

    return codon_extended


def get_next_letter(sequence_tuplelist, current_letter_index):
    letter_index = current_letter_index
    if current_letter_index[2] == 2:
        # last letter of triplet, advance to next codon:
        letter_index = [letter_index[0] + 1, 0, 0]
        try:
            i, j, k = letter_index[0], letter_index[1], letter_index[2]
            letter = (sequence_tuplelist[i][j][k], letter_index)
            return letter
        except:
            raise

    else:
        # letter is not the last in triplet, get the next one:
        letter_index[2] = letter_index[2] + 1
        i, j, k = letter_index[0], letter_index[1], letter_index[2]
        letter = (sequence_tuplelist[i][j][k], letter_index)
        return letter


def get_letter_in_next_triplet(sequence_tuplelist, current_letter_index):
    letter_index = current_letter_index
    letter_index[1] = letter_index[1] + 1
    letter_index[2] = 0  # move to triplet's first letter
    try:
        i, j, k = letter_index[0], letter_index[1], letter_index[2]
        letter = (sequence_tuplelist[i][j][k], letter_index)
        return letter
    except:
        raise


def compare_then_get_letter_recursively(
    host_tuplelist, parasite_tuplelist, host_letter, parasite_letter
):

    is_match = compare_letters(host_letter[0], parasite_letter[0])

    if is_match:
        try:
            next_parasite_letter = get_next_letter(
                parasite_tuplelist, parasite_letter[1]
            )

        except:
            return (
                "Finished parasite sequence, match found! Ending host position:",
                host_letter[1],
            )

        else:
            next_host_letter = get_next_letter(
                host_tuplelist, host_letter[1]
            )  # always OK
            return compare_then_get_letter_recursively(
                host_tuplelist,
                parasite_tuplelist,
                next_host_letter,
                next_parasite_letter,
            )

    else:  # letters do not match, move on to next parasite triplet
        try:
            next_parasite_letter = get_letter_in_next_triplet(
                parasite_tuplelist, parasite_letter[1]
            )

        except:  # no more parasite triplet, get next host
            try:
                next_host_letter = get_letter_in_next_triplet(
                    host_tuplelist, host_letter[1]
                )
            except:
                return "No match for this starting codon position"
            else:
                i = parasite_letter[1][0]  # same codon
                j, k = 0, 0  # reset parasite to first triplet, first letter
                next_parasite_letter = (parasite_tuplelist[i][j][k], [i, j, k])
                return compare_then_get_letter_recursively(
                    host_tuplelist,
                    parasite_tuplelist,
                    next_host_letter,
                    next_parasite_letter,
                )

        else:
            i, j = host_letter[1][0], host_letter[1][1]  # same codon, same triplet
            k = 0  # reset host letter to beginning of triplet because parasite was reset
            next_host_letter = (host_tuplelist[i][j][k], [i, j, k])
            return compare_then_get_letter_recursively(
                host_tuplelist,
                parasite_tuplelist,
                next_host_letter,
                next_parasite_letter,
            )
