from .codontable import (
    codon_to_aa,
    aa_to_codon_extended,
    codon_extended_to_aa,
    complement_table,
    compare_letters,
)


def convert_seq_to_codons(seq):
    """Convert a string (sequence) into a list of 3-letter strings (triplets)."""
    seq_codons = [seq[i : i + 3] for i in range(0, len(seq), 3)]

    return seq_codons


def convert_codonlist_to_tuplelist(seq_codons, codon_to_codon_extended):
    """Convert a list of triplets into a list of tuples, using a swaptable.

    The swaptable is a dict of triplet: triplets, and determines the
    allowed swaps.
    """
    codon_extended = [None] * len(seq_codons)
    for i, codon in enumerate(seq_codons):
        codon_extended[i] = codon_to_codon_extended[codon]

    return codon_extended


def get_next_letter(sequence_tuplelist, current_letter_index):
    """Get next letter in the sequence."""
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
    """Get first letter in the next triplet."""
    letter_index = current_letter_index
    letter_index[1] = letter_index[1] + 1  # triplet position
    letter_index[2] = 0  # move letter position to triplet's first letter
    try:
        i, j, k = letter_index[0], letter_index[1], letter_index[2]
        letter = (sequence_tuplelist[i][j][k], letter_index)
        return letter
    except:
        raise


def compare_then_get_letter_recursively(
    host_tuplelist, parasite_tuplelist, host_letter, parasite_letter
):
    """Compare two letters then get next pair of letters recursively.

    Returns string for match or no match between the sequences.
    """

    # This function works only for in-frame comparisons and is given the first
    # letters of the host and parasite sequences. If the two letters match,
    # then gets the next parasite and host letters, and calls itself.
    # If can't get next parasite letter, then the comparison has finished
    # and a full match has been found.
    # If letters do not match, then gets the next parasite triplet and calls
    # itself; if there are no more parasite triplets, it gets the next host
    # triplet and calls itself. If there are no more host triplets, then there is
    # no match between the sequences.

    is_match = compare_letters(host_letter[0], parasite_letter[0])

    if is_match:
        try:
            next_parasite_letter = get_next_letter(
                parasite_tuplelist, parasite_letter[1]
            )

        except:
            return str(
                "Finished parasite sequence, match found! Ending host position:"
                + str(host_letter[1])
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


def walk_seqstep(seqstep):
    """Compare two sequences by calling `advance_step` until it returns the result."""
    while not seqstep.result:
        seqstep.advance_step()


def compare_sequence_tuplelists(parasite_tuplelist, host_tuplelist, frameshift):
    """Compare two sequence's tuplists for given frame and return list of matches."""
    len_parasite = len(parasite_tuplelist)
    len_host = len(host_tuplelist)
    list_of_matches = []

    for start_host_codon in range(0, len_host - len_parasite + 1):
        seqstep = SeqStep(
            host_tuplelist,
            parasite_tuplelist,
            start_host_codon=start_host_codon,
            frameshift=frameshift,
        )
        walk_seqstep(seqstep)
        if seqstep.match:
            list_of_matches.append(start_host_codon)

    return list_of_matches


def compare_sequence_tuplelists_in_all_frames(
    parasite_tuplelist, host_tuplelist, prefix=""
):
    """Compare two sequence's tuplists for all frames and return dict of matches."""
    results_for_all_frames = dict()
    for frameshift in [0, 1, 2]:
        result = compare_sequence_tuplelists(
            parasite_tuplelist, host_tuplelist, frameshift
        )
        key = prefix + str(frameshift)
        results_for_all_frames[key] = result

    return results_for_all_frames


def find_partial_overlaps(host, parasite, swaptable, verbose=True):
    flank = len(parasite) * "N"
    host_flank = flank + host + flank

    swaptable["NNN"] = ("NNN",)

    result = make_genealloy(host_flank, parasite, swaptable, verbose=True)

    return result


def make_genealloy(host, parasite, swaptable, verbose=True):
    """Compare two sequence strings and return dictionary of matches."""
    host_codons = convert_seq_to_codons(host)
    host_tuplelist = convert_codonlist_to_tuplelist(host_codons, swaptable)

    parasite_codons = convert_seq_to_codons(parasite)
    parasite_tuplelist = convert_codonlist_to_tuplelist(parasite_codons, swaptable)

    forward_results = compare_sequence_tuplelists_in_all_frames(
        parasite_tuplelist, host_tuplelist, prefix="f_"
    )

    reverse_complement_tuplelist = get_reverse_complement_tuplelist(host_tuplelist)

    reverse_complement_results = compare_sequence_tuplelists_in_all_frames(
        parasite_tuplelist, reverse_complement_tuplelist, prefix="rc_"
    )

    result = forward_results.copy()
    result.update(reverse_complement_results)

    if verbose:
        if all(value == [] for value in result.values()):
            print("These sequences cannot be mixed")
        else:
            print("A genealloy can be made using these sequences!")

    return result


def get_complement_tuplelist(codon_tuplelist):
    """Get complement triplets of a sequence tuplelist."""
    complement_tuplelist = []
    for index, codon in enumerate(codon_tuplelist):
        complement_tripletlist = []
        for triplet in codon:
            letter1 = triplet[0]
            letter2 = triplet[1]
            letter3 = triplet[2]

            complement_triplet = (
                complement_table[letter1]
                + complement_table[letter2]
                + complement_table[letter3]
            )
            complement_tripletlist.append(complement_triplet)
        complement_codon = tuple(complement_tripletlist)
        complement_tuplelist.append(complement_codon)

    return complement_tuplelist


def get_reverse_tuplelist(codon_tuplelist):
    """Get reverse of a tuplelist with reversed triplets."""
    reverse_tuplelist = []
    for codon in reversed(codon_tuplelist):
        reverse_tripletlist = []
        for triplet in codon:
            reverse_triplet = triplet[::-1]
            reverse_tripletlist.append(reverse_triplet)
        reverse_codon = tuple(reverse_tripletlist)
        reverse_tuplelist.append(reverse_codon)

    return reverse_tuplelist


def get_reverse_complement_tuplelist(codon_tuplelist):
    """Get reverse complement of a sequence's tuplelist."""
    complement_tuplelist = get_complement_tuplelist(codon_tuplelist)
    reverse_complement_tuplelist = get_reverse_tuplelist(complement_tuplelist)

    return reverse_complement_tuplelist


class Duodon:
    """Class for storing two triplets"""

    def __init__(self, first_triplet, second_triplet):
        self.first_triplet = first_triplet
        self.second_triplet = second_triplet


class SeqStep:
    """Class for keeping track of sequence comparison

    It stores a method that aligns a parasite triplet with two consecutive host
    triplets (duodons), a cursor that marks the position of the comparison process,
    and methods for generating duodons and comparing them with triplets.
    The `advance_step()` method attempts to advance the comparison by one codon
    step. It can (i) advance the cursor or (ii) conclude that there is no match
    between the sequences, or (iii) conclude that there is a match.

    Parameters
    ----------

    host_tuplelist
      A list of tuples. Each tuple stores the allowed triplets for a codon
      position of the host sequence.

    parasite_tuplelist
      A list of tuples. Each tuple stores the allowed triplets for a codon
      position of the parasite sequence.

    frameshift
      An integer (0, 1 or 2) denoting the frameshift between host and parasite.

    start_host_codon
      The host codon position from which the comparison should start.
    """

    def __init__(
        self, host_tuplelist, parasite_tuplelist, frameshift=0, start_host_codon=0
    ):
        self.host_tuplelist = host_tuplelist
        self.parasite_tuplelist = parasite_tuplelist
        self.frameshift = frameshift
        self.start_host_codon = start_host_codon
        self.cursor = 0
        self.len_parasite = len(parasite_tuplelist)
        self.len_host = len(host_tuplelist)
        self.parasite_path = []
        self.host_path = []
        self.match = False
        self.result = None

    def generate_duodons(self):
        self.cursor
        self.start_host_codon
        duodons = []
        host_codon = self.start_host_codon + self.cursor
        if self.host_path == []:  # first time it's made
            for triplet_1 in self.host_tuplelist[host_codon]:
                for triplet_2 in self.host_tuplelist[host_codon + 1]:
                    duodons.append(Duodon(triplet_1, triplet_2))
        else:
            triplet_1 = self.host_path[-1].second_triplet  # of last used duodon
            try:
                for triplet_2 in self.host_tuplelist[host_codon + 1]:
                    duodons.append(Duodon(triplet_1, triplet_2))
            except:
                duodons.append(Duodon(triplet_1, "NNN"))

        return duodons

    def return_all_matching_duodons(self, triplet, duodons, frameshift=0):
        matching_duodons = []
        for duodon in duodons:
            host_letters = duodon.first_triplet + duodon.second_triplet

            # compare 1st letter:
            if not compare_letters(host_letters[0 + frameshift], triplet[0]):
                continue
            if not compare_letters(host_letters[1 + frameshift], triplet[1]):
                continue
            if not compare_letters(host_letters[2 + frameshift], triplet[2]):
                continue

            matching_duodons.append(duodon)

        return matching_duodons

    def compare_triplet_and_duodon(self, triplet, duodon, frameshift=0):
        host_letters = duodon.first_triplet + duodon.second_triplet

        if not compare_letters(host_letters[0 + frameshift], triplet[0]):
            return False
        if not compare_letters(host_letters[1 + frameshift], triplet[1]):
            return False
        if not compare_letters(host_letters[2 + frameshift], triplet[2]):
            return False

        return True

    def advance_step(self):
        if self.result:
            return self.result

        if self.cursor == self.len_parasite:
            self.match = True
            self.result = (
                "Match found between parasite and host sequence. Start codon was: "
                + str(self.start_host_codon)
            )
            return self.result

        parasite_triplets = list(self.parasite_tuplelist[self.cursor])
        host_duodons = self.generate_duodons()

        while True:
            try:
                parasite_triplet = parasite_triplets[0]
            except:
                self.result = "Sequences don't match. Start codon was: " + str(
                    self.start_host_codon
                )
                return self.result
            else:
                host_duodons_for_parasite_triplet = host_duodons[:]
                host_duodons_for_parasite_triplet = self.return_all_matching_duodons(
                    parasite_triplet,
                    host_duodons_for_parasite_triplet,
                    frameshift=self.frameshift,
                )
                self.parasite_path.append(parasite_triplet)

                while True:
                    try:
                        host_doudon = host_duodons_for_parasite_triplet[0]
                    except:
                        self.parasite_path.pop()
                        del parasite_triplets[0]
                        break
                    else:
                        if self.compare_triplet_and_duodon(
                            parasite_triplet, host_doudon, frameshift=self.frameshift
                        ):
                            self.host_path.append(host_doudon)
                            self.cursor += 1
                            return "Codon matched, cursor advanced"
                        else:
                            del host_duodons_for_parasite_triplet[0]
