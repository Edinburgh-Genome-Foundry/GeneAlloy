API reference
=============

.. autoclass:: genealloy.SeqStep
  :members:
  :inherited-members:

.. autoclass:: genealloy.Duodon
  :members:
  :inherited-members:


.. py:data:: genealloy.codon_to_aa
  :annotation: = { "TTT": "F", ...

  Codon to amino acid dictionary.

.. py:data:: genealloy.aa_to_codon_extended
  :annotation: = {"A": ["GCX"], ...

  Amino acid to codon dictionary, using extended nucleotide letters.

.. py:data:: genealloy.codon_extended_to_aa
  :annotation: = {"GCX": "A", ...

  Codon to amino acid dictionary, using extended nucleotide letters.

.. py:data:: genealloy.ambiguity_code_to_nt_set
  :annotation: = {"A": {"A"}, ...

  Extended nucleotide letter to nucleotide letter dictionary.

.. py:data:: genealloy.complement_table
  :annotation: = {"A": "T", ...

  Extended nucleotide letter to complement letter dictionary.

.. py:data:: genealloy.allowed_aa_transitions
  :annotation: = {"A": ["G", "A", "V", "L", "I"], ...

.. autofunction:: genealloy.make_transition_dictionary

.. autofunction:: genealloy.generate_swaptable

.. autofunction:: genealloy.compare_letters

.. autofunction:: genealloy.convert_seq_to_codons

.. autofunction:: genealloy.convert_codonlist_to_tuplelist

.. autofunction:: genealloy.compare_then_get_letter_recursively

.. autofunction:: genealloy.walk_seqstep

.. autofunction:: genealloy.compare_sequence_tuplelists

.. autofunction:: genealloy.compare_sequence_tuplelists_in_all_frames

.. autofunction:: genealloy.find_partial_overlaps

.. autofunction:: genealloy.make_genealloy

.. autofunction:: genealloy.get_complement_tuplelist

.. autofunction:: genealloy.get_reverse_tuplelist

.. autofunction:: genealloy.get_reverse_complement_tuplelist

