GeneAlloy
=========

[![Build Status](https://github.com/Edinburgh-Genome-Foundry/genealloy/actions/workflows/build.yml/badge.svg)](https://github.com/Edinburgh-Genome-Foundry/genealloy/actions/workflows/build.yml)[![Coverage Status](https://coveralls.io/repos/github/Edinburgh-Genome-Foundry/GeneAlloy/badge.svg?branch=master)](https://coveralls.io/github/Edinburgh-Genome-Foundry/GeneAlloy?branch=master)

<p align="center">
<img alt="GeneAlloy logo" title="GeneAlloy" src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/GeneAlloy/master/logo/genealloy.png" width="150">
</p>

**Genealloy** helps designing *overlapping* sequences.

It takes two amino acid coding nucleotide sequences and a codon conversion table of allowed triplet -> triplet transitions, and determines whether one sequence can be inserted into the other one. Note that the package is **under development.**

*Overlapping sequences* are nucleotide sequences that encode different amino acid sequences on the same DNA or RNA region. These sequences are either on the complementary strands (in any frame), or on the same strand as frameshift sequences. This phenomenon is made possible by the redundancy of the genetic code (codon degeneracy).

In the metallurgic terminology used at the genome foundries, the host sequence (into which the shorter sequence is inserted) is called the *matrix* or *solvent,* and the shorter guest (or parasite) is called the *solute;* and a combination sequence is called a *genealloy.*

Install
-------

```bash
pip install genealloy
```

Usage
-----

```python
import genealloy as ga
swaptable = ga.generate_swaptable(ga.codon_to_aa, ga.aa_to_codon_extended)
host = 'TCGTCGTACCAGCCGCAGAGGAGAGCTACTTTT'
parasite =  'GTACCCGCTGCG'  # frameshift 2
ga.make_genealloy(host, parasite, swaptable)
```

Find partial overlaps:

```python
ga.find_partial_overlaps(host, parasite, swaptable)
```

Version
-------

The GeneAlloy project uses the [semantic versioning](https://semver.org) scheme. The package is **under development.**

License = MIT
-------------

Genealloy is [free software](https://www.gnu.org/philosophy/free-sw.en.html), which means the users have the freedom to run, copy, distribute, study, change and improve the software.

Genealloy was written at the [Edinburgh Genome Foundry](https://edinburgh-genome-foundry.github.io/) by [Peter Vegh](https://github.com/veghp) and is released under the MIT license.
