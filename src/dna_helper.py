"""The dna module for damian.

This module contains the dna classes for the damian application. It
includes, among other things, functions to translate ambiguity bases,
reverse complement, reverse, check if the dna strings are valid dna strings.
It also includes a definition of the class dna_pair, which is used to store
tag pair or primer pair.

Attributes
----------
DNA_UNAMBIGUOUS_CHARS : string
    String of upper and lower case chars when no ambiguity is allowed.
DNA_UNAMBIGUOUS_REV : string
    The reverse (dna wise) of the DNA_UNAMBIGUOUS_CHARS string.
DNA_AMBIGUOUS_CHARS: string
    String of upper and lower case chars, including ambiguity bases.
DNA_AMBIGUOUS_REV: string
    The reverse bases for the DNA_AMBIGUOUS_CHARS string.
COMPLEMENT_UNAMBIGOUS_DNA :
    The translation table for complementing dna sequences without ambiuity
    bases.
COMPLEMENT_AMBIGOUS_DNA :
    The translation table for complementing dna sequencing with ambiguity
    bases.
AMB_REGEX_DICT : dict
    The dictionary containing the corresponding regex string for each possible
    base, including ambiguity bases.

"""

import regex

DNA_UNAMBIGUOUS_CHARS = "ACGTUNacgtun"
DNA_UNAMBIGUOUS_REV = "TGCAANtgcaan"
DNA_AMBIGUOUS_CHARS = "ACGTURYSWKMBDHVNacgturyswkmbdhvn"
DNA_AMBIGUOUS_REV = "TGCAAYRWSMKVHDBNtgcaayrwsmkvhdbn"
# COMPLEMENT_UNAMBIGOUS_DNA = maketrans(DNA_UNAMBIGUOUS_CHARS,
#                                       DNA_UNAMBIGUOUS_REV)
# COMPLEMENT_AMBIGOUS_DNA = maketrans(DNA_AMBIGUOUS_CHARS, DNA_AMBIGUOUS_REV)

SUBSET = {"A": "A", "C": "C", "G": "G", "T": "T", "U": "U", "R": "AG",
          "Y": "CT", "S": "GC", "W": "AT", "K": "GT", "M": "AC", "B": "CGTYSK",
          "D": "AGTRWK", "V": "ACTYWM", "H": "ACGRSM", "N": "ACGTURYSWKMBDHVN"}
AMB_REGEX_DICT = {"A": r"A", "a": r"a", "C": r"C", "c": r"c", "G": r"G",
                  "g": r"g", "T": r"T", "t": r"t", "U": r"U", "u": r"u",
                  "R": r"[AG]", "r": r"[ag]", "Y": r"[CT]", "y": r"[ct]",
                  "S": r"[GC]", "s": r"[gc]", "W": r"[AT]", "w": r"[at]",
                  "K": r"[GT]", "k": r"[gt]", "M": r"[AC]", "m": r"[AC]",
                  "B": r"[CGT]", "b": r"[cgt]", "D": r"[AGT]", "d": r"[agt]",
                  "V": r"[ACT]", "v": r"[act]", "H": r"[ACG]", "h": r"[acg]",
                  "N": r"[ACGT]", "n": r"[acgt]"}


class dna_sequence(object):
    """Helper class to manipulate dna sequences.

    This class has helper functions to manipulate dna sequences.

    """

    @staticmethod
    def convert_ambiguity_to_regex(sequence, mismatches=0,
                                   preserve_case=False):
        """Convert the ambigous bases to regular expressions.

        This function takes a dna sequence with ambiguous bases in it, and
        converts the ambiuity bases to regular expressions, so that we can use
        these for string matching.

        Parameters
        ----------
        sequence : string
            dna sequence
        preserve_case : bool
            preserve case of input dna, if false, convert to upper case.

        Returns
        -------
        string
            dna sequence with modified ambiguity bases, ready for use in string
            matching.

        """
        if not preserve_case:   # convert to upper case if not preserve case
            sequence = sequence.upper()
        dna_regex = r""
        for c in sequence:
            dna_regex += AMB_REGEX_DICT[c]
        dna_regex = dna_regex + r"{s<="+mismatches+r"}"
        dna_regex = regex.compile(dna_regex)
        return dna_regex

    @staticmethod
    def find_first_subsequence(matcher_regex, target_dna):
        """Find first subsequnce of given sequence in an unambiguous sequence.

        Given a target sequence, find the first/best match of this target
        sequence in a given unambiguous dna sequence.

        Parameters
        ----------
        matcher_regex : compiled regex pattern
            compiled regular expression object
        target_dna : :obj: Bio.Seq.Seq object
            dna sequence object

        Returns
        -------
        int
            position of best match (first if many), -1 if no matches

        """
        matches = matcher_regex.search(target_dna)
        firstmatch = matches.group(0)
