"""The helper module for damian.

This module contains all the helper function for the damian application. It
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

from string import maketrans

DNA_UNAMBIGUOUS_CHARS = "ACGTUNacgtun"
DNA_UNAMBIGUOUS_REV = "TGCAANtgcaan"
DNA_AMBIGUOUS_CHARS = "ACGTURYSWKMBDHVNacgturyswkmbdhvn"
DNA_AMBIGUOUS_REV = "TGCAAYRWSMKVHDBNtgcaayrwsmkvhdbn"
COMPLEMENT_UNAMBIGOUS_DNA = maketrans(DNA_UNAMBIGUOUS_CHARS,
                                      DNA_UNAMBIGUOUS_REV)
COMPLEMENT_AMBIGOUS_DNA = maketrans(DNA_AMBIGUOUS_CHARS, DNA_AMBIGUOUS_REV)

AMB_REGEX_DICT = {"A": "A", "a": "a", "C": "C", "c": "c", "G": "G", "g": "g",
                  "T": "T", "t": "t", "U": "U", "u": "u", "R": "[AG]",
                  "r": "[ag]", "Y": "[CT]", "y": "[ct]", "S": "[GC]",
                  "W": "[AT]", "w": "[at]", "K": "[GT]", "k": "[gt]",
                  "M": "[AC]", "m": "[AC]", "B": "[CGT]", "b": "[cgt]",
                  "D": "[AGT]", "d": "[agt]", "V": "[ACT]", "v": "[act]",
                  "H": "[ACG]", "h": "[acg]", "N": "[ACGT]", "n": "[acgt]"}


def check_dna(dna_string, is_ambiguous=True):
    """Check validity of dna stringself.

    This function checks if the dna string is a valid dna string, by checking
    if the dna string has characters other than ACGTUN (if unambiguous) and
    ACGTUNRYSWKMBDHVN (if ambiguous).

    Parameters
    ----------
    dna_string : string
        the dna string to be validated.

    Returns
    -------
    bool
        true if dna_string is a valid dna string, false otherwise

    Raises
    ------
    No exceptions are raised in this function.

    """
    dna_string = dna_string.upper()
    if is_ambiguous:    # This checks if ambiguity is allowed.
        return (len(dna_string.translate(None, DNA_AMBIGUOUS_CHARS)) == 0)
    else:
        return (len(dna_string.translate(None, DNA_UNAMBIGUOUS_CHARS)) == 0)


def reverse_complement_dna(dna_string, is_ambiguous=True, preserve_case=False):
    """Reverse complement given dna string.

    This function takes a dna string and returns its reverse complement. Also,
    if the dna sequence has small letters, it converts them to upper case,
    unless the preserve case flag is set to True.

    Parameters
    ----------
    dna_string : string
        the dna sequence to be reverse complemented.
    is_ambiguous : bool
        does the dna sequence contain ambiguous bases?
    preserve_case : bool
        keep the case of the input intact, if false, convert to upper case.

    Returns
    -------
    string
        reverse complemented dna sequence as a string.

    Raises
    ------
    ValueError
        when the dna sequence is an invalid dna sequence.

    """
    if not preserve_case:   # if preserve case is off, then upper case the dna
        dna_string = dna_string.upper()
    if not check_dna(dna_string, is_ambiguous):     # check if dna is valid
        raise ValueError(dna_string+"is not a valid dna string.")
    if is_ambiguous:
        return (dna_string.translate(COMPLEMENT_AMBIGOUS_DNA))
    else:
        return (dna_string.translate(COMPLEMENT_UNAMBIGOUS_DNA))


def convert_ambiguity_to_regex(dna_string, preserve_case=False):
    """Convert the ambigous bases to regular expressions.

    This function takes a dna sequence with ambiguous bases in it, and converts
    the ambiuity bases to regular expressions, so that we can use these for
    string matching.

    Parameters
    ----------
    dna_string : string
        the input dna string that needs regex modification.
    preserve_case : bool
        preserve case of input dna, if false, convert to upper case.

    Returns
    -------
    string
        dna sequence with modified ambiguity bases, ready for use in string
        matching.

    Raises
    ------
    ValueError
        when the input dna_string is not valid dna

    """
    if not check_dna(dna_string, is_ambiguous=True):    # check dna validity
        raise ValueError(dna_string+"is not a valid dna string.")
    if not preserve_case:   # convert to upper case if not preserve case
        dna_string = dna_string.upper()
    dna_regex = ""
    for c in dna_string:
        dna_regex += AMB_REGEX_DICT[c]
    return dna_regex


class dna_pair():
    """Class to store pair of dna sequences, either primer pair or tag pair.

    This class is used to represent pairs of dna sequences. It stores a pair of
    dna sequences. It can return the sequences and their reverse complements on
    demand. *Remember* that both sequences are stored in 5' - 3' orientation.

    """

    def __init__(self, dna1, dna2):
        """Initialize object from dna_pair with the given dna sequences.

        Parameters
        ----------
        dna1 :  string
            first dna sequence
        dna2 : string
            second dna sequence

        Raises
        ------
        ValueError
            when input dna sequences are not valid.

        """
        self.dna1 = dna1
        self.dna2 = dna2

    @property
    def dna1(self):
        """DNA 1 sequence.

        The dna1 sequence should be a valid dna sequence. It will be checked at
        time of setting it. If not a valid sequence, a ValueError will be
        raised in the setter.
        """
        return self._dna1

    @dna1.setter
    def dna1(self, dna1):
        if check_dna(dna1):
            self._dna1 = dna1
        else:
            raise ValueError(dna1+" is not a valid dna sequence.")

    @property
    def dna2(self):
        """DNA 2 sequence.

        The dna2 sequence should be a valid dna sequence. It will be checked at
        time of setting it. If not a valid sequence, a ValueError will be
        raised in the setter.
        """
        return self._dna2

    @dna2.setter
    def dna2(self, dna2):
        if check_dna(dna2):
            self._dna2 = dna2
        else:
            raise ValueError(dna2+" is not a valid dna sequence.")

    @property
    def dna1_reverse_complement(self):
        """Reverse complement of dna1."""
        return reverse_complement_dna(self.dna1, is_ambiguous=True,
                                      preserve_case=False)

    @property
    def dna2_reverse_complement(self):
        """Reverse complement of dna2."""
        return reverse_complement_dna(self.dna2, is_ambiguous=True,
                                      preserve_case=False)
