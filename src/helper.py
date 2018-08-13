"""The helper module for damian.

The module contaning all the helper function for the damian module. It
includes, among other things, functions to translate ambiguity bases,
reverse complement, reverse, check if the dna strings are valid dna strings.
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
