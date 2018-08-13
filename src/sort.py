"""Sorting module for damian.

The sorting module provides the functions required to do the sorting of the
input metabarcoding fastq files. In addition to performing the regular sorting
for merged AdapterRemoved sequences, it allows to sort paired end input data
with 2 output (read1 and read 2) for each sample. This is useful for using
DADA2 downstream.
"""


def read_primer_file(primer_filename):
    """Process the primer file.

    Read the primer file and return a dictionary of primer pairs, with the
    primer name as key and the pair of primers as value.

    Parameters
    ----------
    primer_filename : string
        the name of the file, with primers. One primer pair per line, in the
        format: PrimerName ForwardPrimer ReversePrimer.

    Returns
    -------
    dict
        dictionary with primer name as key and primer pair (forward primer,
        reverse primer) as value.

    Raises
    ------
    IOError
        when primer file is not found, or has incorrect format.
    OtherError
        when an other error

    """
    primer_dict = {}
    try:
        primer_file = open(primer_filename)
        for line in primer_file:
            line = line.strip()
            if line == "":
                continue
            tokens = line.strip().split()
            if len(tokens) != 2:
                raise IOError("File format for primer file is incorrect.")
            primer_dict[tokens[0]] = (tokens[1], tokens[2])
    except Exception as e:
        raise e
