"""The main function of Damian, DAMe's better version.

This script is the master of the Damian metabarcoding data pre-processing
tool. It takes in command line parameters, then figures out which submodule
to call and forward the data to that. The usage of this script is only from
the command line.
"""

import sys

import argparse

# import sort

# import filter


def parse_cl_arguments():
    """Parse command line arguments.

    Function to parser command line arguments.

    Parameters
    ----------
    args_array : array_like
        input array with argument tokens.

    Returns
    -------
    namespace
        namespace with parsed arguments key-value pairs.

    Raises
    ------
    Error
        when required arguments are not passed, or there are superfluous args.

    """
    parser = argparse.ArgumentParser(prog="Damian", description="Damian \
                                     metabarcoding preprocessing tool",
                                     version="0.1", )
    subparser = parser.add_subparsers(help="Damian command to run",
                                      dest="command")

    # Add parser for the sorting step
    sortParser = subparser.add_parser("sort")
    sortParser.add_argument("-1", "--fq1", help="Merged fastq file or \
                            read 1 of paired end", metavar="Fastq1",
                            required=True)
    sortParser.add_argument("-2", "--fq2", help="Read 2 of paired end",
                            metavar="Fastq2", required=True)
    sortParser.add_argument("-p", "--primers", help="File with forward and \
                            reverse primer information\n(Format: PrimerName \
                            ForwardPrimer ReversePrimer)", required=True,
                            metavar="PrimerFile")
    sortParser.add_argument("-t", "--tags", help="File with forward and \
                            reverse tags for each sample\n(Format: TagName \
                            ForwardTag ReverseTag)", required=True,
                            metavar="TagFile")
    sortParser.add_argument("-pm", "--primerMismatches", help="Rate (if < 1) \
                            or number (if >= 1) of mismatches in primer",
                            type=float, required=False, default=1,
                            metavar="#PrimerErrors")
    sortParser.add_argument("-tm", "--tagMismatches", help="Rate (if < 1) or \
                            number (if >= 1) of mismatches in tags",
                            type=float, required=False, default=1,
                            metavar="#TagErrors")

    # Add parser for the filter step
    filterParser = subparser.add_parser("filter")
    filterParser.add_argument("-i", "--infile", help="Input file fastqs",
                              required=True)

    args = parser.parse_args()
    return args


def main():
    """Call main function."""
    args = parse_cl_arguments()
    print args.command


main()
