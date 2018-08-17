"""The main function of Damian, DAMe's better version.

This script is the master of the Damian metabarcoding data pre-processing
tool. It takes in command line parameters, then figures out which submodule
to call and forward the data to that. The usage of this script is only from
the command line.
"""

import argparse

import sort

import logging


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
                            required=True, dest="fastq1")
    sortParser.add_argument("-2", "--fq2", help="Read 2 of paired end",
                            metavar="Fastq2", dest="fastq2", default="")
    sortParser.add_argument("-p", "--primers", help="File with forward and \
                            reverse primer information\n(Format: PrimerName \
                            ForwardPrimer ReversePrimer)", required=True,
                            metavar="PrimerFile")
    sortParser.add_argument("-t", "--tags", help="File with forward and \
                            reverse tags for each sample\n(Format: TagName \
                            ForwardTag ReverseTag)", required=True,
                            metavar="TagFile")
    sortParser.add_argument("-pm", "--primerMismatches", help="Rate of \
                            mismatches in primer: range [0,0.5]", type=float,
                            default=0.0, metavar="PrimerErrorRate",
                            dest="primer_errors")
    sortParser.add_argument("-tm", "--tagMismatches", help="Rate of \
                            mismatches in tags: range [0,0.2]", type=float,
                            default=0.0, metavar="TagErrorRate",
                            dest="tag_errors")
    sortParser.add_argument("-c", "--complexity_bases", help="Complexity \
                            bases present in tags", action="store_true")
    sortParser.add_argument("-mo", "--merge_overlap", help="Merge read1 and \
                            read2 if overlapping by given number of bases or \
                            more (>=5)", type=int, default=0,
                            metavar="MinOverlap")
    sortParser.add_argument("-mm", "--merge_errors", help="Rate of mismatches \
                            allowed in overlap between reads: range [0,0.2]",
                            type=float, default=0.0, metavar="OverlapErrRate")
    sortParser.add_argument("-d", "--output_directory", help="Output \
                            directory", default=".", metavar="OutDirectory")
    sortParser.add_argument("-o", "--output_prefix", help="Prefix for output \
                            files", default="", metavar="OutPrefix",
                            dest="output_prefix")

    # Add parser for the filter step
    filterParser = subparser.add_parser("filter")
    filterParser.add_argument("-i", "--infile", help="Input file fastqs",
                              required=True)

    args = parser.parse_args()
    return args


def main():
    """Call main function."""
    logger = logging.getLogger("main")
    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    # create formatter and add it to the handlers
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - \
                                  %(message)s')
    ch.setFormatter(formatter)
    # add the handlers to the logger
    logger.addHandler(ch)
    args = parse_cl_arguments()
    if args.command == "sort":
        logger.info("Running the sorting module.")
        sorter = sort.sample_sorter(args)
        sorter.read_tag_file(args.tags)
        sorter.read_primer_file(args.primers)
    elif args.command == "filter":
        logger.info("Running the filter module.")


main()
