"""The main function of Damian, DAMe's better version.

This script is the master of the Damian metabarcoding data pre-processing
tool. It takes in command line parameters, then figures out which submodule
to call and forward the data to that. The usage of this script is only from
the command line.
"""

import argparse
import logging

import sort


def parse_cl_arguments():
    """Parse command line arguments.

    Function to parser command line arguments.

    Returns
    -------
    namespace
        namespace with parsed arguments key-value pairs.

    """
    parser = argparse.ArgumentParser(prog="Damian", description="Damian \
                                     metabarcoding preprocessing tool",
                                     version="0.1", )
    subparser = parser.add_subparsers(help="Damian command to run",
                                      dest="command")

    # Add parser for the sorting step
    sortParser = subparser.add_parser("sort")
    sortParser.add_argument("-p", "--primers", help="File with forward and \
                            reverse primer sequence\n(Format: ForwardPrimer \
                            ReversePrimer)", required=False, default="",
                            metavar="PrimerFile")
    sortParser.add_argument("-p1", "--fwdPrimer", help="Sequence of forward \
                            primer", required=False, metavar="FwdPrimer")
    sortParser.add_argument("-p2", "--revPrimer", help="Sequence of reverse \
                            primer", required=False, metavar="RevPrimer")
    sortParser.add_argument("-t", "--tags", help="File with tag name and \
                            sequence \n(Format: TagName FwdTagSequence)",
                            required=True, metavar="TagFile")
    sortParser.add_argument("-s", "--sampleInfo", help="File with tag combo \
                            and pool for each sample\n(Format: Sample \
                            FwdTagName RevTagName PoolName)", required=True,
                            metavar="SampleInformationFile")
    sortParser.add_argument("-l", "--pool", help="File with pool information \
                            \n(Format: Poolname Read1Fastq [Read2Fastq])",
                            required=True, metavar="PoolInformationFile")
    sortParser.add_argument("-m", "--allowMultiplePrimers", help="Allow more \
                            one occurrance of the primer sequence in read. \
                            (Default False)", action="store_true",
                            dest="allow_multiple_primers")
    sortParser.add_argument("-pm", "--primerMismatches", help="Number of \
                            mismatches in primer. (Default 0)", type=int,
                            default=0, metavar="PrimerMismatches",
                            dest="primer_mismatches")
    sortParser.add_argument("-tm", "--tagMismatches", help="Number of allowed\
                            mismatches in tags. (Default 0)", type=int,
                            default=0, metavar="TagMismatches",
                            dest="tag_mismatches")
    sortParser.add_argument("-mo", "--merge_overlap", help="Merge read1 and \
                            read2 if overlapping by given number of bases or \
                            more (>=5) __NOT IMPLEMENTED YET__", type=int,
                            default=0, metavar="MinOverlap")
    sortParser.add_argument("-mm", "--merge_errors", help="Rate of mismatches \
                            allowed in overlap between reads: range [0,0.2]\
                            __NOT IMPLEMENTED YET__", type=float, default=0.0,
                            metavar="OverlapErrRate")
    sortParser.add_argument("-d", "--output_directory", help="Output \
                            directory. (Default: .)", default=".",
                            metavar="OutDirectory")
    sortParser.add_argument("-o", "--output_prefix", help="Prefix for output \
                            files. (Default : '')", default="",
                            metavar="OutPrefix", dest="output_prefix")

    # Add parser for the filter step
    filterParser = subparser.add_parser("filter")
    filterParser.add_argument("-i", "--infofile", help="Information file with\
                              fastq information for each sample.",
                              metavar="Infofile", required=True)
    filterParser.add_argument("-n", "--minNumPCRs", help="Minimum number of \
                              PCR replicates a sequence should be present in.",
                              metavar="minPCRNums", type=int, default=1,
                              required=False)
    filterParser.add_argument("-p", "--propPCRs", help="Minimum proportion of \
                              PCR replicates a sequence should be present in.",
                              metavar="propPCRs", type=float, default=1.0,
                              required=False)
    filterParser.add_argument("-m", "--minOccurence", help="Minimum number \
                              of times a sequence should be present, in a PCR \
                              replicate to be consider a true sequence.",
                              metavar="minTimes", type=int, default=1,
                              required=False)
    filterParser.add_argument("-l", "--minLength", help="Minimum length of \
                              the amplicon sequence - in case of single end \
                              or merged sequences, it is the length of the \
                              sequence, and in case of paired end reads, it \
                              is the sum of the length of the 2 reads.",
                              metavar="minLength", type=int, default=0,
                              required=False)
    filterParser.add_argument("-d", "--output_directory", help="Output \
                              directory", default=".", metavar="OutDirectory")
    filterParser.add_argument("-o", "--output_prefix", help="Prefix for \
                              output files", default="", metavar="OutPrefix",
                              dest="output_prefix")
    args = parser.parse_args()
    return(args)


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
        sorter = sort.sample_sorter(args, logger)
        sorter.read_tag_file(args.tags)
        if args.primers != "":
            sorter.read_primer_file(args.primers)
        elif args.fwdPrimer != "" and args.revPrimer != "":
            sorter.set_primer_seqs(args.fwdPrimer, args.revPrimer)
        else:
            raise IOError("Need to specify either a primer file, \
                           or primer sequences.")
        sorter.read_pool_file(args.pool)
        sorter.read_sample_information_file(args.sampleInfo)
        sorter.process_read_file()
    elif args.command == "filter":
        logger.info("Running the filter module.")


main()
