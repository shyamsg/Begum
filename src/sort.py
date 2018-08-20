"""Sorting module for damian.

The sorting module provides the functions required to do the sorting of the
input metabarcoding fastq files. In addition to performing the regular sorting
for merged AdapterRemoved sequences, it allows to sort paired end input data
with 2 output (read1 and read 2) for each sample. This is useful for using
DADA2 downstream.
"""

import gzip
import logging

import helper


class sample_sorter():
    """Sample sorter class.

    Class contains methods to sort read files into the corresponding sample
    specific files. It can accomodate multiple input formats - single end
    reads, paired end reads and multiple output formats - paired end (as for
    DADA2), single end (if input is single end), and merged (if merge of paired
    end reads is requested).

    Attributes
    ----------
    _tag_dict : dict
        dictionary with tagname and the corresponding tags.
    _primer_pair : :obj: dna_pair
        forward and reverse primers in a dna_pair object.
    _merged : bool
        merge output reads (if paired)
    _is_paired : bool
        are input reads paired end?
    _is_gzipped : bool
        input reads are zipped?
    fastq1 : string
        filename of read 1 (or only file in single end reads).
    fastq2 : string
        filename of read 2 in paired end reads.
    merge_overlap : int
        minimum overlap in 3' end of read1 and 5' end of read2 for merging.
    merge_errors : float
        rate of mismatches allowed in overlap region of the paired end reads.
    tag_errors : float
        rate or number of mismatches tolerated in tags.
    complexity_bases : bool
        are complexity bases allowed in the tags?
    primer_errors : float
        rate or number of mismatches tolerated in primers.
    output_prefix : string
        prefix for output files from sorting.
    output_directory : string
        directory for output files from sorting.

    """

    _tag_dict = {}
    _primer_pair = None
    _merge = False
    _is_paired = False
    _is_gzipped = False
    fastq1 = ""
    fastq2 = ""
    merge = False
    merge_overlap = 10
    tag_errors = 0.0
    complexity_bases = False
    primer_errors = 0.0
    output_directory = ""
    output_prefix = ""
    logger = None

    def __init__(self, sorter_args):
        """Construct a sorter class.

        Class to sort the read files to get sample wise files, using the tag
        and primer information. This constructor initializes all the values
        needed for running the sorting process.

        Parameters
        ----------
        sorter_args : namespace
            namespace from argparse module

        """
        self.logger = logging.getLogger("sample_sorter")
        self.fastq1 = sorter_args.fastq1
        self.fastq2 = sorter_args.fastq2
        self.merge_overlap = sorter_args.merge_overlap
        self.merge_errors = sorter_args.merge_errors
        self.complexity_bases = sorter_args.complexity_bases
        self.tag_errors = sorter_args.tag_errors
        self.primer_errors = sorter_args.primer_errors
        self.input_values_check()

    def input_values_check(self):
        """Check input values for correctness.

        Make sure all the input values are in the correct ranges, and set the
        internal variables accordingly.

        Raises
        ------
        ValueError
            when input values are out of range.

        """
        if self.tag_errors < 0.0 or self.tag_errors > 0.2:
            raise ValueError("Tag mismatch rate should be in [0.0, 0.2].")
        if self.primer_errors < 0.0 or self.primer_errors > 0.5:
            raise ValueError("Primer mismatch rate should be in [0.0, 0.5].")
        if self.merge_errors < 0.0 or self.merge_errors > 0.2:
            raise ValueError("Overlap mismatch rate should be in [0.0, 0.2].")
        if self.merge_overlap > 0 and self.merge_overlap < 5:
            raise ValueError("Minimum overlap between reads for merging is 5.")
        else:
            self._merge = True
        if self.fastq2 != "":
            self._is_paired = True
        if self._merge and not self._is_paired:
            self.logger.info("Merge flags reset since data is not paired end.")
            self._merge = False
            self.merge_overlap = 0
            self.merge_errors = 0
        if self.fastq1[-3:] == ".gz":
            self._is_gzipped = True
            if self.fastq2[-3:] != ".gz":
                self.logger.info("Read 2 file is not zipped while read 1 file \
                                 is (according to extension). Please double \
                                 check and change extensions before \
                                 continuing.")
                raise ValueError("Read file extensions do not match (for \
                                 gzipness)")
        if self.fastq1[-3:] != ".gz" and self.fastq2[-3:] == ".gz":
            self.logger.info("Read 1 file is not zipped while read 2 file is \
                             (according to extension). Please double check \
                             and change extensions before continuing.")
            raise ValueError("Read file extensions do not match (for \
                                 gzipness)")

    def read_primer_file(self, primer_filename):
        """Process the primer file.

        Read the primer file and get a pair of primer sequences.

        Parameters
        ----------
        primer_filename : string
            the name of the file, with primers. Must be in format: primer_name
            forward_primer reverse_primer.

        Raises
        ------
        IOError
            when primer file is not found, or has incorrect format.

        """
        primer_file = open(primer_filename)
        for line in primer_file:
            line = line.strip()
            if line == "":
                continue
            tokens = line.strip().split()
            if len(tokens) != 3:
                raise IOError(line+" does not have 3 tokens.")
            if self._primer_pair is None:
                self._primer_pair = helper.dna_pair(tokens[1], tokens[2])
            else:
                raise IOError("Multiple primer pairs found in primer file.")

    def read_tag_file(self, tag_file_name):
        """Process the tags file.

        Process the tag file, which has the format: TagName ForwardTag
        ReverseTag. Read the tag file, and store a dictionary that contains
        the tagname and a pair containing the forward and reverse tag, for
        all samples.

        Parameters
        ----------
        tag_file_name : string
            input tag file name, with format: TagName ForwardTag ReverseTag,
            one sample per line.

        Raises
        ------
        IOError
            when input file format is incorrect.

        """
        tag_file = open(tag_file_name)
        for line in tag_file:
            line = line.strip()
            if len(line) == 0:
                continue
            tokens = line.split()
            if len(tokens) != 3:
                raise IOError("Line:" + line +
                              "\ndoes not have the correct format.")
            if tokens[0] in self._tag_dict:
                raise IOError(tokens[0] + " already present in file.")
            self._tag_dict[tokens[0]] = helper.dna_pair(tokens[1], tokens[2])
        tag_file.close()

    def process_read_file(self, output_directory, output_prefix):
        """Process the read file to create per sample read files.

        Read the paired-end or single end fastq (possible gzipped) file and
        process it to separate it into the different samples, as defined in the
        tag dictionary.

        Parameters
        ----------
        output_directory : string
            output directory path
        output_prefix : string
            output prefix for sample files

        """
        if self._is_paired:
            self._process_paired_end(output_directory, output_prefix)
        else:
            self._process_single_end(output_directory, output_prefix)

    def _process_single_end(self, output_directory, output_prefix):
        """Process single end read files.

        Read and separate single end read file into sample read files.

        Parameters
        ----------
        output_directory : string
            output directory path
        output_prefix : string
            output prefix for sample files

        """
        if self._is_gzipped:
            infile = gzip.open(self.fastq1)
        else:
            infile = open(self.fastq1)
        for line in infile:
            toks = line.strip().split()

        def _find_best_tag_match(self, sequence):
            """Find the best tag for the given sequence.

            Function to find the best tag match given, the sequence.

            Parameters
            ----------
            seqeuence : :obj: dna_string
                dna sequence of a read

            Returns
            -------
            best_tag_name : string
                name of the best string match, under the constraints
            best_tag_mismatch : int
                mismatch with the best tag sequence

            """
