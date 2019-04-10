"""Sorting module for damian.

The sorting module provides the functions required to do the sorting of the
input metabarcoding fastq files. In addition to performing the regular sorting
for merged AdapterRemoved sequences, it allows to sort paired end input data
with 2 output (read1 and read 2) for each sample. This is useful for using
DADA2 downstream - maybe.
"""

import gzip
import logging

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq


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
    _primer_pair : list
        forward and reverse primers in a list.
    _merged : bool
        merge output reads (if paired)
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
    _primer_pair_rc = None
    _samp_info = {}
    _pool_info = {}
    _merge = False
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
        ###################################################
        # The next 2 options are not user since           #
        # they are not implemented                        #
        ###################################################
        self.merge_overlap = sorter_args.merge_overlap
        self.merge_errors = sorter_args.merge_errors
        ###################################################
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
        if self.primer_errors < 0:
            raise ValueError("Primer mismatches should be >= 0.")
        if self.merge_errors < 0.0 or self.merge_errors > 0.2:
            raise ValueError("Overlap mismatch rate should be in [0.0, 0.2].")
        if self.merge_overlap > 0 and self.merge_overlap < 5:
            raise ValueError("Minimum overlap between reads for merging is 5.")
        else:
            self._merge = True

    def set_primer_seqs(self, fwd_sequence, rev_sequence):
        """Set the primer sequences.

        Set the primer sequences from the given forward and reverse sequences.

        Parameters
        ----------
        fwd_sequence : string
            forward primer sequence - ambiguities allowed.
        rev_sequence : string
            reverse primer sequence - ambiguities allowed.

        """
        fwd_primer = Seq(fwd_sequence, IUPAC.IUPACAmbiguousDNA())
        rev_primer = Seq(rev_sequence, IUPAC.IUPACAmbiguousDNA())
        self.logger.info("Setting foward primer to " + fwd_sequence)
        self.logger.info("Setting reverse primer to " + rev_sequence)
        self._primer_pair = (fwd_primer, rev_primer)
        self._primer_pair_rc = (fwd_primer.reverse_complement(),
                                rev_primer.reverse_complement())

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
            if len(tokens) != 2:
                self.logger.error("Primer file does not have correct number \
                                  of tokens.")
                raise IOError(line+" does not have 2 tokens.")
            if self._primer_pair is None:
                self.set_primer_seqs(tokens[0], tokens[1])
            else:
                self.logger.error("More than one primer pair found in file.")
                raise IOError("Multiple primer pairs found in primer file.")
        primer_file.close()

    def read_tag_file(self, tag_file_name):
        """Process the tags file.

        Process the tag file, which has the format: TagName ForwardTag
        ReverseTag. Read the tag file, and store a dictionary that contains
        the tagname and a pair containing the forward and reverse tag, for
        all samples.

        Parameters
        ----------
        tag_file_name : string
            input tag file name, with format:
            TagName ForwardTag ReverseTag, one sample per line.

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
                self.logger.error("Repeat tag name: " + tokens[0])
                raise IOError(tokens[0] + " already present in file.")
            forwardTag = Seq(tokens[1], IUPAC.IUPACUnambiguousDNA())
            reverseTag = Seq(tokens[2], IUPAC.IUPACUnambiguousDNA())
            self._tag_dict[tokens[0]] = (forwardTag, reverseTag)
        self.logger.info("Read " + len(self._tag_dict) + " valid tag \
                         combinations.")
        tag_file.close()

    def read_sample_information_file(self, sample_info_filename):
        """Process the sample information file.

        Read and process the sample information about which tag combinations
        for each of the samples.

        Parameters
        ----------
        sample_info_filename : string
            file with sample and their tag combinations and the pool name

        Raises
        ------
        IOError
            Sample tag combination file has incorrect format
        KeyError
            Pool or tag names not in the pool info or tag info files.

        """
        self._samp_info = {}
        sample_count = {}
        samp_file = open(sample_info_filename)
        for line in samp_file:
            toks = line.strip().split()
            if len(toks) != 4:
                raise IOError("Incorrect number of tokens in sample file:" +
                              line)
            (sample, fwd_tagname, rev_tagname, pool_name) = toks
            if pool_name not in self._pool_info:
                self.logger.error("Pool name " + pool_name + "not found in " +
                                  "pool file.")
                raise KeyError("Pool name " + pool_name + "not found in " +
                               "pool file.")
            # If this sample is in a new pool, then make that key
            if pool_name not in self._samp_info:
                self._samp_info[pool_name] = {}
            # check if tags names exist
            if fwd_tagname not in self._tag_dict:
                self.logger.error("Tag name " + fwd_tagname + " not found in" +
                                  " tag file.")
                raise KeyError("Tag name " + fwd_tagname + " not found in " +
                               "tag file.")
            if rev_tagname not in self._tag_dict:
                self.logger.error("Tag name " + rev_tagname + " not found in" +
                                  " tag file.")
                raise KeyError("Tag name " + rev_tagname + " not found in " +
                               "tag file.")
            # check if this tag combo already used in this pool
            tag_pair = (fwd_tagname, rev_tagname)
            if tag_pair in self._pool_info[pool_name]:
                self.logger.error("Tag combination " + tag_pair + " already " +
                                  "used in pool " + pool_name)
                raise ValueError("Tag combination " + tag_pair + " already " +
                                 "used in pool " + pool_name)
            # update replicate number
            if sample not in sample_count:
                sample_count[sample] = 1
            else:
                sample_count[sample] += 1
            replicate_num = sample_count[sample]
            # update pool info
            self._pool_info[pool_name][tag_pair] = [sample, replicate_num]
        samp_file.close()

    def read_pool_file(self, pool_info_filename):
        """Process the pool information file.

        Read and process the pool information file with pool names and their
        respective fastq files.

        Parameters
        ----------
        pool_info_filename : string
            file with pool name and fastq file(s)

        Raises
        ------
        IOError
            when file format is incorrect
        KeyError
            when duplicate pool names
        ValueError
            when mixture of paired and single end reads across pools

        """
        pool_file = open(pool_info_filename)
        for line in pool_file:
            toks = line.strip().split()
            if len(toks) < 2 or len(toks) > 3:
                self.logger.error("Pool information file has incorrect number"
                                  + " of columns on line: " + line)
                raise IOError("Pool information file has incorrect number"
                              + " of columns on line: " + line)
            if toks[0] in self._pool_info:
                self.logger.error("Pool name " + toks[0] + " already exists.")
                raise KeyError("Pool name " + toks[0] + " already exists.")
            if len(toks) == 2:  # single end reads
                self._pool_info[toks[0]] = {"read1": toks[1],
                                            "read2": "",
                                            "zipped": (toks[1][-3:] == ".gz")}
            elif len(toks) == 3:  # paried end reads
                first_gzipped = (toks[1][-3:] == ".gz")
                second_gzipped = (toks[2][-3:] == ".gz")
                if first_gzipped ^ second_gzipped:
                    self.logger.error("Only one of the 2 fastq files seems " +
                                      "be gzipped:" + toks[1] + " and " +
                                      toks[2])
                    raise IOError("Only one of the 2 fastq files seems be " +
                                  "gzipped:" + toks[1] + " and " + toks[2])
                self._pool_info[toks[0]] = {"read1": toks[1],
                                            "read2": toks[2],
                                            "zipped": first_gzipped}
        # check that all the reads are the same type, across the pools
        # either all single end or all paired end.
        paired_status = [(self._pool_info[pool_name]["read2"] != "")
                         for pool_name in self._pool_info]
        if not all(x == paired_status[0] for x in paired_status):
            self.logger.error("A combination of paired end and single end \
                               reads are give in the pool information file. \
                               Currently, we only allow one type of read \
                               across all pools.")
            raise ValueError("A combination of paired end and single end reads \
                              are give in the pool information file. We only \
                              allow one type of read across all pools.")
        pool_file.close()

    def process_read_file(self, output_directory, output_prefix):
        """Process the read file to create per sample read files.

        Read the paired-end or single end fastq (possibly gzipped) file and
        process it to separate it into the different samples, as defined in the
        sample information, tag and pool information data in this class.

        Parameters
        ----------
        output_directory : string
            output directory path
        output_prefix : string
            output prefix for sample files

        """
        # figure out for each pool, what is going on.
        for pool_name in self._pool_info:
            self.logger.debug("Processing reads for pool " + pool_name)
            read1_filename = self._pool_info[pool_name]["read1"]
            read2_filename = self._pool_info[pool_name]["read2"]
            is_gzipped = self._pool_info[pool_name]["zipped"]
            if read2_filename == "":
                self.__process_single_end(read1_filename, is_gzipped,
                                          output_directory, output_prefix)
            else:
                self.__process_paired_end(read1_filename, read2_filename,
                                          is_gzipped, output_directory,
                                          output_prefix)

    def __process_single_end(self, read_filename, is_gzipped, output_directory,
                             output_prefix):
        """Process single end read files.

        Read and separate single end read file into sample read files.

        Parameters
        ----------
        read_filename : string
            fastq name for read file
        is_gzipped : bool
            is input file gzipped?
        output_directory : string
            output directory path
        output_prefix : string
            output prefix for sample files

        """
        if is_gzipped:
            infile = gzip.open(read_filename)
        else:
            infile = open(read_filename)
        for record in SeqIO.parse(infile, "fastq"):
            (best_fprimer_pos, best_fprimer_mis, best_rprimer_pos,
             best_rprimer_mis) = self.__find_primer_positions(record)
            (bestftag, bestrtag) = self.__find_best_tag_match(ftag, rtag)

    def __process_paired_end(self, read1_filename, read2_filename, is_gzipped,
                             output_directory, output_prefix):
        """Process paired end read files.

        Read and separate single end read file into sample read files.

        Parameters
        ----------
        read1_filename : string
            fastq_file for the first read
        read2_filename : string
            fastq file for the second read
        is_gzipped : bool
            are the inputs gzipped?
        output_directory : string
            output directory path
        output_prefix : string
            output prefix for sample files

        """
        if is_gzipped:
            infile1 = gzip.open(read1_filename)
            infile2 = gzip.open(read2_filename)
        else:
            infile1 = open(read1_filename)
            infile2 = open(read2_filename)
        counter = 0
        for (line1, line2) in zip(infile1, infile2):
            counter += 1
            if counter % 4 != 1:
                continue
            
        pass

    def __find_primer_positions(self, fq_read1, fq_read2):
        """Find primer positions for using the 5' and 3' sequences.

        Parameters
        ----------
        fq_read1 : :obj: seqio record
            dna sequence from fastq file - read 1
        fq_read2 : :obj: seqio record
            dna sequence from fastq file - read 2 - or reverse of read 1, if
            single end sequence

        Returns
        -------
        best_fprimer_position
            Position on the sequence that is matching the forward primer
        best_fprimer_mismatches
            Number of mismatches on the sequence with best forward primer match
        best_rprimer_position
            Position on the sequence that is matching the reverse primer_errors
        best_rprimer_mismatches
            Number of mismatches on the sequence with best reverse primer match

        """
        # First find the forward primer in read 1, and reverse in read 2.
        forwar_primer_match

    def __find_best_tag_match(self, ftag, rtag):
        """Find the best tag for the given sequence.

        Function to find the best tag match given, the sequence.

        Parameters
        ----------
        ftag : string
            forward tag dna sequence of a read
        rtag : string
            reverse tag dna sequence of a read

        Returns
        -------
        best_tag_name : string
            name of the best string match, under the constraints
        best_tag_mismatch : int
            mismatch with the best tag sequence

        """
        pass
