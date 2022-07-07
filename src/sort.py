"""Sorting module for damian.

The sorting module provides the functions required to do the sorting of the
input metabarcoding fastq files. In addition to performing the regular sorting
for merged AdapterRemoved sequences, it allows to sort paired end input data
with 2 output (read1 and read 2) for each sample. This is useful for using
DADA2 downstream - maybe.
"""

import gzip

from Bio import SeqIO
#from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

from dna_helper import dna_utility as DU


class sample_sorter():
    """Sample sorter class.

    Class contains methods to sort read files into the corresponding sample
    specific files. It can accomodate multiple input formats - single end
    reads, paired end reads and multiple output formats - paired end (as for
    DADA2), single end (if input is single end), and merged (if merge of paired
    end reads is requested).

    Attributes
    ----------
    merge_overlap : int
        minimum overlap in 3' end of read1 and 5' end of read2 for merging.
    merge_errors : float
        rate of mismatches allowed in overlap region of the paired end reads.
    tag_errors : float
        rate or number of mismatches tolerated in tags.
    primer_errors : float
        rate or number of mismatches tolerated in primers.
    output_prefix : string
        prefix for output files from sorting.
    output_directory : string
        directory for output files from sorting.

    """

    _tag_dict = {}
    _primer_pair = None
    _primers_rgx = None
    _samp_info = {}
    _pool_info = {}
    _merge = False
    # Primer type counts details
    # 0: no primer, 1: F-R good, 2: F-noR, 3:noF-R
    # 4: R-F, 5: R-noF, 6: noR-F, 7: ampliconempty
    # 8: F-R multi, 9: R-F multi
    _primer_type_counts = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    _tag_type_counts = [0, 0, 0, 0, 0, 0, 0, 0]
    merge = False
    merge_overlap = 10
    tag_errors = 0.0
    primer_errors = 0
    output_directory = ""
    output_prefix = ""
    logger = None

    def __init__(self, sorter_args, logger):
        """Construct a sorter class.

        Class to sort the read files to get sample wise files, using the tag
        and primer information. This constructor initializes all the values
        needed for running the sorting process.

        Parameters
        ----------
        sorter_args : namespace
            namespace from argparse module

        """
        self.logger = logger
        ###################################################
        # The next 2 options are not user since           #
        # they are not implemented yet, but will they     #
        # ever be.                                        #
        ###################################################
        self.merge_overlap = sorter_args.merge_overlap
        self.merge_errors = sorter_args.merge_errors
        ###################################################
        self. allowMultiPrims = sorter_args.allow_multiple_primers
        self.tag_errors = sorter_args.tag_mismatches
        self.primer_errors = sorter_args.primer_mismatches
        self.output_directory = sorter_args.output_directory
        self.output_prefix = sorter_args.output_prefix
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
        if self.tag_errors < 0:
            raise ValueError("Tag mismatch rate should be > 0.")
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
        fwd_primer = Seq(fwd_sequence)#, IUPAC.IUPACAmbiguousDNA())
        rev_primer = Seq(rev_sequence)#, IUPAC.IUPACAmbiguousDNA())
        self.logger.info("Setting foward primer to " + fwd_sequence)
        self.logger.info("Setting reverse primer to " + rev_sequence)
        self._primer_pair = (fwd_primer, rev_primer)

    def _conv_primers_regex(self):
        """Convert primer seqs to regexes."""
        self._primers_rgx = (DU.conv_ambig_regex(str(self._primer_pair[0]),
                                                 mismatches=self.primer_errors,
                                                 preserve_case=False),
                             DU.conv_ambig_regex(str(self._primer_pair[1]),
                                                 mismatches=self.primer_errors,
                                                 preserve_case=False))

    def __log_in_details(self):
        """Print the information about the sorter.

        Print detailed information on the dictionaries in this class.

        """
        self.logger.debug("Sorter details")
        self.logger.debug("--------------")
        self.logger.debug("Primers:")
        self.logger.debug("  Fwd: " + self._primer_pair[0])
        self.logger.debug("  Rev: " + self._primer_pair[1])
        self.logger.debug("  # mismatches allowed: " + str(self.primer_errors))
        self.logger.debug("Tags:")
        for tag_name, tag_seq in self._tag_dict.items():
            self.logger.debug("  " + tag_name + ": " + tag_seq)
        self.logger.debug("  # of mismatches allowed: " + str(self.tag_errors))
        self.logger.debug("Pools and samples:")
        for pool_name, fastqs in self._pool_info.items():
            log = "  " + pool_name + ": "
            log += " ".join([(x + ":" + str(y)) for x, y in fastqs.items()])
            self.logger.debug(log)
            cur_samples = self._samp_info[pool_name]
            for tag_pair, samp in cur_samples.items():
                log = ("    " + "\t".join(tag_pair) + "\t")
                log += ("\t".join([str(x) for x in samp]))
                self.logger.debug(log)
        self.logger.debug("Output details:")
        self.logger.debug("  Directory: " + self.output_directory)
        self.logger.debug("  Prefix: " + self.output_prefix)

    def __log_out_details(self):
        """Print primer and tag type count details.

        Print all the details on the primer types founds and the tag type
        found in the pool.
        """
        self._primer_type_counts[1] += self._primer_type_counts[8]
        self._primer_type_counts[4] += self._primer_type_counts[9]
        ptc = [str(x) for x in self._primer_type_counts]
        self.logger.info("Primer match type details")
        self.logger.info("-------------------------")
        self.logger.info("  Neither primer       :" + ptc[0])
        self.logger.info("  Both F-R primer      :" + ptc[1])
        self.logger.info("     with multi-primers:" + ptc[8])
        self.logger.info("  F primer, no R primer:" + ptc[2])
        self.logger.info("  no F primer, R primer:" + ptc[3])
        self.logger.info("  Both R-F primer      :" + ptc[4])
        self.logger.info("     with multi-primers:" + ptc[9])
        self.logger.info("  R primer, no F primer:" + ptc[5])
        self.logger.info("  no R primer, F primer:" + ptc[6])
        self.logger.info("  No barcode (SE only) :" + ptc[7])
        ttc = [str(x) for x in self._tag_type_counts]
        self.logger.info("Tag match type details")
        self.logger.info("----------------------")
        self.logger.info("  F-R primer orientation")
        self.logger.info("    Both tags            :" + ttc[0])
        self.logger.info("    no F tag, R tag found:" + ttc[1])
        self.logger.info("    F tag, no R tag found:" + ttc[2])
        self.logger.info("    Neither tag found    :" + ttc[3])
        self.logger.info("  R-F primer orientation")
        self.logger.info("    Both tags            :" + ttc[4])
        self.logger.info("    no F tag, R tag found:" + ttc[5])
        self.logger.info("    F tag, no R tag found:" + ttc[6])
        self.logger.info("    Neither tag found    :" + ttc[7])
        self.__reset_counts()

    def __reset_counts(self):
        """Reset counts of tag and primer matches."""
        for index in range(len(self._primer_type_counts)):
            self._primer_type_counts[index] = 0
        for index in range(len(self._tag_type_counts)):
            self._tag_type_counts[index] = 0

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
                self._conv_primers_regex()
            else:
                self.logger.error("More than one primer pair found in file.")
                raise IOError("Multiple primer pairs found in primer file.")
        primer_file.close()

    def read_tag_file(self, tag_file_name):
        """Process the tags file.

        Process the tag file, which has the format: TagName ForwardTag. Read
        the tag file, and store a dictionary that contains the tagname and a
        pair containing the forward and reverse tag, for all samples.

        Parameters
        ----------
        tag_file_name : string
            input tag file name - one sample per line, with format:
            TagName ForwardTag

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
            if len(tokens) != 2:
                self.logger.error("Line does not have the correct format.")
                raise IOError("Line:" + line +
                              "\ndoes not have the correct format.")
            if tokens[0] in self._tag_dict:
                self.logger.error("Repeat tag name: " + tokens[0])
                raise IOError(tokens[0] + " already present in file.")
            forwardTag = Seq(tokens[1])#, IUPAC.IUPACUnambiguousDNA())
            self._tag_dict[tokens[0]] = forwardTag
        self.logger.info("Read " + str(len(self._tag_dict)) + " valid tag \
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
                self.logger.error("Incorrect file format in sample file.")
                raise IOError("Incorrect number of tokens in sample file:" +
                              line)
            (sample, fwd_tagname, rev_tagname, pool_name) = toks
            if pool_name not in self._pool_info:
                self.logger.error("Pool name " + pool_name + " not found in " +
                                  "pool file.")
                raise KeyError("Pool name " + pool_name + " not found in " +
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
            if tag_pair in self._samp_info[pool_name]:
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
            self._samp_info[pool_name][tag_pair] = [sample, replicate_num]
        self.logger.info("Information for " + str(len(self._pool_info)) +
                         " pools in file.")
        self.logger.info("Read details for " + str(len(sample_count)) +
                         " samples.")
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
        self.logger.info("Read information for " + str(len(self._pool_info)) +
                         " pools.")
        pool_file.close()

    def process_read_file(self):
        """Process the read file to create per sample read files.

        Read the paired-end or single end fastq (possibly gzipped) file and
        process it to separate it into the different samples, as defined in the
        sample information, tag and pool information data in this class.

        """
        # First print out the details :)
        self.__log_in_details()
        # figure out for each pool, what is going on.
        for pool_name in self._pool_info:
            self.logger.debug("Processing reads for pool " + pool_name)
            read1_filename = self._pool_info[pool_name]["read1"]
            read2_filename = self._pool_info[pool_name]["read2"]
            is_gzipped = self._pool_info[pool_name]["zipped"]
            outname = self.output_directory + "/" + self.output_prefix
            outname += "_" + pool_name
            self.logger.debug("Read1 file:" + read1_filename)
            self.logger.debug("Read2 file:" + read2_filename)
            self.logger.debug("Zip status:" + str(is_gzipped))
            if read2_filename == "":
                self.logger.debug("Processing single end files.")
                haps = self.__process_single_end(read1_filename, is_gzipped)
                self.__write_out_files(haps, outname, pool_name,
                                       single_end=True)
                self.__log_out_details()
            else:
                self.logger.debug("Processing paired end files.")
                haps = self.__process_paired_end(read1_filename,
                                                 read2_filename, is_gzipped)
                self.__write_out_files(haps, outname, pool_name,
                                       single_end=False)
                self.__log_out_details()

    def __write_out_files(self, haps, outprefix, pool_name, single_end):
        """Write the amplicons to the output file.

        Write the amplicon information  to the output file.

        Parameters
        ----------
        haps : dict
            Dictionary with amplicon data
        outprefix : string
            Output file prefix
        pool_name : string
            Pool name
        single_end : bool
            Is this single end?

        """
        # Make a set of tags used in this pool.
        self.logger.debug("Writing output files.")
        pool_tags = set()
        pool_tag_pairs = []
        for tag_pair in self._samp_info[pool_name]:
            pool_tags.update(tag_pair)
            pool_tag_pairs.append(tag_pair)
        # Decode key :)
        # C => correct pair
        # B => Both used, but not the correct pair
        # F => Forward used, not reverse
        # R => Reverse used, not forward
        # N => Neither used
        summary_lines = {"C": "", "B": "", "F": "", "R": "", "N": ""}
        tag_out = open(outprefix + ".tagInfo", "w")
        if single_end:
            tag_out.write("FTag\tRTag\tSeq\tCount\tType\n")
        else:
            tag_out.write("FTag\tRTag\tFSeq\tRSeq\tCount\tType\n")
        for tag_pair, current_haps in haps.items():
            (ftag, rtag) = tag_pair
            ftag_in_pool = (ftag in pool_tags)
            rtag_in_pool = (rtag in pool_tags)
            tag_type = "N"
            if ftag_in_pool and rtag_in_pool:
                if tag_pair in pool_tag_pairs:
                    tag_type = "C"
                else:
                    tag_type = "B"
            elif ftag_in_pool:
                tag_type = "F"
            elif rtag_in_pool:
                tag_type = "R"
            total_seqs = 0
            header = ftag + "\t" + rtag + "\t"
            footer = "\t" + tag_type + "\n"
            # if single_end:
            #     for amplicon, amplicon_count in current_haps.items():
            #         tag_out.write(header + amplicon + "\t")
            #         tag_out.write(str(amplicon_count) + footer)
            #         total_seqs += amplicon_count
            # else:
            for amplicon, amplicon_count in current_haps.items():
                tag_out.write(header + amplicon + "\t")
                tag_out.write(str(amplicon_count) + footer)
                total_seqs += amplicon_count
            summary_lines[tag_type] += (header + str(len(current_haps)) + "\t")
            summary_lines[tag_type] += (str(total_seqs) + footer)
        tag_out.close()
        self.logger.debug("Finished writing output files.")
        # Write the summary file.
        self.__write_summary_file(outprefix, summary_lines)

    def __write_summary_file(self, outprefix, summary_lines):
        """Write the summary file.

        Write summary file based on the summary lines from the write out files
        function.

        Parameters
        ----------
        outprefix : string
            Prefix for output summary file
        summary_lines : dict
            Info on different types of output tag combinations

        """
        self.logger.debug("Writing summary file.")
        tag_summary = open(outprefix + ".summaryCounts", "w")
        tag_summary.write("Tag1\tTag2\tUniqSeqs\tTotalSeqs\tType\n")
        tag_summary.write("Correct combination of tags used in pool\n")
        tag_summary.write("----------------------------------------\n")
        tag_summary.write(summary_lines["C"])
        tag_summary.write("Both tags used in pool, but not this combination\n")
        tag_summary.write("------------------------------------------------\n")
        tag_summary.write(summary_lines["B"])
        tag_summary.write("Only forward tag used in pool, reverse not found\n")
        tag_summary.write("------------------------------------------------\n")
        tag_summary.write(summary_lines["F"])
        tag_summary.write("Only reverse tag used in pool, forward not found\n")
        tag_summary.write("------------------------------------------------\n")
        tag_summary.write(summary_lines["R"])
        tag_summary.write("Neither tag used in this pool\n")
        tag_summary.write("-----------------------------\n")
        tag_summary.write(summary_lines["N"])
        tag_summary.close()
        self.logger.debug("Finished writing summary file.")

    def __process_single_end(self, read_filename, is_gzipped):
        """Process single end read files.

        Read and separate single end read file into sample read files.

        Parameters
        ----------
        read_filename : string
            fastq name for read file
        is_gzipped : bool
            is input file gzipped?

        Returns
        -------
        hap : dict
            with amplicon and tag combination information

        """
        haps = {}
        if is_gzipped:
            infile = gzip.open(read_filename)
        else:
            infile = open(read_filename)
        for record in SeqIO.parse(infile, "fastq"):
            seq = record.seq
            seq_rc = record.seq.reverse_complement()
            (fstart, fend, rstart, rend, match_type) = self.__find_primer_pos(
                seq, seq_rc)
            self._primer_type_counts[match_type] += 1
            if (self.allowMultiPrims and match_type == 8) or match_type == 1:
                ftag = str(seq)[0:fstart]
                rtag = str(seq_rc)[0:rstart]
                amplicon = str(seq)[fend:(len(seq)-rend)]
                self.logger.debug("FT: " + str(fstart) + " " + ftag)
                self.logger.debug("RT: " + str(rstart) + " " + rtag)
                self.logger.debug("match: " + str(match_type))
            elif (self.allowMultiPrims and match_type == 9) or match_type == 4:
                ftag = str(seq_rc)[0:rstart]
                rtag = str(seq)[0:fstart]
                amplicon = str(seq_rc)[rend:(len(seq)-fend)]
                self.logger.debug("FT: " + str(rstart) + " " + ftag)
                self.logger.debug("RT: " + str(fstart) + " " + rtag)
                self.logger.debug("match: " + str(match_type))
            else:
                continue
            if len(amplicon) == 0:
                self._primer_type_counts[7] += 1
                continue
            best_ftag = self.__find_best_tag_match(ftag)
            best_rtag = self.__find_best_tag_match(rtag)
            tag_type = (best_ftag == "") + 2*(best_rtag == "")
            self._tag_type_counts[4*(match_type == 4) + tag_type] += 1
            # At least one of the two tags was not found :(
            if tag_type:
                continue
            # Both tags found
            tag_combo = (best_ftag, best_rtag)
            if tag_combo not in haps:
                haps[tag_combo] = {}
                haps[tag_combo][amplicon] = 0
            elif amplicon not in haps[tag_combo]:
                haps[tag_combo][amplicon] = 0
            haps[tag_combo][amplicon] += 1
        infile.close()
        return (haps)

    def __process_paired_end(self, read1_filename, read2_filename, is_gzipped):
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

        Returns
        -------
        haps : dict
            with amplicon and tag combination information

        """
        haps = {}
        if is_gzipped:
            infile1 = gzip.open(read1_filename)
            infile2 = gzip.open(read2_filename)
        else:
            infile1 = open(read1_filename)
            infile2 = open(read2_filename)
        for (record1, record2) in zip(SeqIO.parse(infile1, "fastq"),
                                      SeqIO.parse(infile2, "fastq")):
            seq1 = record1.seq
            seq2 = record2.seq
            (fstart, fend, rstart, rend, match_type) = self.__find_primer_pos(
                seq1, seq2)
            self._primer_type_counts[match_type] += 1
            if match_type == 1:
                ftag = str(seq1)[0:fstart]
                rtag = str(seq2)[0:rstart]
                amplicon = (str(seq1)[fend:] + "\t" + str(seq2)[rend:])
            elif match_type == 4:
                ftag = str(seq2)[0:rstart]
                rtag = str(seq1)[0:fstart]
                amplicon = (str(seq2)[rend:] + "\t" + str(seq1)[fend:])
            else:
                continue
            best_ftag = self.__find_best_tag_match(ftag)
            best_rtag = self.__find_best_tag_match(rtag)
            tag_type = (best_ftag == "") + 2*(best_rtag == "")
            self._tag_type_counts[4*(match_type == 4) + tag_type] += 1
            # At least one of the two tags was not found :(
            if tag_type:
                continue
            # Both tags found
            tag_combo = (best_ftag, best_rtag)
            if tag_combo not in haps:
                haps[tag_combo] = {}
                haps[tag_combo][amplicon] = 0
            elif amplicon not in haps[tag_combo]:
                haps[tag_combo][amplicon] = 0
            haps[tag_combo][amplicon] += 1
        infile1.close()
        infile2.close()
        return (haps)

    def __find_primer_pos(self, fq_read1, fq_read2):
        """Find primer positions for using the 5' and 3' sequences.

        Find the primer positions, where both times the primers are found
        starting at the 5' of the read. In case of single end reads, it is
        such that both primers (F and R) are searched for from the 5' end of
        read. For paired end, in read 1 we search for it from the 5' end and in
        read2 from the 3' end (since it is already complemented). In our case,
        since we give read1 and RC(read1) for single end, it works out that the
        procedure is the same.

        Parameters
        ----------
        fq_read1 : :obj: Bio seq object
            dna sequence from fastq file - read 1
        fq_read2 : :obj: Bio seq object
            dna sequence from fastq file - read 2 - or rc of read 1, if
            single end sequence

        Returns
        -------
        int
            Start of forward primer match
        int
            End of forward primer match
        int
            Start of reverse primer match
        int
            End of reverse primer match
        int
            Type of match

        """
        read1_seq = str(fq_read1)
        read2_seq = str(fq_read2)
        match_type = 0  # no primers found
        # First find the forward primer in read 1, and reverse in read 2.
        # F in read 1 and R' in read2
        (fstart, fend, num_fwd) = DU.find_first_match(self._primers_rgx[0],
                                                      read1_seq)
        (rstart, rend, num_rev) = DU.find_last_match(self._primers_rgx[1],
                                                     read2_seq)
        if fstart != -1 and rstart != -1:
            if num_fwd > 1 or num_rev > 1:
                match_type = 8
            else:
                match_type = 1
        elif fstart != -1 and rstart == -1:
            match_type = 2
        elif fstart == -1 and rstart != -1:
            match_type = 3
        if fstart != -1 or rstart != -1:
            return((fstart, fend, rstart, rend, match_type))

        # R in read 1 and F' in read2
        (rstart, rend, num_rev) = DU.find_first_match(self._primers_rgx[1],
                                                      read1_seq)
        (fstart, fend, num_fwd) = DU.find_last_match(self._primers_rgx[0],
                                                     read2_seq)
        if fstart != -1 and rstart != -1:
            if num_fwd > 1 or num_rev > 1:
                match_type = 9
            else:
                match_type = 4
        elif rstart != -1 and fstart == -1:
            match_type = 5
        elif rstart == -1 and fstart != -1:
            match_type = 6
        return((fstart, fend, rstart, rend, match_type))

    def __find_best_tag_match(self, tag):
        """Find the best tag for the given sequence.

        Function to find the best tag match given, the sequence.

        Parameters
        ----------
        tag : string
            tag dna sequence of a read

        Returns
        -------
        best_tag : string
            name of the best string match

        """
        # First figure out forward tag
        best_dist = len(tag) + 100
        best_tag = ""
        num_matches = 0
        for tag_name in self._tag_dict:
            tag_seq = self._tag_dict[tag_name]
            cur_dist = DU.find_hamming_distance(tag_seq, tag)
            if cur_dist < best_dist:
                num_matches = 1
                best_tag = tag_name
                best_dist = cur_dist
            elif cur_dist == best_dist:
                if len(tag_seq) > len(self._tag_dict[best_tag]):
                    best_tag = tag_name
                    num_matches = 1
                elif len(tag_seq) == len(self._tag_dict[best_tag]):
                    num_matches += 1
        if best_dist > self.tag_errors or num_matches > 1:
            best_tag = ""
        return(best_tag)
