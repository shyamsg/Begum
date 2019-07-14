"""Filtering module for damian.

The filtering module provides the functions required to filter the sorted reads
or read pairs by length, number of occurrences, % of PCRs they show up in, etc.
It is mostly just a simple read the output of sorting and apply simple filters.
"""


class filter_sorted():
    """Filter class.

    Class contains methods for filtering the output of the sort function from
    Damian.

    Attributes
    ----------
    input_file : string
        output file from sort
    output_prefix : string
        prefix for output files from sorting.
    output_directory : string
        directory for output files from sorting.

    """

    prop_pcr = 0
    min_count = 0
    min_length = 0
    input_prefix = ""
    output_directory = ""
    output_prefix = ""
    logger = None
    _samp_info = None
    _num_pools = 0
    _haps_info = {}
    _rep_info = {}

    def __init__(self, filter_args, logger):
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
        self.prop_pcr = filter_args.propPCRs
        self.min_count = filter_args.minOccurence
        self.min_length = filter_args.minLength
        self.input_prefix = filter_args.inputPrefix
        self.output_directory = filter_args.output_directory
        self.output_prefix = filter_args.output_prefix
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
        if self.prob_pcr < 0 or self.prob_pcr > 1:
            raise ValueError("PCR proportion should be between 0 and 1.")
        if self.min_count < 1:
            raise ValueError("Minimum read count per PCR should be >= 1.")
        if self.min_length < 1:
            raise ValueError("Minimum length should be greater than 0.")

    def __log_in_details(self):
        """Print the information about the filter.

        Print detailed information on the members of this class.

        """
        self.logger.debug("Filter details")
        self.logger.debug("--------------")
        self.logger.debug("Options:")
        self.logger.debug("  Proportion of PCR:   " + str(self.prop_pcr))
        self.logger.debug("  Mininum # of reads:  " + str(self.min_count))
        self.logger.debug("  Minimum read length: " + str(self.min_length))
        self.logger.debug("  Input name prefix:   " + str(self.input_prefix))
        self.logger.debug("Output details:")
        self.logger.debug("  Directory: " + self.output_directory)
        self.logger.debug("  Prefix: " + self.output_prefix)

    # def __log_out_details(self):
    #     """Print primer and tag type count details.
    #
    #     Print all the details on the primer types founds and the tag type
    #     found in the pool.
    #     """
    #     self._primer_type_counts[1] += self._primer_type_counts[8]
    #     self._primer_type_counts[4] += self._primer_type_counts[9]
    #     ptc = [str(x) for x in self._primer_type_counts]
    #     self.logger.info("Primer match type details")
    #     self.logger.info("-------------------------")
    #     self.logger.info("  Neither primer       :" + ptc[0])
    #     self.logger.info("  Both F-R primer      :" + ptc[1])
    #     self.logger.info("     with multi-primers:" + ptc[8])
    #     self.logger.info("  F primer, no R primer:" + ptc[2])
    #     self.logger.info("  no F primer, R primer:" + ptc[3])
    #     self.logger.info("  Both R-F primer      :" + ptc[4])
    #     self.logger.info("     with multi-primers:" + ptc[9])
    #     self.logger.info("  R primer, no F primer:" + ptc[5])
    #     self.logger.info("  no R primer, F primer:" + ptc[6])
    #     self.logger.info("  No barcode (SE only) :" + ptc[7])
    #     ttc = [str(x) for x in self._tag_type_counts]
    #     self.logger.info("Tag match type details")
    #     self.logger.info("----------------------")
    #     self.logger.info("  F-R primer orientation")
    #     self.logger.info("    Both tags            :" + ttc[0])
    #     self.logger.info("    no F tag, R tag found:" + ttc[1])
    #     self.logger.info("    F tag, no R tag found:" + ttc[2])
    #     self.logger.info("    Neither tag found    :" + ttc[3])
    #     self.logger.info("  R-F primer orientation")
    #     self.logger.info("    Both tags            :" + ttc[4])
    #     self.logger.info("    no F tag, R tag found:" + ttc[5])
    #     self.logger.info("    F tag, no R tag found:" + ttc[6])
    #     self.logger.info("    Neither tag found    :" + ttc[7])
    #     self.__reset_counts()

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
            # If this sample is in a new pool, then make that key
            if pool_name not in self._samp_info:
                self._samp_info[pool_name] = {}
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
        self.logger.info("Information for " + str(len(self._samp_info)) +
                         " pools in file.")
        self.logger.info("Read details for " + str(len(sample_count)) +
                         " samples.")
        self._num_pools = len(self._samp_info)
        samp_file.close()
        self.__initialize_rep_info()

    def __initialize_rep_info(self):
        """Initialize the haps info dictionary."""
        for pool in self._samp_info:
            cur_pool = self._samp_info[pool]
            for tag_pair in cur_pool:
                (sample, replicate) = cur_pool[tag_pair]
                if sample not in self._rep_info:
                    self._rep_info[sample] = {}
                if replicate in self._rep_info[sample]:
                    self.logger.error("Replicate found more than once for \
                                       sample in pool.")
                    raise KeyError("Replicate found more than once for sample \
                                    in pool.")
                temp = tag_pair[0] + "." + tag_pair[1]
                self._rep_info[sample][replicate] = temp

    def process_sort_output_files(self):
        """Process the sort output file to create sample files.

        Read the output file from the sorter class and process it for quality,
        including minimum length, proportion of PCRs and number of pools. Then
        output the reads in a file for each sample, with the tag written on the
        fasta header line.

        Raises
        ------
        IOError
            when the input file does not exist

        """
        # First print out the details :)
        self.__log_in_details()
        # figure out for each pool, what is going on.
        pool_names = self._samp_info.keys()
        self._paired_end = None
        for pool in pool_names:
            sort_file = open(self.input_prefix + "_" + pool + ".tagInfo")
            current_pool = self._samp_info[pool]
            headerline = sort_file.readline()
            headerline = headerline.strip().split()
            ntoks = len(headerline)
            if self._paired_end is None:
                self._paired_end = (ntoks == 6)
            elif self._paired_end != (ntoks == 6):
                self.logger.error("Some sort outputs are paired end, whereas \
                                  others are not.")
                raise IOError("Some sort outputs are paired end, whereas \
                              others are not.")
            type_index = (5 if self._paired_end else 4)
            for line in sort_file:
                toks = line.strip().split()
                # Check if correct combi else continue on
                if toks[type_index] != "C":
                    continue
                tag_pair = (toks[0], toks[1])
                # raise error if tag pair not found in pool - unlikely since
                # this should have been sorted out in sorting.
                if tag_pair not in current_pool:
                    self.logger.error("The tag pair " + toks[0] + " and " +
                                      toks[1] + "not found in pool " + pool)
                    raise KeyError("The tag pair " + toks[0] + " and " +
                                   toks[1] + "not found in pool " + pool + ".")
                (sample, replicate) = current_pool[tag_pair]
                # process it, one way for single end and another for paired end
                if self._paired_end:
                    seq = toks[2:4]
                    count = int(toks[4])
                else:
                    seq = toks[2]
                    count = int(toks[3])
                if sample not in self._haps_info:
                    self._haps_info[sample] = {}
                if seq not in self._haps_info[sample]:
                    self._haps_info[sample][seq] = {}
                    for rep in self._rep_info[sample]:
                        temp_tp = self._rep_info[sample][rep]
                        self._haps_info[sample][seq][rep] = (temp_tp, 0)
                self._haps_info[sample][seq][replicate][1] = count
            sort_file.close()

    def process_haps_info(self):
        """Process the amplicon information and write output files.

        Function to take the processed haps_info dict and filter it to retain
        only the reads that pass these filters and then write out these fasta
        file with only these reads.

        Raises
        ------
        IOError
            when files are missing or cannot be written.

        """
        out_name = self.output_directory + "/" + self.output_prefix + ".fna"
        out_file = open(out_name, "w")
        for samp in self._haps_info:
            sample_haps = self._haps_info[samp]
            for seq in sample_haps:
                counts = []
                tag_pairs = []
                for rep in sample_haps[seq]:
                    (tag_pair, count) = sample_haps[seq][rep]
                    counts.append(count)
                    tag_pairs.append(tag_pair)
                nreps = len(sample_haps[seq])
                valid_counts = sum([1 for x in counts if x >= self.min_count])
                if ((valid_counts*1.0)/nreps) < self.prop_pcr:
                    continue
                fasta_header = ">" + samp + "\t"
                fasta_header += "_".join(tag_pairs) + "\t"
                fasta_header += "_".join([str(x) for x in counts])
                if self._paired_end:
                    if (len(seq[0]) + len(seq[1])) >= self.min_length:
                        out_file.write(fasta_header + "\tread1\n")
                        out_file.write(seq[0] + "\n")
                        out_file.write(fasta_header + "\tread2\n")
                        out_file.write(seq[1] + "\n")
                else:
                    if len(seq) >= self.min_length:
                        out_file.write(fasta_header + "\n")
                        out_file.write(seq + "\n")
        out_file.close()
