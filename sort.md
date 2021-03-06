### sort
The _sort_ function demultiplexes the pooled sequence data into barcodes corresponding to each of the tag combinations that were used in the pool. The usage of sort is given below.

```
usage: Begum sort [-h] [-p PrimerFile] [-p1 FwdPrimer] [-p2 RevPrimer] -t
                   TagFile -s SampleInformationFile -l PoolInformationFile
                   [-m] [-pm PrimerMismatches] [-tm TagMismatches]
                   [-mo MinOverlap] [-mm OverlapErrRate] [-d OutDirectory]
                   [-o OutPrefix]

optional arguments:
  -h, --help            show this help message and exit
  -p PrimerFile, --primers PrimerFile
                        File with forward and reverse primer sequence (Format:
                        ForwardPrimer ReversePrimer)
  -p1 FwdPrimer, --fwdPrimer FwdPrimer
                        Sequence of forward primer
  -p2 RevPrimer, --revPrimer RevPrimer
                        Sequence of reverse primer
  -t TagFile, --tags TagFile
                        File with tag name and sequence (Format: TagName
                        FwdTagSequence)
  -s SampleInformationFile, --sampleInfo SampleInformationFile
                        File with tag combo and pool for each sample (Format:
                        Sample FwdTagName RevTagName PoolName)
  -l PoolInformationFile, --pool PoolInformationFile
                        File with pool information (Format: Poolname
                        Read1Fastq [Read2Fastq])
  -m, --allowMultiplePrimers
                        Allow more one occurrance of the primer sequence in
                        read. (Default False)
  -pm PrimerMismatches, --primerMismatches PrimerMismatches
                        Number of mismatches in primer. (Default 0)
  -tm TagMismatches, --tagMismatches TagMismatches
                        Number of allowed mismatches in tags. (Default 0)
  -mo MinOverlap, --merge_overlap MinOverlap
                        Merge read1 and read2 if overlapping by given number
                        of bases or more (>=5) __NOT IMPLEMENTED YET__
  -mm OverlapErrRate, --merge_errors OverlapErrRate
                        Rate of mismatches allowed in overlap between reads:
                        range [0,0.2] __NOT IMPLEMENTED YET__
  -d OutDirectory, --output_directory OutDirectory
                        Output directory. (Default: .)
  -o OutPrefix, --output_prefix OutPrefix
                        Prefix for output files. (Default : '')
```
The _sort_ function requires 4 files, which are now described in detail.
#### Primer file
The primer file contains just one line, with two dna sequences on it, seperated by whitespace. The first sequence is the forward primer, and the second sequence is the reverse primer. Instead of specifying the primers in a file, one can specify the primers on the command line directly, using teh `-p1` and `-p2` arguments. Primer sequences need to be specified for __begum__ to function properly. Primers sequences are allowed to have IUPAC ambiguity in them.

#### Tag file
The tag file contains a list of all the tag sequences used in the pooling of multiple samples in one pool. The format of the tag file is a two column file, with each line representing a tag sequence. The first column is the tag name whereas the second column is the tag sequence. Note that, unlike the primer sequence, the tag sequences are ___NOT___ allowed to have IUPAC ambiguity codes in them.

#### Pool information file
The pool information file contains information on the fastq sequence files for each of the pools in a study. The file contains 2 or 3 columns, with each line containing information on a single pool. The first column is the name of the pool (case sensitive), the second column is the name of the single end fastq file (Read 1 if paired end) and the third column contains the name of the second read in a paired end sequencing experiment (absent if single end). Note that all the pools must have the same type of sequencing, i.e., either all the pools must use single end sequencing (or collapsed using AdapterRemoval/cutadapt), or paired end sequencing.

#### Sample information file
The sample information file contains information on the various replicates for each sample. This is a 4 column file, with the first column being the sample name, the second column being the name of the forward tag, the third column being the name of the reverse tag, and the last column being the pool that this sample-tag combination is present in. It is important to note that the order of the samples is ___NOT___ important, and that the replicate number for a sample is assigned automatically. Further, the name of the samples, pools and tags are case sensitive and must be exactly the same across differnt files.
