# Damian: DAMe's evil twin

This program is a updated and impoved version of [DAMe](https://github.com/shyamsg/dame) and its predecessor [DAMe](https://github.com/lisandracady/DAMe). Like DAMe, Damian is built to metabarcoding sequence data and make it ready for processing using downstream processing tools such as [obitools](https://metabarcoding.org/obitools), and the [DADA2](https://benjjneb.github.io/dada2/) pipeline.

## Usage
Currently, __Damian__ is a wrapper script which encompasses two core functions, viz., _sort_ and _filter_. 

### sort
The _sort_ function demultiplexes the pooled sequence data into barcodes corresponding to each of the tag combinations that were used in the pool. The usage of sort is given below.

```
usage: Damian sort [-h] [-p PrimerFile] [-p1 FwdPrimer] [-p2 RevPrimer] -t
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
The primer file contains just one line, with two dna sequences on it, seperated by whitespace. The first sequence is the forward primer, and the second sequence is the reverse primer. Instead of specifying the primers in a file, one can specify the primers on the command line directly, using teh `-p1` and `-p2` arguments. Primer sequences need to be specified for __damian__ to function properly. Primers sequences are allowed to have IUPAC ambiguity in them. 

#### Tag file
The tag file contains a list of all the tag sequences used in the pooling of multiple samples in one pool. The format of the tag file is a two column file, with each line representing a tag sequence. The first column is the tag name whereas the second column is the tag sequence. Note that, unlike the primer sequence, the tag sequences are ___NOT__ allowed to have IUPAC ambiguity codes in them. 

#### Pool information file
The pool information file contains information on the fastq sequence files for each of the pools in a study. The file contains 2 or 3 columns, with each line containing information on a single pool. The first column is the name of the pool (case sensitive), the second column is the name of the single end fastq file (Read 1 if paired end) and the third column contains the name of the second read in a paired end sequencing experiment (absent if single end). Note that all the pools must have the same type of sequencing, i.e., either all the pools must use single end sequencing (or collapsed using AdapterRemoval/cutadapt), or paired end sequencing. 


