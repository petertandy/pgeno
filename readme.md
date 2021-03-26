![pgeno](logo.svg)
## pgeno
Pgeno is a multithreaded genotyping program based on k-mer detection in FASTQ data. Given a small, specially annotated fasta file for reference, Pgeno will compile RegEx phrases to search for k-mers of a user-defined length within reads of the provided fastq-formatted dataset.

The resulting output will be the determined genotype of each sample based on the number of detected sequences between the k-mers. This output can then be interpreted to determine the presence or absence of a specific allele or haplotype.

#### What's a k-mer?
A `k-mer` is a sequence of DNA of length `k`. For example, `GATC` is a "4-mer" because it's a sequence of length 4.

Any given sequence of DNA can be arbitrarily broken into several combinations of sequences of varying lengths.

`GATC`'s' k-mers:

k | k-mers
- | ------
1 | G, A, T, C
2 | GA, AT, TC
3 | GAT, ATC
4 | GATC

As demonstrated in the table above, longer k-mers in a given sequence are less abundant than shorter k-mers.

*Generally speaking*, a k-mer of sufficient length is enough to identify a certain location in a genome. A PCR primer is a good example of this.


### Requirements

Pgeno is written with `Python 3.8` in mind. A more recent version of Python should work, but older versions likely will not.

Pgeno has two third-party dependencies in addition to the Python standard library:

* Both third-party requirements are available with `pip`, and obtaining them is covered in the Setup section below.
    1. [BioPython](https://github.com/biopython/biopython)
    1. [regex](https://bitbucket.org/mrabarnett/mrab-regex/src/hg/)
        * *Not to be confused with the standard library's `re`*

Pgeno has been tested and developed on Windows 10 and Ubuntu 20. It *ought to* run on any platform capable of installing Python 3.8 and the dependencies, but is untested on other platforms.

### Setup

* Clone or download this repository.
* Open a terminal and install the requirements with `pip install -r /path/to/pgeno/requirements.txt` (Please replace "/path/to/pgeno/" with the path where you have cloned or downloaded this repo).

* You will need two things to use pgeno:
    1. A FASTQ format dataset.
    1. A specially annotated reference file.

Your FASTQ data should be raw single-read DNA data.

The specially annotated reference file you can make from either a genomic reference you already have, or a known sequence you're looking for.
It must be formatted as such:

Filename: `locus.fa`

    >my_locus
    GGAATTCCGGAATTCC[GATC]GGAATTCCGGAATTCC
    

*Note! Other characters like hyphens and asterisks are unsupported and will lead to undefined behavior.*

This is a modified FASTA format. The square braces are the most important part. The square braces denote the area you would like pgeno to genotype.

Pgeno will look for these braces in this modified reference file (or "locus file") and will use their position to determine k-mers from the surrounding sequence.

The sequence within the square braces is ignored.

Once you have a FASTQ dataset and a modified reference file, you can go on to using pgeno in the Usage section below.

### Usage

    pgeno.py [-h] -d ./path/to/DATA -l ./path/to/locus.fa [-o ./path/to/OUTPUT] [-k N] [-m N] [-g N] [-t N] [-z N] [--serial]

    required arguments:
      -d ./path/to/DATA, --data ./path/to/DATA
                            path to the folder or file which holds the samples to genotype
      -l ./path/to/locus.fa, --locus ./path/to/locus.fa
                            path to the locus file.

    optional arguments:
      -h, --help            show this help message and exit
      -o ./path/to/OUTPUT, --out ./path/to/OUTPUT
                            path to the folder where the processed information should be written
      -k N, --kmerlength N  The length of k-mer to be created based on the supplied FASTA files.
      -m N, --minreads N    Minimum number of reads that must be discovered to consider genotyping.
      -g N, --threshold N   Fractional representation a sensed haplotype must reach to be considered for genotyping.
      -t N, --threads N     Number of threads to use. 0 to let pgeno decide.
      -z N, --fuzzy N       Number of fuzzy errors (substitutions) allowed when matching k-mers. A higher number is more permissive, but computationally
                            slower.
      --serial              Specify this flag to calculate genotypes in serial instead of parallel (will be very slow).

#### Example (Windows 10)

* Pgeno is located at `C:\pgeno`
* Your dataset is located at `C:\Documents\DATASET`
* Your locus file is located at `C:\Documents\my_locus.fa`

You can open a terminal and run `python3 C:\pgeno\pgeno.py -d C:\Documents\DATASET -l C:\Documents\my_locus.fa -o C:\Documents\my_genotypes -k 10 -m 40 -g 0.15 -t 8 -z 1`

This command will:

* Run the pgeno script with Python 3
* With the `-d` flag, use the data in the folder `C:\Documents\DATASET`
* With the `-l` flag, use the locus file `C:\Documents\my_locus.fa`
* With the `-o` flag, place the resulting data in the `C:\Documents\my_genotypes` folder
    * The folder will be created if it doesn't already exist
* With the `-k` flag, 10-mers will be extracted from the locus file and used in processing
* With the `-m` flag, a bare minimum of 40 reads will have to be detected in order to genotype a given sample
* With the `-g` flag, a "genotyping threshold" of 0.15 (15%) is set. Only haplotypes which make up at least 0.15 of the total detected will be considered a part of the true genotype for a given sample
* With the `-t` flag, pgeno will use 8 threads to process the data in parallel
* With the `-z` flag, one "fuzzy match" substitution is allowed, meaning that k-mers which are one letter different will be allowed to match in addition to the k-mers which match perfectly.
    * These are *only* substitutions (i.e. `GATC`->`CATC`), *not* insertions and deletions (`GATC`->`G-TC` or `GATC`->`GAATC`)

### Output
In your output folder, pgeno will create a `genotypes.csv` file. This file contains the final genotypes determined per sample as well as the counts of each of the displayed alleles found in the sample.

The alleles displayed are the most abundant haplotypes that were determined during processing the FASTQ data.

*For example:*

Filename: `genotypes.csv`

Dataset: `sample1.fastq`, `sample2.fastq`, `sample3.fastq`

sample | my_locus_genotypes | my_locus_distributions
------ | ------------------ | ----------------------
sample1 | GATC/- | 12000/7000
sample2 | GATC/GATC | 9750
sample3 | None |

This table shows:

1. The "sample1" sample was genotyped as heterozygous for `GATC` and `-`. 12000 instances of `GATC` were found, and 7000 `-`s were found.
    * A hyphen means that flanking k-mers were positively identified, but between them was nothing.
1. "sample2" was genotyped as homozygous `GATC`. 9750 instances of `GATC` were found.
1. "sample3" was not genotyped because there were not enough reads which positively contained the k-mers.


In addition to the genotypes file, there will be an `intermediate_data` folder, which contains one csv file for every sample. These files are very messy and are primarily included for debugging and sanity checking.

The header fields in these files contain *all* of the detected haplotypes across all of the samples, even if they were only found once in only one sample.

The next line is the count for each of those haplotypes found in the sample, an empty value representing 0.
