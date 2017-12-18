[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [](#lang-au)

# clermonTyping : Phylogroup typing of E. coli contigs fasta files

## Contents

- [Introduction](#introduction)
- [Installation & Dependencies](#installation)
- [Command line options](#command-line-options)
- [Output Files](#output-files)

## Introduction

Clermont PCR method In-Silico. See abstract.
Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.

## Installation & Dependencies

### Dependencies

* [NCBI BLAST+ blastn](https://www.ncbi.nlm.nih.gov/books/NBK279671/) 
* [Biopython](http://biopython.org/) under [python3](https://www.python.org/downloads/)
* [R](https://cran.r-project.org/) version 3.4.2 or higher

### Instalation

## Command line options

### Main script usage

```
% clermonTyping.sh
Script usage :
	-h					    : print this message and exit
	--fasta					: fasta contigs name(s). Can be separated by an At sign (@) value
	--name					: name for this analysis (optional)
	--treshold			: Option for ClermontTyping, do not use contigs under this size (optional)
```

This script will execute the pipeline blast, mash and python to give the full output (html file). If you need to analyse several fasta files you can list them with a @ sign:
```
% clermonTyping.sh --fasta my_ecoli1.fasta@my_ecoli2.fasta@my_ecoli3.fasta
```

### Clermont Typing without mash and R

If you do not want to use mash analysis and/or R you can independently launch any part of the pipeline.

#### blastn launch

To launch blast you will need to locate the primers.fasta file in the data folder from clermonTyping's installation directory. This contains the essentials primers for PCR amplification.
You will need to format the output in XML format in order to use the clermontyping script.

```
% makeblastdb -in my_fasta.fasta -input_type fasta -out my_fasta -dbtype nucl
% blastn -query ./data/primers.fasta -perc_identity 90 -task blastn -outfmt 5 -db my_fasta -out my_fasta.xml
```
#### Clermontyping launch

The python script will use the output of blastn only in xml format (option -outfmt 5 ).

```
% bin/clermont.py -x my_fasta.xml
```
If you really want to, there are several options for filtering the output.

```
-m/--mismatch <integer> : The maximum number of mismatches in hits. Default = 2.
-l/--length <integer> : The length of the crucial hybridation fragment (seed). Default = 5.
-s/--min_size <integer>: Minimum size for a hit to be counted. This avoid finding primers in smalls contigs.
```

## Output Files

Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.
