[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [](#lang-au)

# ClermonTyping: an easy-to-use and accurate *in silico* method for *Escherichia* genus strain phylotyping

## Contents

- [Introduction](#introduction)
- [Dependencies & Installation](#dependencies--installation)
- [Command line options](#command-line-options)
- [Output Files](#output-files)
- [Citing](#citing)

## Introduction

The genus Escherichia is composed of Escherichia albertii, E. fergusonii, five cryptic Escherichia clades and E. coli sensu stricto.
Furthermore, the E. coli species can be divided into seven main phylogroups termed A, B1, B2, C, D, E and F. As specific
lifestyles and/or hosts can be attributed to these species/phylogroups, their identification is meaningful for epidemiological
studies. Classical phenotypic tests fail to identify non-sensu stricto E. coli as well as phylogroups. Clermont and colleagues
have developed PCR assays that allow the identification of most of these species/phylogroups, the triplex/quadruplex PCR
for E. coli phylogroup determination being the most popular. With the growing availability of whole genome sequences, we
have developed the ClermonTyping method and its associated web-interface, the ClermonTyper, that allows a given strain
sequence to be assigned to E. albertii, E. fergusonii, Escherichia clades I–V, E. coli sensu stricto as well as to the seven main E.
coli phylogroups. The ClermonTyping is based on the concept of in vitro PCR assays and maintains the principles of ease of
use and speed that prevailed during the development of the in vitro assays. This in silico approach shows 99.4 % concordance
with the in vitro PCR assays and 98.8 % with the Mash genome-clustering tool. The very few discrepancies result from
various errors occurring mainly from horizontal gene transfers or SNPs in the primers. We propose the ClermonTyper as a
freely available resource to the scientific community at:

__http://clermontyping.iame-research.center/.__

## Dependencies Installation

### Dependencies

* [NCBI BLAST+ blastn](https://www.ncbi.nlm.nih.gov/books/NBK279671/) 
* [Biopython](http://biopython.org/) under [python3](https://www.python.org/downloads/)
* [pandoc](https://pandoc.org/) 
* [R](https://cran.r-project.org/) version 3.4.2 or higher with the following packages installed:
	* readr
	* dplyr
	* tidyr
	* stringr
	* knitr

### Installation
...
## Command line options

### Main script usage

```
% clermonTyping.sh
Script usage :
	-h				: print this message and exit
	--fasta				: fasta contigs file(s). If multiple files, they must be separated by an arobase (@) value
	--name				: name for this analysis (optional)
	--threshold			: Option for ClermontTyping, do not use contigs under this size (optional)
```

This script will execute the pipeline blast, mash and python to give the full output (html file). If you need to analyse several fasta files you can list them with a @ sign (absolute path required):
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

The default analysis name is analysis_date and every results are stored in the corresponding folder.

* analysis.html : final output with the main script pipeline. Gives informations about phylogroups with 2 differents methods (mash and clermontyping).
* analysis_phylogroups.txt : final output of clermontyping
* analysis.R : intermediate file for producing the html output. You can run this Rscript alone.
* strain.xml : intermediate file. Goes with the "db" folder. Output of blastn.
* strain_mash_screen.tab : intermediate file. Output of mash.

### HTML output

This is a table with each line is a fasta file you analyzed.
* The file name
* The quadruplex: you will find a representation of presence(+)/absence(-) for the 4 genes described in [Clermont O, Christenson JK, Denamur E, Gordon DM. The Clermont Escherichia coli phylo-typing method revisited: improvement of specificity and detection of new phylo-groups](https://www.ncbi.nlm.nih.gov/pubmed/23757131). Respectively, arpA, chuA, yjaA and TspE4.C2 .
* The supp column: When the quadruplex PCR give an ambiguity, this column show the alleles specific for groups C and E.
* The phylogroup detected on your sample
* The mash prediction. This prediction is only informative and is not used for phylogroup assignation.

### Phylogroup output

analysis_phylogroups.txt is a TSV file containing every fasta file analyzed with blastn + clermont method.
Exemple:

```
ROAR344_fergusonii.fasta	['trpA', 'trpBA', 'aesI']	['-', '-', '-', '-']	[]	Fergusonii
```
* 1st column : Name of the sample analyzed.
* 2nd column : All genes detected with PCR (including verification genes that may be not usefull). 
* 3rd column : In the quadruplex column you will find a representation of presence(+)/absence(-) for the 4 genes described in [Clermont O, Christenson JK, Denamur E, Gordon DM. The Clermont Escherichia coli phylo-typing method revisited: improvement of specificity and detection of new phylo-groups](https://www.ncbi.nlm.nih.gov/pubmed/23757131). Respectively, arpA, chuA, yjaA and TspE4.C2 .
* 4th column : When the quadruplex PCR give an ambiguity, this column show the alleles specific for groups C and E.
* 5th column : final conclusion about phylogroup.

## Citing
Please cite:

[Beghain, J., Bridier-Nahmias, A., Le Nagard, H., Denamur, E. & Clermont, O. ClermonTyping: an easy-to-use and accurate in silico method for Escherichia genus strain phylotyping. Microbial Genomics (2018). doi:10.1099/mgen.0.000192](http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000192)

