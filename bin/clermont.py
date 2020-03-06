#!/usr/bin/env python3
"""
Current Version - 1.4.0: Fév. 2019
        1.4.0:
            - Adding an essential primer for G group
        1.3.0:
            - Adding an essential primer for fergusonii amplification
            - Changine behavior for other fergusonii to Unkown
        1.2.0:
            - Changing conclusion about trpA missing
            - Profile --+- is now fergusonii
        1.1.0:
			- Adding disambiguation in profile +-+- with clade I
		1.0.0:
			- Final version for github
			- Removing option -c
		0.5.0:
			- Update outputs
		0.4.0:
			- Update decision making for phylogroups and clades
		0.3.0:
			- option min-size for filtering small contigs
		0.2.0:
			- Clade markers method added
			- Control do discrimine C and E groups
			
Clermont method from a blast results. This program check the presence of primers on a fasta sequence (whole genome or contigs).
This is an 'in silico' PCR that sould give the same results as a regular PCR in lab (i.e. amplifications of fragments)

Required: 	- Python 3.5.3 or later
			- BioPython [1]
			- (For input dataset files) BLAST 2.6.0+ [2]
			
Clermont O, Christenson JK, Denamur E, Gordon DM. The Clermont Escherichia coli phylo-typing method revisited: improvement of specificity and detection of new phylo-groups.
Environ Microbiol Rep. 2013 Feb;5(1):58-65. doi: 10.1111/1758-2229.12019. Epub 2012 Dec 24.

Clermont O, Gordon DM, Brisse S, Walk ST, Denamur E. Characterization of the cryptic Escherichia lineages: rapid identification and prevalence.
Environ Microbiol. 2011 Sep;13(9):2468-77. doi: 10.1111/j.1462-2920.2011.02519.x. Epub 2011 Jun 8.

[1] http://biopython.org
[2] https://www.ncbi.nlm.nih.gov/books/NBK279690/
"""
import argparse
import ast
import logging
import math
import multiprocessing as mp
import subprocess
import os
import builtins
from Bio.Blast import NCBIXML

###############################################################################
#############                CLASS                            #################
###############################################################################

#Definition of Primers names and size of PCR product (bp)
class Primers:
    def __init__(self):
        self.names = {"trpA": 783, "trpBA": 489, "chuAalbertii": 136, "citPferg": 300, "chuA": 288, "yjaA": 211, "TspE4.C2": 152, "arpA": 400, "ArpAgpE": 301, "trpAgpC": 219, "aesI": 315, "aesII": 125, "chuIII": 183, "chuIV": 461, "chuV": 600, "ybgD": 177}

###############################################################################
#############                FUNCTIONS                            #############
###############################################################################

def ArgumentsParser():
  parser=argparse.ArgumentParser(description = '''In sillico Clermont method.''',
                                 epilog = '''Need a blastn result in xml format.''')
  parser.add_argument('-x', '--xml', type = str, required = True,
                      help = 'The blastn results in xml format only.')
  parser.add_argument('-m', '--mismatch', type = int, required = False,
                      help = 'The maximum number of mismatches in hits. Default = 2')
  parser.add_argument('-l', '--length', type = int, required = False,
                      help = 'The length of the crucial hybridation fragment (seed). Default = 5.')
  parser.add_argument('-s', '--min_size', type = int, required = False,
                      help = 'Minimum size for a hit to be counted. This avoid finding primers in smalls contigs.')
  logger_args = parser.add_argument_group('Logger configuration')
  logger_args.add_argument( '--log_level',
                         type = str,
                         nargs = '?',
                         default = 'INFO',
                         help = 'log level',
                         choices = ['ERROR', 'error', 'WARNING', 'warning', 'INFO', 'info', 'DEBUG', 'debug'],
                         required = False )
  logger_args.add_argument( '--log_file',
                         type = str,
                         nargs = '?',
                         help = 'log file (use the stderr by default)',
                         required = False )
  return(parser.parse_args())

def blast_parser(xml_content, min_size):
    blast_records = NCBIXML.parse(xml_content, 0)
    blast_records = list(blast_records)
    quadruplex = dict()
    for blast_record in blast_records:
        logger.debug(blast_record.query)
        #Check the name
        all_primers = Primers()
        query_name = blast_record.query.split("_")
        if query_name[0] not in all_primers.names:
            logger.error("Error in primer name: "+blast_record.query+" . Sould be one of "+str(all_primers.names.keys()))
            break
        if query_name[1] not in ['F', 'R']:
            logger.error("Error in primer name: "+blast_record.query+" . Sould be _F or _R for forward and reverse.")
            break
        query_length = blast_record.query_length
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                #Check number of mismatches
                num_mismatches = query_length - hsp.identities
                if num_mismatches <= max_mismatches:
                    #Check if last nucleotides are presents
                    #i.e the end of the query is at least equals to the query_length minus cutoff
                    query_length_limit = query_length-cutoff
                    if hsp.query_end < query_length_limit:
                        logger.debug("Mismatch found at the begining of seed: query end at "+str(hsp.query_end))
                        continue
                    #Check if any mismatches occurs in 3' (/!\ Primers have to be 5' 3' oriented - strand 1)
                    if hsp.match[len(hsp.match)-seed_length:].count('|') >= seed_length-cutoff:
                        #Check size of hit
                        if min_size > 0 and alignment.length < min_size:
                            logger.debug("Matching sequence too small to be counted.")
                            continue
                        #Record. If it already exists, separate results by the '$' symbol.
                        #This mean that there is more than one good hit on sbjct.
                        #Also store the frame of the sbjct (this is the second one in the tuple hsp frame)
                        logger.debug("Found one primer for "+blast_record.query+" on "+alignment.title)
                        if query_name[0] in quadruplex:
                            if query_name[1] in quadruplex[query_name[0]]:
                                quadruplex[query_name[0]][query_name[1]] = quadruplex[query_name[0]][query_name[1]]+'$'+str(hsp.sbjct_start)+'&'+alignment.title+'&'+str(hsp.frame[1])
                            else:
                                quadruplex[query_name[0]][query_name[1]] = str(hsp.sbjct_start)+'&'+alignment.title+'&'+str(hsp.frame[1])
                        else:
                            quadruplex[query_name[0]] = {}
                            quadruplex[query_name[0]][query_name[1]] = str(hsp.sbjct_start)+'&'+alignment.title+'&'+str(hsp.frame[1])
                    else:
                        logger.debug("One mismatch found in seed for primer "+blast_record.query+". Dicard it.")
                else:
                    logger.debug("Too much mismatches for primer "+blast_record.query+" : "+str(num_mismatches))
    return(quadruplex)

#Check if the primers are ok for Clermont PCR
#Check the frames and size of fragments
#Return all PCR products names
def pcr_parser_groups(quadruplex):
    all_primers = Primers()
    primers = all_primers.names
    #Valid pcr_products to be returned
    pcr_products = []
    #Iterate over all primers (chuA, yjaA, etc.)
    for primer in quadruplex:
        #For each strand (_F or _R i.e. forward or reverse)
        positions = {"F": [], "R": []}
        frames = {"F": [], "R": []}
        names = {"F": [], "R": []}
        for strand in quadruplex[primer]:
            #Get all hits
            hits = quadruplex[primer][strand].split("$")
            for hit in hits:
                hit_informations = hit.split('&')
                #Get all hits position
                positions[strand].append(hit_informations[0])
                #Get all hits name
                names[strand].append(hit_informations[1])
                #Get all hits frame
                frames[strand].append(hit_informations[2])
        #Check the sizes
        f = 0
        for position_f in positions['F']:
            r = 0
            for position_r in positions['R']:
                #Check if positions are on the same fragment
                if names['F'][f] != names['R'][r]:
                    logger.debug("hit names differents for "+primer)
                    r=r+1
                    continue
                insert_size = math.fabs(int(position_f)-int(position_r))
                #Size of insert must be near the given insert size +/- 20% bp
                minimum_size = primers[primer]-(primers[primer]*10/100)
                maximum_size = primers[primer]+(primers[primer]*10/100)
                if (insert_size >= minimum_size) and (insert_size <= maximum_size):
                    #Insert size is correct but do the frame is correct?
                    #frame must be one at -1 and one at 1
                    frame_shift = int(frames['F'][f])+int(frames['R'][r])
                    if frame_shift == 0:
                        pcr_products.append(primer)
                        logger.debug(primer+" found on "+names['F'][f]+" .At positions: "+str(position_f)+" - "+str(position_r)+" .Insert size is "+str(insert_size)+"bp.")
                    else:
                        logger.debug("Frames of primers in wrong way for amplification: "+frames['F'][f]+" and "+frames['R'][r])
                else:
                    logger.debug("insert size for "+primer+" out of allowed length ("+str(primers[primer])+"): found "+str(insert_size)+"bp.")
                r=r+1
            f=f+1
    #return the quadruplex markers in the presence/absence matrix
    #[arpA, chuA, yjaA, TspE4.C2]
    quadruplex = []
    quadruplex_names = ['arpA', 'chuA', 'yjaA', 'TspE4.C2']
    for name in quadruplex_names:
        if name in pcr_products:
            quadruplex.append('+')
        else:
            quadruplex.append('-')
    #return the specific primers (grp E or C disambiguity)
    specific = []
    specific_names = ["chuAalbertii", "citPferg", "aesI", "aesII", "chuIII", "chuIV", "chuV", 'ArpAgpE', 'trpAgpC']
    for name in specific_names:
        if name in pcr_products:
            specific.append(name)
    return(pcr_products, quadruplex, specific)

#Find the phylogroup from the markers (clermont method)
#Take a list of presents markers
def find_phylo_group(markers):
    all_primers = Primers()
    primers = all_primers.names
    #Internal control
    if "chuAalbertii" in markers:
        return("albertii")
    if "citPferg" in markers:
        return("fergusonii")
    if "arpA" in markers:
        if "chuA" in markers:
            if "yjaA" in markers:
                if "TspE4.C2" in markers:
                    return("Unknown")
                else:
                    if "aesI" in markers:
                        return("cladeI")
                    else:
                        #Double check for clade E
                        if "ArpAgpE" in markers:
                            return("E")
                        else:
                            return("E or cladeI")
            else:
                if "TspE4.C2" in markers:
                    if "ArpAgpE" in markers:
                        return("E")
                    else:
                        return("D")
                else:
                    if "ArpAgpE" in markers:
                        return("E")
                    else:
                        return("D")
        else:
            if "yjaA" in markers:
                if "TspE4.C2" in markers:
                    return("Unknown")
                else:
                    if "aesI" in markers:
                        return("cladeI")
                    else:
                        if "trpAgpC" in markers:
                            return("C")
                        else:
                            return("A")
            else:
                if "TspE4.C2" in markers:
                    return("B1")
                else:
                    return("A")
    else:
        if "chuA" in markers:
            if "yjaA" in markers:
                if "TspE4.C2" in markers:
                    return("B2")
                else:
                    return("B2")
            else:
                if "TspE4.C2" in markers:
                    if "ybgD" in markers:
                        return("G")
                    else:
                        return("B2")
                else:
                    if "ybgD" in markers:
                        return("G")
                    else:
                        return("F")
        else:
            if "yjaA" in markers:
                if "TspE4.C2" in markers:
                    return("B2")
                else:
                    if "aesI" in markers:
                        return("cladeI")
                    elif "aesII" in markers:
                        return("cladeII")
                    else:
                        return("Unknown")
            else:
                if "TspE4.C2" in markers:
                    return("B1")
                else:
                    clades = []
                    if "aesII" in markers:
                        clades.append("cladeII")
                    if "chuIII" in markers:
                        clades.append("cladeIII")
                    if "chuIV" in markers:
                        clades.append("cladeIV")
                    if "chuV" in markers:
                        clades.append("cladeV")
                    if clades:
                        return(' '.join(clades))
                    else:
                        if "trpA" not in markers and "trpBA" not in markers:
                            return("Non Escherichia")
                        return("Unknown")

###############################################################################
#############                MAIN                                 #############
###############################################################################
if __name__ == "__main__":
    args=ArgumentsParser()
    # Logger config
    logging_std_format = '[%(levelname)s] %(message)s'
    logging_debug_format = '%(asctime)s [%(levelname)s] [%(threadName)s - %(name)s] %(message)s'
    log_level = args.log_level.upper()
    if ( log_level == 'DEBUG' ):
        logging_std_format = logging_debug_format
    logging_datefmt = '%d/%m/%Y - %H:%M:%S'
    if ( args.log_file != None ):
        logging.basicConfig( format = logging_std_format,
                             datefmt = logging_datefmt,
                             filename = args.log_file,
                             filemode = 'w',
                             level = log_level )
    else:
        logging.basicConfig( format = logging_std_format,
                             datefmt = logging_datefmt,
                             level = log_level )
    logger = logging.getLogger( 'main' )
    
    #Default number of mismatches
    max_mismatches = 6
    if args.mismatch:
        max_mismatches = args.mismatch
    #Default length of seed (the last N nucleotides in the primer are essentials)
    seed_length = 3
    if args.length:
        seed_length = args.length
    #Default number of mismatches allowed in the seed.
    cutoff = 0
    #Default minimum contig size = 0 (disabled)
    min_contig_size = 0
    if args.min_size:
        min_contig_size = args.min_size

    #Check xml file
    try:
        result_handle = open(args.xml, 'r')
        #Launch Blast parsing
        results = blast_parser(result_handle, min_contig_size)
        #print(results)
    except EnvironmentError:
        logger.info("Fail to read xml file "+args.xml)
        exit(1)
    pcr_products, quadruplex, specific = pcr_parser_groups(results)
    print(pcr_products, end="\t")
    print(quadruplex, end="\t")
    print(specific, end="\t")
    phylo_group = find_phylo_group(pcr_products)
    print(phylo_group)
