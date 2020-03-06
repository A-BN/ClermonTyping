#!/usr/local/bin/python3
# coding: utf-8

############################
# Import
############################
import os
import sys
import argparse
import subprocess
from subprocess import Popen, PIPE
import shutil

############################
# Global Variables
############################

############################
# Functions
############################


def check_file(file):
    """
    Check if a file exist and is readable
    """
    if not os.path.isfile(file) or not os.access(file, os.R_OK):
        print("the file " + file + " is missing or not readable.\n")
        PARSER.print_help()
        sys.exit(1)

def List_to_dic(filename):
    """
    Send a file (three columns : sample, phylogroup and mashgroup with header) to a list
    """

    dico={}
    with open(filename, "r") as f:
        f.readline()
        for line in f:
            try:
                sample, phylo, mash=line.strip().split("\t")
                dico[sample]={}
                dico[sample]['phylo']=phylo
                dico[sample]['mash']=mash
            except ValueError:
                pass
            
    return dico

def Html_to_dic(filename):
    """
    Send a file (three columns : sample, phylogroup and mashgroup with header) to a list
    """

    dico={}
    with open(filename, "r", encoding="utf-8") as f:
        for line in f:
            if '<tr class="odd">' in line or '<tr class="even">' in line:
                sample=f.readline().strip()
                sample_tmp=sample.replace('<td align="left">', '')
                sample=sample_tmp.replace('</td>', '')
                f.readline()
                f.readline()
                phylo=f.readline().strip()
                phylo_tmp=phylo.replace('<td align="left">', '')
                phylo=phylo_tmp.replace('</td>', '')
                mash=f.readline().strip()
                mash_tmp=mash.replace('<td align="left">', '')
                mash=mash_tmp.replace('</td>', '')
                dico[sample]={}
                dico[sample]['phylo']=phylo
                dico[sample]['mash']=mash
            
    return dico

def comparaison(dico_old, dico_new):
    
    with open('check.log', "w") as w:
        
        for sample in dico_old.keys():
            if dico_old[sample]['phylo'] != dico_new[sample]['phylo'] or dico_old[sample]['mash'] != dico_new[sample]['mash']:
                sentence='Error with '+ sample 
                sentence+="\tWaited :" + dico_old[sample]['phylo'] + '/' + dico_old[sample]['mash']
                sentence+="\tNew : " + dico_new[sample]['phylo'] + '/' + dico_new[sample]['mash']
                
                print(sentence)
                w.write(sentence+"\n")
    
    
    

############################
# Main
############################
if __name__ == "__main__":
    
    fastafile='fastafile'
    
    # Clermontyping
    if os.path.isdir('res_check'):
        abso_dir = os.path.abspath('res_check')
        shutil.rmtree(abso_dir)
    cmd = ["../clermonTyping.sh", "--fastafile", fastafile, "--name", "res_check"]
    subprocess.call(cmd)

    ## 
    inputList='ref_resultat'
    check_file(inputList)
    res_waited=List_to_dic(inputList)
    
    inputHtml='res_check/res_check.html'
    check_file(inputHtml)
    res_new=Html_to_dic(inputHtml)
    
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    comparaison(res_waited, res_new)
    if os.stat('check.log').st_size == 0:
        print("Congratulation! Database is updated succesfull.")
    else:
        print("ERROR!!! Look at the check.log file.")
    
    
    
    
    
    
