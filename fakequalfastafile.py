#!/usr/bin/python
# -*- coding: utf-8 -*-

# Villain Adrien
# 2015
# villain.adrien@gmail.com

import os
import re
import sys
import argparse
import warnings
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast import NCBIXML

def generatequals(fasta):
    with open(fasta, 'r') as filin:
        out=''
        seq=''
        for l in filin:
            if l[:1]==">":
                if seq:
                    for i in range(len(seq)):
                        out+="60 "
                    out=out[:-1]
                    out+="\n"
                out+=l
                seq=''
            else:
                seq+=l.strip()
    for i in range(len(seq)):
        out+="60 "
    out=out[:-1]
    return out

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generates a .qual file')
    parser.add_argument("fasta", help='DNA fasta file')
    
    if not sys.argv[1:] :
       sys.argv.append('-h')
    
    prog_name = os.path.basename(sys.argv[0])
    args = parser.parse_args()
    sout=generatequals(args.fasta)
    print sout

