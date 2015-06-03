#!/usr/bin/python
# -*- coding: utf-8 -*-

# Villain Adrien
# 2015
# villain.adrien@gmail.com

import os
import sys
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parses a .tblout hmm3r file')
    parser.add_argument("tblout", help='tabular hmm3r output file')
    parser.add_argument("-e","--evalue", help='maximum e-value for hit',
                        default=1e-5)
    
    if not sys.argv[1:] :
       sys.argv.append('-h')
    
    prog_name = os.path.basename(sys.argv[0])
    args = parser.parse_args()
    
    if not args.tblout :
        parser.error("\ntblout file missing.")
    else :
        if not os.access(args.tblout, os.R_OK):
            parser.error("Unable to read file %s " %args.tblout)
    
    with open(args.tblout,"r") as filtblout:
        for lin in filtblout:
            if lin[:1]!='#':
                fields=lin.split()
                evalglob=float(fields[4])
                evaldom=float(fields[7])
                if evalglob < float(args.evalue) and evaldom < float(args.evalue):
                    print "%s\t%s\t%s" %(fields[0][:-2],fields[4],fields[7])


