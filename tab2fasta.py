#!/usr/bin/python
# -*- coding: utf-8 -*-

# Villain Adrien
# 2015
# villain.adrien@gmail.com

import os
import sys
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Converts a .tsv file to .fa')
    parser.add_argument("tsv", help='tabular sequence file')
    
    if not sys.argv[1:] :
       sys.argv.append('-h')
    
    prog_name = os.path.basename(sys.argv[0])
    args = parser.parse_args()
    
    if not args.tsv :
        parser.error("\ntsv file missing.")
    else :
        if not os.access(args.tsv, os.R_OK):
            parser.error("Unable to read tsv file %s " %args.tsv)
    
    with open(args.tsv,"r") as filtab:
        flag=0
        for lin in filtab:
            if flag:
                if len(lin.split())==2:
                    print ">%s\n%s" %(lin.split()[0],lin.split()[1])
                else:
                     print 2> "SPLIT INCORRECT line %s" %lin
            elif lin[:8]=="ContigID":
                flag=1

