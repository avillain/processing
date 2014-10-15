#!/mount/biotools/bin/python2.7
# -*- coding: utf-8 -*-

# Villain Adrien
# 2014     
# avillain@pasteur.fr

import os
import sys
import argparse
import re

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Takes a vcf file and extracts position and length of deletion/insertion')
    parser.add_argument("vcf", help='vcf file')

    if not sys.argv[1:] :
       sys.argv.append('-h')

    prog_name = os.path.basename(sys.argv[0])
    args = parser.parse_args()

    if not args.vcf :
        parser.error("\nVcf file missing.")
    else :
        if not os.access(args.vcf, os.R_OK):
            parser.error("Unable to read vcf file %s " %args.vcf)
    
    with open(args.vcf,'r') as filvcf:
	for line in filvcf:
		if line[:1]!='#':
			fields=line.split("\t")
			delet=len(fields[3])
			inser=len(fields[4])
			if delet>inser:
			    leng=-(delet-inser)
			else:
			    leng=inser-delet
			print "%s\t%s\t%d"  %(fields[0],fields[1],leng)

