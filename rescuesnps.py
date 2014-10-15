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
    parser = argparse.ArgumentParser(description='Takes a vcf file and rescues SNPs that have been filtered out if they are known as confident SNPs and pass a less stringent filtering')
    parser.add_argument("known", help='confident sites vcf file')
    parser.add_argument("filtered", help='filtered vcf file')
    parser.add_argument("raw", help='raw vcf file')

    if not sys.argv[1:] :
       sys.argv.append('-h')

    prog_name = os.path.basename(sys.argv[0])
    args = parser.parse_args()

    if not args.known or not args.filtered or not args.raw :
        parser.error("\nVcf file missing.")
    else :
        if not os.access(args.known, os.R_OK):
            parser.error("Unable to read vcf file %s " %args.known)
    	if not os.access(args.filtered, os.R_OK):
            parser.error("Unable to read vcf file %s " %args.filtered)
	if not os.access(args.raw, os.R_OK):
            parser.error("Unable to read vcf file %s " %args.raw)

    knowns=[]
    filts=[]
    toadd=[]
    with open(args.known,'r') as filvcf:
	for line in filvcf:
		if line[:1]!='#':
			fields=line.split("\t")
			knowns.append([fields[0],fields[1]])

    with open(args.filtered,'r') as filtered:
	for line in filtered:
		if line[:1]!='#':
                        fields=line.split("\t")
                        filts.append([fields[0],fields[1]])

    with open(args.raw,'r') as filraw:
	for line in filraw:
		if line[:1]!='#':
                        fields=line.split("\t")
			if [fields[0], fields[1]] in knowns and not [fields[0], fields[1]] in filts:# known snp filtered out in this strain
				#do filtering
				if int(float(fields[5]))>500 and 100*(int(fields[9].split(':')[1].split(',')[1])+1.0)/(int(fields[9].split(':')[1].split(',')[0])+int(fields[9].split(':')[1].split(',')[1]))>60 and int(fields[9].split(':')[1].split(',')[1])>10 :
					toadd.append(line[:-1])

    with open(args.filtered,'a') as filtered:
	for a in toadd:
		filtered.write("%s\n"%a)


