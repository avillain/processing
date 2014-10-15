#!/mount/biotools/bin/python2.7
# -*- coding: utf-8 -*-

import sys
import argparse
from pprint import pprint
import re

#parser = argparse.ArgumentParser(description='Extracts flanking regions of SNPS.')
#parser.add_argument("reference", help='Reference file in FASTA format')
#parser.add_argument("snpfile", help='SNP file in VCF format')
#args = parser.parse_args()


combine=sys.argv[1]
snpfile=sys.argv[2]
pileupfile=sys.argv[3]

covdict={}
with open(pileupfile, 'r') as refpil:
        for line in refpil:
		fields=line.split("\t")
        	covdict[(fields[0],fields[1])]=int(fields[3])

#pprint(covdict)
#exit(0)
#print covdict.keys()[0]

posdict={}
with open(snpfile, 'r') as snp:
	#print "Type\tContig\tPos\tRef\tVar\tFlanking seq in ref"
        for line in snp:
		if line[:1]!="#":
        		fields=line.split("\t")
			contig=fields[0]
                	pos=fields[1]
                	var=fields[4]
                	posdict[(contig,pos)]=var

#print posdict.keys()[0]
#pprint(posdict)
#pprint(strains)
#exit(0)

def getallele(posdict,covdict,ref,tupl):
	try:
		cov=covdict[tupl]
	except:
		cov=0
	if tupl in posdict.keys():
		#print posdict[tupl]
		return posdict[tupl]
		#si la position est trouvée dans les snps de la souche
	else:
		if cov>50:
			return ref
			#renvoie l'allèle de référence si (contig,position) couvert ET absent dans la souche
		else:
			return "N"
			#si position non couverte, supposer qu'elle est variante

#print "Chr\tPos\tRef\t%s" %snpfile
ou=">%s\n"%snpfile
with open(combine, 'r') as comb:
	junk=comb.readline()
	for line in comb:
		fields=line.split('\t')
		contig=fields[0]
		pos=fields[1]
		ref=fields[2]
		ou+="%s"%getallele(posdict,covdict,ref,(contig,pos))
		#print "%s\t%s\t%s\t%s" %(contig,pos,ref,getallele(posdict,covdict,ref,(contig,pos)))

print ou
