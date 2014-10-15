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



with open(sys.argv[1],'r') as filfasta:
	seq=''
	seqs={}
	for l in filfasta:
		if l[:1]==">":
			if seq:
				seqs[name]=seq
			name=l
			seq=''
		else:
			seq+=l[:-1]
	seqs[name]=seq

#pprint(seqs)

with open(sys.argv[2],'r') as filord:
	for i in filord:
		print "%s\n%s" %(i[:-1],seqs[i])

exit(0)

strains={}
ref={}
nb=len(sys.argv)-1
string="Chr\tPos\tRef\t"
lim="\t"
seq=[">ref\n"]
for i in range(nb):
	with open(sys.argv[i+1], 'r') as strain:
		try:
			name = re.search('(.+?)flanksnp\.xls', sys.argv[i+1]).group(1)
		except:
			name = sys.argv[i+1]
		seq.append(">%s\n"%(name))
		for line in strain:
			fields=line.split('\t')
			if fields[0]=='SNP':
				strains[(sys.argv[i+1],fields[1],fields[2])]=fields[4]
				if (fields[1],fields[2]) not in ref.keys():
					ref[(fields[1],fields[2])]=fields[3]
	if i+1==nb:
		lim=""
	string+=name+lim

print string
#pprint(strains)

for k in ref.keys():
#k=ref.keys()[0]
	string=k[0]+"\t"+k[1]+"\t"+ref[k]+"\t"
	seq[0]+=ref[k]
	for s in range(nb):
		try:
			base=strains[(sys.argv[s+1],k[0],k[1])]
		except:
			base=ref[k]
		string+=base
		seq[s+1]+=base
		if s+1==nb:
			string+=""
		else:
			string+="\t"
	print string

#for s in seq:
	#print s
