#!/mount/biotools/bin/python2.7
# -*- coding: utf-8 -*-

import sys
from pprint import pprint
import re


description="Takes multiple flanksnps files and outputs a tabulated file with every mutation (lines) for every strain (columns) to a file (last argument) and a fasta file with all the positions for every strain to stdout\nExample : ./combineflanks.py strain1_flanksnps.xls strain2_flanksnps.xls recap.xls > recap.fa"

if not sys.argv[1:] or sys.argv[1] in ["-h","-help","--help"]:
    sys.exit(description)

#parser = argparse.ArgumentParser(description='Extracts flanking regions of SNPS.')
#parser.add_argument("reference", help='Reference file in FASTA format')
#parser.add_argument("snpfile", help='SNP file in VCF format')
#args = parser.parse_args()

strains={}
ref={}
nb=len(sys.argv)-2
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

#print string
#pprint(strains)
with open(sys.argv[-1],'w') as filout:
	filout.write("%s\n" %string)
	for k in sorted(ref.keys(),key=lambda r: int(r[1])):
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
		filout.write("%s\n" %string)

for s in seq:
	print s
