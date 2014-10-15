#!/mount/biotools/bin/python2.7
# -*- coding: utf-8 -*-

import sys
import argparse
from pprint import pprint

parser = argparse.ArgumentParser(description='Extracts flanking regions of SNPS.')
parser.add_argument("reference", help='Reference file in FASTA format')
parser.add_argument("snpfile", help='SNP file in VCF format')
args = parser.parse_args()

fastadict={}

with open(args.reference, 'r') as refasta:
	seq=""
	name=""
	for line in refasta:
		if line[:1]=='>':
			if seq!="":
				fastadict[name]=seq
			name=line[1:-1]
			seq=""
		else:
			seq+=line[:-1]
	fastadict[name]=seq.upper()

#pprint(fastadict)

with open(args.snpfile, 'r') as snp:
	print "Type\tContig\tPos\tRef\tVar\tFlanking seq in ref"
	for line in snp:
		if line[:1]!="#":
			fields=line.split("\t")
			id=fields[0]
			pos=int(fields[1])
			#print "SNP : %s %d %s seq in ref : %s" %(id, pos, fields[3], fastadict[id][pos-1:pos])
			##### we have to be sure that the contig is not too short to extract
			if(len(fastadict[id])<=80) : # then just output the whole thing
				extr=fastadict[id][:pos-1].lower()
				extr+=fastadict[id][pos-1:pos].upper()
				extr+=fastadict[id][pos:].lower()
			elif(len(fastadict[id][:pos-1])<40): # print more seq after the snp
				extr=fastadict[id][:pos-1].lower()
				extr+=fastadict[id][pos-1:pos].upper()
				extr+=fastadict[id][pos:pos+(80-len(fastadict[id][:pos-1]))-1].lower()
			elif(len(fastadict[id][pos-1:])<40): # print more seq before the snp
				extr=fastadict[id][pos-(81-len(fastadict[id][pos-1:])):pos-1].lower()
				b=len(extr)
				extr+=fastadict[id][pos-1:pos].upper()
				extr+=fastadict[id][pos:].lower()
			else: # print as much seq after as before the snp
				extr=fastadict[id][pos-40:pos-1].lower()
				extr+=fastadict[id][pos-1:pos].upper()
				extr+=fastadict[id][pos:pos+40].lower()
			print "SNP\t%s\t%d\t%s\t%s\t%s" %(id, pos, fields[3], fields[4], extr)
