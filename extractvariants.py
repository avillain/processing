#!/mount/biotools/bin/python2.7
# -*- coding: utf-8 -*-

import sys
import argparse
from pprint import pprint
import re

parser = argparse.ArgumentParser(description='Extracts flanking regions of SNPS/indels.')
parser.add_argument("reference", help='Reference file in FASTA format')
parser.add_argument("snpfile", help='SNP/indel file in VCF format')
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

with open(args.snpfile, 'r') as snp:
	print "Type\tContig\tPos\tRef\tVar\tFlanking seq in ref"
	for line in snp:
		if line[:1]!="#":
			fields=line.split("\t")
			id=fields[0]
			pos=int(fields[1])
			leng=len(fields[3])-1
			typ="SNP" if leng==0 else "INDEL"
			##### we have to be sure that the contig is not too short to extract
			if(len(fastadict[id])<=80) : # then just output the whole thing
				extr=fastadict[id][:pos-1].lower()
				extr+=fastadict[id][pos-1:pos+leng].upper()
				extr+=fastadict[id][pos+leng:].lower()
				print "cas 1 :",len(extr)
			elif(len(fastadict[id][:pos-1])<40): # print more seq after the snp
				extr=fastadict[id][:pos-1].lower()
				extr+=fastadict[id][pos-1:pos+leng].upper()
				extr+=fastadict[id][pos+leng:pos+(80-len(fastadict[id][:pos-1])-leng)-1].lower()
				print "cas 2 :",len(extr)
			elif(len(fastadict[id][pos-1:])<40): # print more seq before the snp
				extr=fastadict[id][pos-(81-len(fastadict[id][pos-1:])-leng):pos-1].lower()
				extr+=fastadict[id][pos-1:pos+leng].upper()
				extr+=fastadict[id][pos+leng:].lower()
				print "cas 3 :",len(extr)
			else: # print as much seq after as before the snp
				extr=fastadict[id][pos-40:pos-1].lower()
				extr+=fastadict[id][pos-1:pos+leng].upper()
				extr+=fastadict[id][pos+leng:pos+40-leng].lower()
				print "cas 4 :",len(extr)
			print "%s\t%s\t%d\t%s\t%s\t%s" %(typ,id, pos, fields[3], fields[4], extr)

