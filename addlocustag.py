#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-

# Villain Adrien
# 2014     
# avillain@pasteur.fr

import re
import sys
import random
import string

def randomword(length):
   return ''.join(random.choice(string.lowercase) for i in range(length))

def addlocus(fil): # fichier .gbk
	recds= re.compile('CDS.*([0-9]+)\.\.([0-9]+)')
	relen= re.compile('^LOCUS.*')
	num=0
	with open(fil,'r') as filin:
    		for line in filin: # parcours du contenu du fichier
			print "%s" %line[:-1]
			if relen.findall(line):
				#id=line.split()[1]
				id=randomword(4).upper()
				num=0
			elif recds.findall(line):
				new="/locus_tag=\"%s\"\n" %(id+str(num))
				num+=1
				print "%s" %new[:-1].rjust(20+len(new))

if __name__ == '__main__':
	addlocus(sys.argv[1])
