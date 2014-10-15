#!/mount/biotools/bin/python2.7
# -*- coding: utf-8 -*-

import sys


try:
	maxdis=int(sys.argv[2])
except:
	maxdis=20

lastpos=[[0,""]]
with open(sys.argv[1], 'r') as snpfile:
	for line in snpfile:
			flag=0
			fields=line.split('\t')
			if fields[0]=='SNP':
				if lastpos[-1][0]==0:
					flag=1
					lastpos=[[int(fields[2]),line]]
				elif int(fields[2])-lastpos[-1][0] <= maxdis:
					lastpos.append([int(fields[2]),line])
				else:
					flag=1
					if len(lastpos)==1:
						print lastpos[0][1],
				if flag==1:
					lastpos=[[int(fields[2]),line]]
				

if len(lastpos)==1:
	print lastpos[0][1],

