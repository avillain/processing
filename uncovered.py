#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

from pprint import pprint
import sys
import numpy

with open(sys.argv[1]) as filseq:
	print "coucou"

bases=[[]]
coverages=[[]]
chromosomes=[]
with open(sys.argv[2]) as filcov:
	for line in filcov:
		fields=line.split("\t")
		if not chromosomes:
			chromosomes.append(fields[0])
		elif fields[0]!=chromosomes[-1]:
			chromosomes.append(fields[0])
			bases.append([])
			coverages.append([])
		bases[len(chromosomes)-1].append(fields[2])
		coverages[len(chromosomes)-1].append(int(fields[3]))

means=[]
for i in range(len(chromosomes)):
	means.append(sum(coverages[i])/len(coverages[i]))

#pprint(means)
uncovered=[]
for j in range(len(chromosomes)):
	uncovered.append([])
	for k in range(len(coverages[j])):
		if coverages[j][k]<means[j]*1.0/10:
			uncovered[-1].append(k)

pprint(uncovered)
#pprint(chromosomes)
#pprint(coverages)
