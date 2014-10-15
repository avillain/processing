#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

"""
fasta-stats.py: Utility script to count number of nucleotides/Aminoacids in a FASTA files.

Version 1.0
"""
from __future__ import with_statement
import sys
import argparse
import gzip
#Counter is used for per-character statistics
from collections import Counter

__author__    = "Uli Koehler & Anton Smirnov"
__copyright__ = "Copyright 2013 Uli Koehler"
__license__   = "Apache v2.0"

#def printSequenceStats(fileName, sequenceName, charOccurrenceMap, totalCharCount):
#    """
#    Print details about a sequence to stdout.
#    Called from inside parseFile().
#    
#    Keyword arguments:
#        sequenceName: The name of the sequence to print
#        charOccurrenceMap: A dict-like object that contains a char -->
#                            number of occurrences mapping
#        totalCharCount: The overall character count of the sequence
#    """
def printSequenceStats(index,sequenceName,totalCharCount):
    #print "Sequence '%s' from FASTA file '%s' contains %d sequence characters:" % (sequenceName, fileName, totalCharCount)
    print "%d\t%s\t%s" %(index,sequenceName,totalCharCount)
    #for char in sorted(charOccurrenceMap.iterkeys()):
    #    charCount = charOccurrenceMap[char]
    #    relativeFrequency = charCount * 100.0 / totalCharCount
    #    print "\t%s : %d = %f%%" % (char, charCount, relativeFrequency)
    #For nucleotide sequences (ATGC only), also print A+T vs G+C count
    #if sorted(charOccurrenceMap.iterkeys()) == ["A","C","G","T"]:
    #    #Print A+T count
    #    atCount = charOccurrenceMap["A"] + charOccurrenceMap["T"]
    #    atRelFrequency = atCount * 100.0 / totalCharCount
    #    print "\tA+T : %d = %f%%" % (atCount, atRelFrequency)
    #    #Print G+C count
    #    gcCount = charOccurrenceMap["G"] + charOccurrenceMap["C"]
    #    atRelFrequency = gcCount * 100.0 / totalCharCount
    #    print "\tG+C : %d = %f%%" % (gcCount, atRelFrequency)
        

def parseFile(index,filename):
    """
    Parse a FASTA fil and call printRe
    """
    #Set to the header line, with ">" removed from the beginning
    sequenceName = None
    #Key: character, value: Number of occurrences
    charOccurrenceMap = Counter()
    #The number of characters in the current sequence, not including \n
    charCount = 0
    #Keep track of consecutive comments, because they are appended
    previousLineWasComment = False
    #Open and iterate the file, auto- detect gzip
    openFunc = gzip.open if filename.endswith(".gz") else open
    with openFunc(filename, "r") as infile:
        for line in infile:
            line = line.strip()
            #Be super-compatible with the original specification
            if line.startswith(">"):
		index+=1
                #Process previous sequence, if any
                if sequenceName is not None:
                    printSequenceStats(index, sequenceName, charCount)
                    charOccurrenceMap = Counter()
                    charCount = 0
                #Take the entire comment line as (new) sequence ID (with ">" stripped)
                #Concatenate consecutive sequence lines
                if previousLineWasComment: #Append -- add one space between to normalize whitespace count
                    sequenceName += " " + line[1:].strip()
                else: 
                    sequenceName = line[1:].strip()
                previousLineWasComment = True
            else: #Line belongs to the sequence
                previousLineWasComment = False
                #Line has been stripped before, so we can count directly
                #Increment per-character stats (character occurrences)
                for char in line:
                    #Skip any character not in the whitelist, if whitelist (--only) is enabled
                    #if charWhitelist is not None and not char in charWhitelist:
                    #    continue
                    #We can only count after whitelite filter
                    charCount += 1
                    #In case-insensitive mode (default) count uppercased chars only
                    char = char.upper()
                    charOccurrenceMap[char] += 1
        #The last line has been read, print the last record, if any
        if sequenceName is not None:
	    index+=1
            printSequenceStats(index, sequenceName, charCount)
	return index

if __name__ == "__main__":
    #Allow single or multiple files to be specified
    parser = argparse.ArgumentParser(description='Compute simple statistics for FASTA files.')
    parser.add_argument('infiles', nargs='+', help='FASTA files (.fa, .fa.gz) to generate statistics for')
    args = parser.parse_args()
    #Process all FASTA files
    index=-1
    for filin in args.infiles:
        index=parseFile(index, filin)

