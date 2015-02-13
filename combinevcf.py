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
    parser = argparse.ArgumentParser(description='Takes vcf files from Pindel SAOPindel and PRISM and combines their output while checking if the position is covered in the control')
    parser.add_argument("pindel", help='Pindel vcf file')
    parser.add_argument("soap", help='SOAPindel vcf file')
    parser.add_argument("prism", help='PRISM vcf file')
    parser.add_argument("mpileup", help='Coverage mpileup file from the control')
    parser.add_argument("output",help='Output snpEffable file')

    if not sys.argv[1:] :
       sys.argv.append('-h')

    prog_name = os.path.basename(sys.argv[0])
    args = parser.parse_args()

#    if not args.known or not args.filtered or not args.raw :
#        parser.error("\nVcf file missing.")
#    else :
#        if not os.access(args.known, os.R_OK):
#            parser.error("Unable to read vcf file %s " %args.known)
#    	if not os.access(args.filtered, os.R_OK):
#            parser.error("Unable to read vcf file %s " %args.filtered)
#	if not os.access(args.raw, os.R_OK):
#            parser.error("Unable to read vcf file %s " %args.raw)

    cov={}
    with open(args.mpileup,"r") as mpifile:
	for line in mpifile:
		fields=line.split("\t")
		if fields[0] in cov.keys():
			cov[fields[0]][fields[1]]=fields[3]
		else:
			cov[fields[0]]={fields[1]:fields[3]}

    header=""
    vcfrecords={}
    with open(args.pindel,'r') as pinfile:
	for line in pinfile:
	    if line[:1]=='#':
		if line[1:6]!="CHROM":
		    header+=line
	    else:
		flds=line.split("\t") # info = field 8
		fields=[s.replace("\n", "") for s in flds]
		vcfrecords[tuple(fields[:2])]=[fields[2:],1,0] # code is +1 for pindel, +2 for soap and +4 for Prism ; eg 1=p 2=s 3=ps 4=P 5=pP 6=sP 7=psP, last field is coverage in control
    
    
    with open(args.soap,'r') as soafile:
	for line in soafile:
	    if line[:1]=='#':
		if line[1:6]!="CHROM":
		    header+=line
	    else:
		flds=line.split("\t")
		fields=[s.replace("\n", "") for s in flds]
		id=tuple(fields[:2])
		if tuple(id) in vcfrecords.keys():
		    vcfrecords[id][1]+=2
		    vcfrecords[id][0][1]=vcfrecords[id][0][1] if vcfrecords[id][0][1]==fields[3] else vcfrecords[id][0][1]+";%s"%(fields[3])
		    vcfrecords[id][0][2]=vcfrecords[id][0][2] if vcfrecords[id][0][2]==fields[4] else vcfrecords[id][0][2]+";%s"%(fields[4])
		    vcfrecords[id][0][3]="%s"%(fields[5]) if vcfrecords[id][0][3]=="." else vcfrecords[id][0][3]+";%s"%(fields[5])
		    vcfrecords[id][0][4]=vcfrecords[id][0][4] if fields[6]==vcfrecords[id][0][4] else vcfrecords[id][0][4]+";%s"%(fields[6])
		    vcfrecords[id][0][5]+=";%s"%(fields[7])
		    vcfrecords[id][0][6]+=";%s"%(fields[8])
		    vcfrecords[id][0][7]+=";%s"%(fields[9])
		else:
		    vcfrecords[tuple(fields[:2])]=[fields[2:],2,0]
    with open(args.prism,'r') as prifile:
	for line in prifile:
            if line[:1]=='#':
		if line[1:6]!="CHROM":
                    header+=line
            else:
                flds=line.split("\t")
		fields=[s.replace("\n", "") for s in flds]
                id=tuple(fields[:2])
                if tuple(id) in vcfrecords.keys():
                    vcfrecords[id][1]+=4
                    vcfrecords[id][0][1]=vcfrecords[id][0][1] if vcfrecords[id][0][1]==fields[3] else vcfrecords[id][0][1]+";%s"%(fields[3])
                    vcfrecords[id][0][2]=vcfrecords[id][0][2] if vcfrecords[id][0][2]==fields[4] else vcfrecords[id][0][2]+";%s"%(fields[4])
		    vcfrecords[id][0][5]+=":%s"%(fields[7])
                else:
                    vcfrecords[tuple(fields[:2])]=[fields[2:]+['.','.'],4,0]
    header+="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
    header2=header[:-1]+"\tPROG\tCOV\n"
    print header2,
    with open(args.output,"w") as filout:
	filout.write(header)
	for i in sorted(vcfrecords, key = lambda x: (x[0], int(x[1]))):
	    try:
                vcfrecords[i][2]=cov[i[0]][i[1]]
            except:
                pass
            print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"%(i[0],i[1],vcfrecords[i][0][0],vcfrecords[i][0][1],vcfrecords[i][0][2],vcfrecords[i][0][3],vcfrecords[i][0][4],vcfrecords[i][0][5],vcfrecords[i][0][6],vcfrecords[i][0][7],vcfrecords[i][1],vcfrecords[i][2])
            filout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(i[0],i[1],vcfrecords[i][0][0].split(";",1)[0],vcfrecords[i][0][1].split(";",1)[0],vcfrecords[i][0][2].split(";",1)[0],vcfrecords[i][0][3].split(";",1)[0],vcfrecords[i][0][4].split(";",1)[0],vcfrecords[i][0][5].split(";",1)[0],vcfrecords[i][0][6].split(";",1)[0],vcfrecords[i][0][7].split(";",1)[0]))
#	for i in vcfrecords:
#	    try:
#	        vcfrecords[i][2]=cov[i[0]][i[1]]
#	    except:
#	        pass
#	    print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"%(i[0],i[1],vcfrecords[i][0][0],vcfrecords[i][0][1],vcfrecords[i][0][2],vcfrecords[i][0][3],vcfrecords[i][0][4],vcfrecords[i][0][5],vcfrecords[i][0][6],vcfrecords[i][0][7],vcfrecords[i][1],vcfrecords[i][2])
#    	    filout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(i[0],i[1],vcfrecords[i][0][0],vcfrecords[i][0][1],vcfrecords[i][0][2],vcfrecords[i][0][3],vcfrecords[i][0]
#[4],vcfrecords[i][0][5],vcfrecords[i][0][6],vcfrecords[i][0][7]))
#    sys.exit(0)    
    
