#!/usr/bin/python
# -*- coding: utf-8 -*-

# Villain Adrien
# 2015
# villain.adrien@gmail.com

import os
import re
import sys
import argparse
import warnings
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast import NCBIXML

def importdb(fildb):
    """
     Checks that a file is a DNA fasta file and its protein translation
     is ready to be used for homologous protein search
    """
    with open(fildb,'r') as f:
        l=f.readline()
        assert(l[:1]=='>'), "Database %s is not a fasta file.\n"%fildb
        seq=''
        while(1):
            l=f.readline()
            if l[:1]=='>':
                break
            else:
                seq+=l[:-1]
        assert(validate(seq)), "Database %s is not dna fasta.\n"%fildb
        protdb=os.path.splitext(fildb)[0]+'.prot'
        if os.path.isfile(protdb):
            print "6 frame translation already done.\n"
        else:
            print "Doing 6 frame translation using transeq.\n"
            try:
                return_code = subprocess.call("transeq -sformat pearson "
                                              "-frame=6 %s %s" %(fildb, 
                                              protdb), shell=True)
            except:
                print("Error calling transeq\n")
                raise
    return fildb, protdb


def validate(seq, alphabet='dna'):
## Source : https://www.biostars.org/p/102/ Giovanni M Dall'Olio 
    """ Check that a sequence only contains values from an alphabet """
    alphabets = {'dna': re.compile('^[acgtn]*$', re.I), 
             'protein': re.compile('^[acdefghiklmnpqrstvwy]*$', re.I)}
    if alphabets[alphabet].search(seq) is not None:
         return True
    else:
         return False

def filtertblout(filtblout, emaxglob, emaxdom):
    """ Filters significant hmmsearch hits """
    filout=os.path.splitext(filtblout)[0]+"_filtered.list"
    if os.path.isfile(filout):
        print "hmmsearch ouptut %s already filtered : %s\n" %(filtblout, filout)
    else:
        print "Filtering hmmsearch output %s in %s\n" %(filtblout, filout)
        try:
            with open(filout,"w") as fout:
                with open(filtblout,"r") as fin:
                    for lin in fin:
                        if lin[:1]!='#':
                            fields=lin.split()
                            eglob=float(fields[4])
                            edom=float(fields[7])
                            if eglob < emaxglob and edom < emaxdom:
                                fout.write("%s\t%s\t%s\n" %(fields[0][:-2],
                                                    fields[4],fields[7]))
        except:
            print "Error filtering hmmsearch output %s\n" %filtblout
            raise
    return filout

def fastaindex(filfasta):
    """ Calling fastaindex from exonerate """
    fastindex=os.path.splitext(filfasta)[0]+'.index'
    if os.path.isfile(fastindex):
        print "Fasta indexation already done : %s\n" %fastindex
    else:
        try:
            return_code = subprocess.call("fastaindex %s %s"%(
                                          filfasta,fastindex), shell=True)
        except:
            print "Error calling fastaindex on %s\n" %filfasta
    return fastindex

def ishmm3r(fil):
    """ Checks if a file is a HMM3R hmm profile file """
    with open(fil,'r') as f:
        li=f.readline()
        if li[:8]=="HMMER3/f":
            return True
        else:
            return False

def importprofiles(filprof):
    """ Reads a HMM3R file or list of files """
    if ishmm3r(filprof):
            return [filprof]
    else:
        print("File %s not a HMM3R hmm profile. Trying to read it "
              "as a list of HMM3R hmm profile files.\n" %filprof)
        listfiles=[]
        with open(filprof, 'r') as f:
            for line in f:
                if not os.access(line.strip(), os.R_OK):
                    print "Cannot read file %s\n" %line.strip()
                elif not ishmm3r(line.strip()):
                    print "File %s is not a HMM3R hmm profile.\n" %line.strip()
                else:
                    listfiles.append(line.strip())
        if listfiles:
            return listfiles
        else:
            sys.exit("No HMM3R hmm profile detected in %s. Aborting.\n"%filprof)

def search(db, profile, nthreads=1):
    """ Performs hmmsearch for given hmm profile on given protein fasta file """
    tblout=os.path.basename(os.path.splitext(profile)[0])+"_vs_"+\
                            os.path.basename(os.path.splitext(db)[0])+".tblout"
    out=os.path.basename(os.path.splitext(profile)[0])+"_vs_"+\
                         os.path.basename(os.path.splitext(db)[0])+".out"
    if os.path.isfile(tblout) and os.path.isfile(out):
        print "hmmsearch %s versus %s already performed\n" %(profile, db)
    else:
        print "Performing hmmsearch %s versus %s\n" %(profile, db)
        try:
            return_code = subprocess.call("hmmsearch --cpu %d --tblout %s %s %s"
                                          " > %s" %(nthreads, tblout, profile,
                                          db, out), shell=True)
        except:
            print "Error calling hmmsearch\n"
            raise
    return out, tblout

def fetchcandidate(cand, fasta, index, outfa):
    """ Fetching fasta record from file """
    print "Fetching fasta record %s from file %s\n" %(cand, fasta)
    try:
        return_code = subprocess.call("fastafetch %s %s %s >> %s"%(
                                      fasta,index,cand,outfa), shell=True)
    except:
        print "Error calling fastafetch on id %s\n"%cand
        raise
    return outfa

def blastorfcandidates(orfasta, nthreads=1):
    """ Runs blast to validate hmmsearch hits """
    blout=os.path.basename(os.path.splitext(orfasta)[0])+"_blast.xml"
    if os.path.isfile(blout):
        print "Validation of candidate ORFs by blast already performed\n"
    else:
        print "Running blast on %s to validate candidate ORFs\n" %orfasta
        try:
            return_code = subprocess.call(
          "blastp -outfmt 5 -max_target_seqs 1 -num_threads %d -evalue 0.1 -db "
          "uniref90.fasta -query %s > %s" %(
                                          nthreads, orfasta, blout), shell=True)
        except:
            print "Error running blast on file %s\n" %orfasta
            raise
    return blout

def getorfs(dnafafile):
    """ Writes open reading frames from dna fasta file to ORF fasta file """
    orfile=os.path.basename(os.path.splitext(dnafafile)[0])+"_ORFs.fa"
    if os.path.isfile(orfile):
        print "Retrieval of ORFs already performed : %s\n" %orfile
    else:
        print "Retrieving ORFs from  %s to validate candidates\n" %dnafafile
        try:
            table=1
            min_cds_len=150
            with open(orfile, 'w') as filout:
                for s in SeqIO.parse(dnafafile, "fasta"):
                    l=len(s.seq)
                    for strand, nuc in [(+1, s.seq),
                                        (-1, s.seq.reverse_complement())]:
                        for frame in xrange(3):
                            seq=str(nuc[frame:].translate(table))
                            for m in re.finditer(
                                        '[ARNDBCEQZGHILKMFPSTWYV]+', seq):
                                candcds=seq[m.start():m.end()]
                                if len(candcds) >= min_cds_len:
                                    fr=frame+1 if strand==1 else frame+4
                                    st, en = prot2dnacoords(seq, m.start(),
                                                            m.end(), fr)
                                    filout.write(">%s_%d_%d_%d\n%s\n"
                                                 %(s.id, fr, st, en, candcds))
        except:
            print "Error retrieving ORFs from file %s\n" %dnafafile
            raise
    return orfile

def prot2dnacoords(seq, start, end, frame): 
    """ Calculates dna coordinates from dna ORF and fully translated sequence
        frame = 1, 2, 3 (+1, +2, +3), 4, 5, 6 (-1, -2, -3) """
    if frame > 3:
        frame-=3
        st=3*(len(seq)-end+1)
        en=3*(len(seq)-start+1)
    else:
        st=3*start+frame
        en=3*(end+1) if end==len(seq) else 3*end
    en=en if en>=3*len(seq) else en+3
    return st, en


def parseblastorf(xmlout):
    """ Parses xml blastp ouptput to identify proteins """
    output=os.path.basename(os.path.splitext(xmlout)[0])+"_parsed.out"
    if os.path.isfile(output):
        print "Blast xml parsing already performed : %s\n" %output
    else:
        print "Parsing blast xml file  %s\n" %xmlout
        try:
            out=""
            with open(xmlout, 'r') as xmlin:
                blast_records = NCBIXML.parse(xmlin)
                lastcontig=""
                n=0
                for record in blast_records:
                    n+=1
                    res=re.match("([^_]+)_(\d+)_(\d+)_(\d+)", record.query)
                    contig=res.groups()[0]
                    st="+" if int(res.groups()[1]) <= 3 else "-"
                    s=res.groups()[2]
                    e=res.groups()[3]
                    if lastcontig!=contig:
                        out+="%s" %contig
                        lastcontig=contig
                        n=0
                    if record.alignments:
                        out+="\tORF%d\t%s\t%s\t%s\t%s\t%s\n" %(
              n+1, st, s, e, record.alignments[0].hit_def.replace("\n"," "),
              str(record.alignments[0].hsps[0].expect))
                    else:
                        out+="\tORF%d\t%s\t%s\t%s\tNo hits found\t\n" %(
                                                           n+1, st, s, e)
            with open(output, "w") as filout:
                filout.write(out)
        except:
            print "Error parsing blast xml output file %s\n" %xmlout
            raise
    return output


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                   description='Metagenomic search of marker genes')
    parser.add_argument("database", help='DNA fasta file to search')
    parser.add_argument("hmm", help='HMM3R hmm profile or list of hmm profiles')
    parser.add_argument("-e","--evalue", help='maximum e-value for hit',
                        default=1e-5)
    parser.add_argument("-n","--cpu", type=int, help='Number of threads',
                        default=1)
    
    if not sys.argv[1:] :
       sys.argv.append('-h')
    
    prog_name = os.path.basename(sys.argv[0])
    args = parser.parse_args()
    
    if not args.database :
        parser.error("\nDatabase file missing.")
    else :
        if not os.access(args.database, os.R_OK):
            parser.error("Unable to read file %s " %args.database)
    
    if not args.hmm :
        parser.error("\nHMM file missing.")
    else :
        if not os.access(args.hmm, os.R_OK):
            parser.error("Unable to read file %s " %args.hmm)
    
    ### check fasta, hmm, e-values
    evaldomain=10*args.evalue
    fadna,fapro=importdb(args.database)
    faind=fastaindex(args.database)
    hmmprofiles=importprofiles(args.hmm)

    ### search fasta database, filter hits
    tblouts=[]
    for hmm in hmmprofiles:
        o,tblo=search(fapro, hmm, args.cpu)
        tblouts.append(tblo)
    listmatches=[]
    for tb in tblouts:
        match=filtertblout(tb, args.evalue, evaldomain)
        listmatches.append(match)

    ### extract candidate contigs, ORFs
    candidates=[]
    for m in listmatches:
        with open(m, 'r') as fcand:
            for line in fcand:
                idee=line.split()[0]
                if idee not in candidates:
                    candidates.append(idee)
    candfa=os.path.basename(os.path.splitext(args.database)[0])+"_candidates.fa"
    if not os.path.isfile(candfa):
        for c in candidates:
            res=fetchcandidate(c, fadna, faind, candfa)
    else:
        print "Fasta fetch for %s candidates already done\n" %args.database
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        filorfs=getorfs(candfa)
    blorfs=blastorfcandidates(filorfs, args.cpu)
    results=parseblastorf(blorfs)
    print "Results stored in %s\n" %results 

