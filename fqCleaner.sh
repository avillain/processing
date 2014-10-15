#!/bin/bash
#$ -S /bin/bash 

#################################################################################
#                                                                               #
# fqCleaner : a BASH script to filter reads from FASTQ formatted files          #
# Copyright (C) 2011  Alexis Criscuolo, Ghislaine Guigon and Thierry Hieu       #
#                                                                               #
  VERSION=0.5.0
#                                                                               #
# fqCleaner is  free software;  you can  redistribute it  and/or modify it      #
# under the terms of  the GNU General Public  License as published  by the      #
# Free Software Foundation;  either version 2 of the License,  or (at your      #
# option) any later version.                                                    #
#                                                                               #
# fqCleaner is distributed in the hope that it will be useful, but WITHOUT      # 
# ANY WARRANTY;  without even the  implied warranty of  MERCHANTABILITY or      #
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for      #
# more details.                                                                 #
#                                                                               #
# You should  have received  a copy  of the  GNU  General  Public  License      #
# along   with  this   program;   if  not,  write  to  the  Free  Software      #
# Foundation Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307            #
#                                                                               # 
# Contact:                                                                      #
#  PF8 - Genotypage des Pathogenes et Sante Publique                            #
#  INSTITUT  PASTEUR                                                            #
#  28 rue du Dr Roux - 75724 Paris  (France)                                    #
#                                                                               #
# alexis.criscuolo@pasteur.fr                                                   #
#                                                                               #
#################################################################################

#################################################################################
#                                                                               #
# ================                                                              #
# = INSTALLATION =                                                              #
# ================                                                              #
#                                                                               #
# Prior to any launch, first verify that the following pathes to each required  #
# binary are correct:                                                           #
#
# sickle : to   trim   off    low   quality   bases   in    both   read   ends  #
#          (https://github.com/najoshi/sickle)                                  #

  SICKLE_TRIMMER=/local/gensoft/scripts/sickle

# fastq_quality_filter : to filter out non-confident reads (part of the FASTX-  #
#                        Toolkit: hannonlab.cshl.edu/fastx_toolkit/)            #

  FASTQ_QUALITY_FILTER=/local/gensoft/scripts/fastq_quality_filter

# fastx_artifacts_filter : to filter out artifactual reads (part of the FASTX-  #
#                          Toolkit: hannonlab.cshl.edu/fastx_toolkit/)          #

  FASTX_ARTIFACTS_FILTER=/local/gensoft/scripts/fastx_artifacts_filter

# AlienTrimmer : to  trim  off  alien  oligonucleotides   in  both  read  ends  #
#                (ftp://ftp.pasteur.fr/pub/GenSoft/projects/AlienTrimmer/);     #
#                by default,  alien oligonucleotide sequences are written into  #
#                file $ALIEN_SEQUENCES,  but it is  possible to input  its own  #
#                alien sequence file with option -c                             #

  ALIEN_TRIMMER=/local/gensoft/scripts/AlienTrimmer
  ALIEN_SEQUENCES=/pasteur/projets/specific/PFBAG_ngs/Analyse/avillain/scripts/alienTrimmerPF8contaminants.fasta

# fqconvert : to convert any fastq-formatted file into sanger-encoded file      #
#             (ftp://ftp.pasteur.fr/pub/GenSoft/projects/fqtools/)              #

  FQ2SANGER=/local/gensoft/scripts/fqconvert
                                                                               
# fqextract : to  extract  specified   reads  from  a   fastq  formatted  file  #
#             (ftp://ftp.pasteur.fr/pub/GenSoft/projects/fqtools/)              #

  FQEXTRACT=/local/gensoft/scripts/fqextract

# fqduplicate : to obtain the read names that are present several times within  #
#               a FASTQ formated file                                           #
#               (ftp://ftp.pasteur.fr/pub/GenSoft/projects/fqtools/)            #

  FQDUPLICATE=/local/gensoft/scripts/fqduplicate

# Secondly,  give the execute  permission on the  script fqCleaner.sh by using  #
# the following command:                                                        #
#                                                                               #
#   chmod +x fqCleaner.sh                                                       #
#                                                                               #
#################################################################################
 
#################################################################################
#                                                                               #
# =============                                                                 #
# = EXECUTION =                                                                 #
# =============                                                                 #
#                                                                               #
# You can launch fqCleaner with the following command:                          #
#                                                                               #
#   ./fqCleaner.sh [options]                                                    #
#                                                                               #
# Available options:                                                            #
#                                                                               #
# -f <infile>  FASTQ formatted input file name (mandatory option)               #
#                                                                               #
# -r <infile>  when using  paired-ends data,  this option  allows inputing the  #
#              second file (i.e. reverse reads)                                 #
#                                                                               #
# -q <int>     quality score  threshold (default: 20);  all bases with quality  #
#              score below this threshold are considered as non-confident       #
#                                                                               #
# -l <int>     minimum  required length for a read  (default: 30); all trimmed  #
#              reads with length below this threshold are filtered out          #
#                                                                               #
# -p <int>     minimum percent of bases (default: 80) that must have a quality  #
#              score higher than  the fixed threshold  (set by the option -q);  #
#              all reads  that do  not verify  this  confidence  criterion are  #
#              filtered out                                                     #
#                                                                               #
# -s <string>  a sequence of tasks to be iteratively performed, each being set  #
#              by one of the following letters:                                 #
#                Q: each read is trimmed off to  remove non-confident bases at  #
#                   5' and 3' ends (quality score threshold is set with option  #
#                   -q); all short reads are discarded (minimum read length is  #
#                   set with option -l)                                         #
#                F: each read containing  too few confident  bases is filtered  #
#                   out (the minimum percentage of confident bases is set with  #
#                   option -p)                                                  #
#                A: each artefactual read is discarded                          #
#                C: contaminant oligonucleotide sequences are trimmed off when  #
#                   occuring in either  5' or 3' ends of each read;  all short  #
#                   reads  are  discarded  (minimum  read  length is  set with  #
#                   option -l);  user contaminant  sequences can  be set  with  #
#                   option -c                                                   #
#                D: all duplicated reads are removed                            #
#                d: (only for  paired-ends data)  all duplicated  reads within  #
#                   each input file are removed; this leads to unpaired reads   #
#                N: the number of (remaining) reads is displayed                #
#              (default sequence of tasks: NQNCNFNANDN)                         #
#                                                                               #
# -c <infile>  sequence file containing contaminant nucleotide sequences to be  #
#              trimmed off during step 'C'                                      #
#                                                                               #
# -x <file>    single-end: output file name;  paired-ends: forward read output  #
#              file name (default: <option -f>.<option -s>.fq )                 #
#                                                                               #
# -y <file>    paired-ends data only: reverse read output file name  (default:  #
#              <option -r>.<option -s>.fq)                                      #
#                                                                               #
# -z <file>    paired-ends data only:  single read output file name  (default:  #
#              <option -f>.<option -s>.sgl.fq)                                  #
#                                                                               #
#################################################################################


#################################################################################
# Help displaying                                                               #
#################################################################################

if [ "$1" = "-?" ] || [ $# -le 1 ]
then
  echo "" ;
  echo " fqCleaner v.$VERSION" ;
  echo "" ;
  echo " USAGE :" ;
  echo "    fqCleaner.sh [options]" ;
  echo "  where 'options' are :" ;
  echo "   -f <infile>  FASTQ formatted input file name (mandatory option)               "
  echo "   -r <infile>  when using  paired-ends data,  this option  allows inputing the  "
  echo "                second file (i.e. reverse reads)                                 "
  echo "   -q <int>     quality score  threshold (default: 20);  all bases with quality  "
  echo "                score below this threshold are considered as non-confident       "
  echo "   -l <int>     minimum required length for a read (default: 30)                 "
  echo "   -p <int>     minimum percent of bases (default: 80) that must have a quality  "
  echo "                score higher than the fixed threshold (set by option -q)         "
  echo "   -s <string>  a sequence of tasks to be iteratively performed, each being set  "
  echo "                one of the following letters:                                    "
  echo "                  Q: each read is trimmed off to  remove non-confident bases at  "
  echo "                     5' and 3' ends (quality score threshold is set with option  "
  echo "                     -q); all short reads are discarded (minimum read length is  "
  echo "                     set with option -l)                                         "
  echo "                  F: each read containing  too few confident  bases is filtered  "
  echo "                     out (the minimum percentage of confident bases is set with  "
  echo "                     option -p)                                                  "
  echo "                  A: each artefactual read is filtered out                       "
  echo "                  C: contaminant oligonucleotide sequences are trimmed off when  "
  echo "                     occuring in either  5' or 3' ends of each read;  all short  "
  echo "                     reads  are  discarded  (minimum  read  length is  set with  "
  echo "                     option -l);  user contaminant  sequences can  be set  with  "
  echo "                     option -c                                                   " 
  echo "                  D: all duplicated single- or paired-ends reads are removed     "
  echo "                  d: (only for  paired-ends data)  all duplicated  reads within  "
  echo "                     each input file are removed                                 "
  echo "                  N: the number of (remaining) reads is displayed                "
  echo "                (default sequence of tasks: NQCFNDN)                             "
  echo "   -c <infile>  sequence file containing contaminant nucleotide sequences  (one  "
  echo "                per line) to be trimmed off during step 'C'                      "
  echo "   -x <file>    single-end: output file name;  paired-ends: forward read output  "
  echo "                file name (default: <option -f>.<option -s>.fq )                 "
  echo "   -y <file>    paired-ends data only: reverse read output file name  (default:  "
  echo "                <option -r>.<option -s>.fq)                                      "
  echo "   -z <file>    paired-ends data only:  single read output file name  (default:  "
  echo "                <option -f>.<option -s>.sgl.fq)                                  "
  echo "" ;
  echo "  Examples :                                                                     " 
  echo "" ;
  echo "   Single-ends  reads inside  a file  named 'reads.fq',  default  options,  and  "
  echo "   output file named 'reads.clean.fq':                                           "
  echo "     fqCleaner.sh  -f reads.fq  -x reads.clean.fq                                "
  echo "" ;
  echo "   Paired-ends reads  inside files  named  'fwd.fq' and 'rev.fq',  minimum read  "
  echo "   length of 50 bps, quality score threshold of 30:                              "
  echo "     fqCleaner.sh  -f fwd.fq  -r rev.fq  -l 50  -q 30                            "
  echo "" ;
  exit
fi
 

#################################################################################
# Data                                                                          #
#################################################################################

PHRED=("!"  "\"" "#"  "$"  "%"  "&"  "'"  "("  ")"  "*"        # 0-9
       "+"  ","  "-"  "."  "/"  "0"  "1"  "2"  "3"  "4"        # 10-19
       "5"  "6"  "7"  "8"  "9"  ":"  ";"  "<"  "="  ">"        # 20-29
       "?"  "@"  "A"  "B"  "C"  "D"  "E"  "F"  "G"  "H"  "I"); # 30-40
# for q in $(seq 0 40); do echo "$q ${PHRED[$q]}"; done


#################################################################################
# Functions                                                                     #
#################################################################################

# randomfile  $INFILE
#
#   INFILE: file name, this will returns the randomly generated file name INFILE.RANDOM
#
randomfile() {
  rdmf=$1.$RANDOM; while [ -e $rdmf ]; do rdmf=$1.$RANDOM ; done
  echo $rdmf ;
}

# fqsize  $INFILE
#
#   INFILE: input file name in fastq format
#
fqsize() {
  if [ -e $1 ] 
  then echo $(( $(wc -l $1 | sed 's/ .*//g') / 4 ));
  else echo "0";
  fi
}

# fq2sanger  $INFILE  $OUTFILE  $LOGFILE
#
#   INFILE:  input file name in fastq format (with extension .fq)
#   OUTFILE: output file name converted from INFILE with sanger quality scores
#   LOGFILE: output file name containing INFILE format info (append)
#
fq2sanger() {
  fq2sanger_tmp=$(randomfile $1);
  $FQ2SANGER -o $2 $1 2>$fq2sanger_tmp ; 
  head -1 $fq2sanger_tmp | sed 's/Found//g' >> $3; rm $fq2sanger_tmp ;
  if [ ! -e $2 ]; then cp $1 $2; fi
}

# fqtrim_se  $INFILE  $QUALITY_THRESHOLD  $MIN_READ_LENGTH  $OUTFILE
#
#   INFILE:            input file name in fastq and sanger format
#   QUALITY_THRESHOLD: integer between 0 and 40
#   MIN_READ_LENGTH:   integer higher than 0
#   OUTFILE:           output file name containing trimmed single-end reads
#
fqtrim_se() {
  $SICKLE_TRIMMER se --quiet -f $1 -t sanger -q $2 -l $3 -o $4 ;
}

# fqtrim_pe  $INFILE_FWD  $INFILE_REV  $QUALITY_THRESHOLD  $MIN_READ_LENGTH  $OUTFILE_FWD  $OUTFILE_REV  $OUTFILE_SGL
#
#   INFILE_FWD:        input file name containing forward reads in fastq and sanger format
#   INFILE_REV:        input file name containing reverse reads in fastq and sanger format
#   QUALITY_THRESHOLD: integer between 0 and 40
#   MIN_READ_LENGTH:   integer higher than 0
#   OUTFILE_FWD:       output file name containing trimmed forward reads
#   OUTFILE_REV:       output file name containing trimmed reverse reads
#   OUTFILE_SGL:       output file name containing single reads
#
fqtrim_pe() {
  $SICKLE_TRIMMER pe --quiet -f $1 -r $2 -t sanger -q $3 -l $4 -o $5 -p $6 -s $7 ;
}

# fqfilter  $INFILE  $QUALITY_THRESHOLD  $PERCENT_CONFIDENT_BPS  $OUTFILE
#
#   INFILE:                input file name in fastq and sanger format
#   QUALITY_THRESHOLD:     integer between 0 and 40
#   PERCENT_CONFIDENT_BPS: integer between 0 and 100
#   OUTFILE:               output file name containing filtered single-end reads
#
fqfilter() {
  $FASTQ_QUALITY_FILTER -Q 33 -i $1 -q $2 -p $3 -o $4 ;
}

# fqartfilter  $INFILE  $OUTFILE
#
#   INFILE:  input file name in fastq or fasta format
#   OUTFILE: output file name containing single-end reads without artifacts
#
fqartfilter() {
  $FASTX_ARTIFACTS_FILTER -Q 33 -i $1 -o $2 ;
}

# fqalienclipper_se  $INFILE  $ALIEN_FILE  $QUALITY_THRESHOLD  $PERCENT_CONFIDENT_BPS  $MIN_READ_LENGTH  $OUTFILE
#
#   INFILE:                input file name in fastq and sanger format
#   ALIEN_FILE:            infile name containing alien oligonucleotide sequences
#   QUALITY_THRESHOLD:     integer between 0 and 40
#   PERCENT_CONFIDENT_BPS: integer between 0 and 100
#   MIN_READ_LENGTH:       integer higher than 0
#   OUTFILE:               output file name containing trimmed single-end reads
#
fqalienclipper_se() {
  at_tmp=$(randomfile $1);
  $ALIEN_TRIMMER -i $1 -c $2 -q ${PHRED[$3]} -p $4 -l $5 -k 10 -o $6 > $at_tmp ;
  rm $at_tmp ;
}

# fqalienclipper_pe  $INFILE_FWD  $INFILE_REV  $ALIEN_FILE  $QUALITY_THRESHOLD  $PERCENT_CONFIDENT_BPS  $MIN_READ_LENGTH  $OUTFILE_FWD  $OUTFILE_REV  $OUTFILE_SGL
#
#   INFILE_FWD:            input file name containing forward reads in fastq and sanger format
#   INFILE_REV:            input file name containing reverse reads in fastq and sanger format
#   ALIEN_FILE:            file name containing alien oligonucleotide sequences
#   QUALITY_THRESHOLD:     integer between 0 and 40
#   PERCENT_CONFIDENT_BPS: integer between 0 and 100
#   MIN_READ_LENGTH:       integer higher than 0
#   OUTFILE_FWD:           output file name containing trimmed forward reads
#   OUTFILE_REV:           output file name containing trimmed reverse reads
#   OUTFILE_SGL:           output file name containing single reads
#
fqalienclipper_pe() {
  at_tmp=$(randomfile $1);
  $ALIEN_TRIMMER -if $1 -ir $2 -cf $3 -cr $3 -q ${PHRED[$4]} -p $5 -l $6 -k 10 -of $7 -or $8 -os $9 > $at_tmp ;
  rm $at_tmp ;
}

# fqnoduplicate_se  $INFILE  $OUTFILE
#
#   INFILE:  input file name in fastq format
#   OUTFILE: output file name containing single-end reads with no duplicate
#
fqnoduplicate_se() {
  fqnoduplicate_se_tmp=$(randomfile $1);
  $FQDUPLICATE $1 > $fqnoduplicate_se_tmp ;
  $FQEXTRACT -px -l $fqnoduplicate_se_tmp $1 > $2 ; rm $fqnoduplicate_se_tmp ;
}

# fqnoduplicate_pe  $INFILE_FWD  $INFILE_REV  $OUTFILE_FWD  $OUTFILE_REV
#
#   INFILE_FWD:  input file name containing forward reads in fastq and sanger format
#   INFILE_REV:  input file name containing reverse reads in fastq and sanger format
#   OUTFILE_FWD: output file name containing forward reads with no paired duplicate
#   OUTFILE_REV: output file name containing reverse reads with no paired duplicate
#
fqnoduplicate_pe() {
  fqnoduplicate_pe_tmp=$(randomfile $1);
  $FQDUPLICATE $1 $2 > $fqnoduplicate_pe_tmp ;
  $FQEXTRACT -px -l $fqnoduplicate_pe_tmp $1 > $3 ;
  $FQEXTRACT -px -l $fqnoduplicate_pe_tmp $2 > $4 ;
  rm $fqnoduplicate_pe_tmp
}

# fqintersect  $INFILE_FWD  $INFILE_REV  $OUTFILE_FWD  $OUTFILE_REV  $OUTFILE_SGL  [RA]
#
#   INFILE_FWD:  input file name containing forward reads in fastq and sanger format
#   INFILE_REV:  input file name containing reverse reads in fastq and sanger format
#   OUTFILE_FWD: output file name containing each forward read for which the reverse read exists
#   OUTFILE_REV: output file name containing each reverse read for which the forward read exists
#   OUTFILE_SGL: output file name containing single reads 
#   [RA]:        if 'R', then all outfiles are replaced; if 'A', all outfiles are appended
#
fqintersect() {
  sed -n '1~4p' $1 > $1.tmp ; sed 's/^@//g' $1.tmp > $1.idf ; rm $1.tmp ; grep -o ":.*[/# ]" $1.idf > $1.id ; sort $1.id -o $1.ids ; rm $1.id ;
  sed -n '1~4p' $2 > $2.tmp ; sed 's/^@//g' $2.tmp > $2.idf ; rm $2.tmp ; grep -o ":.*[/# ]" $2.idf > $2.id ; sort $2.id -o $2.ids ; rm $2.id ;
  comm -2 -3 $1.ids $2.ids > $1.idu ; comm -1 -3 $1.ids $2.ids > $2.idu ; rm $1.ids $2.ids ;
  grep -F -f $1.idu $1.idf > $1.idfsgl ; rm $1.idu $1.idf ; 
  grep -F -f $2.idu $2.idf > $2.idfsgl ; rm $2.idu $2.idf ;
  if [ "$6" = "R" ]; then if [ -e $3 ]; then rm $3 ; fi; if [ -e $4 ]; then rm $4 ; fi; if [ -e $5 ]; then rm $5 ; fi; fi
  $FQEXTRACT -px -l $1.idfsgl $1 >> $3 ;
  $FQEXTRACT -px -l $2.idfsgl $2 >> $4 ;
  $FQEXTRACT -p -l $1.idfsgl $1 >> $5 ;
  $FQEXTRACT -p -l $2.idfsgl $2 >> $5 ;
  rm $1.idfsgl $2.idfsgl ;
}

# gettime() $START
#
#   START: the starting time in seconds
#
gettime() {
  t=$(( $SECONDS - $1 )); sec=$(( $t % 60 )); min=$(( $t / 60 ));
  if [ $sec -lt 10 ]; then sec="0$sec"; fi
  if [ $min -lt 10 ]; then min="0$min"; fi
  echo "[$min:$sec]" ;
}


#################################################################################
# Beginning fqCleaner                                                           #
#################################################################################

###############################
#  Reading options            #
###############################
FASTQ_FWD="XXX"
FASTQ_REV="XXX"
OUTFASTQ_FWD="XXX"
OUTFASTQ_REV="XXX"
OUTFASTQ_SGL="XXX"
STEPS="NQCFNDN"; # STEPS="NQNCNFNANDN"
QUALITY_THRESHOLD=20
MIN_READ_LENGTH=30
PERCENT_CONFIDENT_BPS=80;
while getopts f:r:c:x:y:z:s:q:l:p: option
do
  case $option in
  f)
    FASTQ_FWD="$OPTARG"
    if [ ! -e $FASTQ_FWD ]; then echo "   problem with input file (option -f): '$FASTQ_FWD' does not exist" ; exit ; fi
   ;;
  r)
    FASTQ_REV="$OPTARG"
    if [ ! -e $FASTQ_REV ]; then echo "   problem with input file (option -r): '$FASTQ_REV' does not exist" ; exit ; fi
    PAIRED_ENDS="true"
   ;;
  c)
    ALIEN_SEQUENCES="$OPTARG"
    if [ ! -e $ALIEN_SEQUENCES ]; then echo "   problem with input file (option -c): '$ALIEN_SEQUENCES' does not exist" ; exit ; fi
   ;;
  x)
    OUTFASTQ_FWD="$OPTARG"
    #if [ -e $OUTFASTQ_FWD ]
    #then echo "   [WARNING] the specified output file $OUTFASTQ_FWD exists and will be removed"; rm -i $OUTFASTQ_FWD ; if [ -e $OUTFASTQ_FWD ]; then exit; fi
    #fi
   ;;
  y)
    OUTFASTQ_REV="$OPTARG"
    #if [ -e $OUTFASTQ_REV ]
    #then echo "   [WARNING] the specified output file $OUTFASTQ_REV exists and will be removed"; rm -i $OUTFASTQ_REV ; if [ -e $OUTFASTQ_REV ]; then exit; fi
    #fi
   ;;
  z)
    OUTFASTQ_SGL="$OPTARG"
    #if [ -e $OUTFASTQ_SGL ]
    #then echo "   [WARNING] the specified output file $OUTFASTQ_SGL exists and will be removed"; rm -i $OUTFASTQ_SGL ; if [ -e $OUTFASTQ_SGL ]; then exit; fi
    #fi
   ;;
  s)
    STEPS="$OPTARG"
   ;;
  q)
    QUALITY_THRESHOLD=$OPTARG
    if [ $QUALITY_THRESHOLD -lt 0 ] || [ $QUALITY_THRESHOLD -gt 40 ]; then echo "   the quality score threshold must range from 0 to 40 (option -q)" ; exit ; fi
   ;;
  l)
    MIN_READ_LENGTH=$OPTARG
    if [ $MIN_READ_LENGTH -le 0 ]; then echo "   the read length threshold must be a positive integer (option -l)" ; exit ; fi
   ;;
  p)
    PERCENT_CONFIDENT_BPS=$OPTARG
    if [ $PERCENT_CONFIDENT_BPS -lt 0 ] || [ $PERCENT_CONFIDENT_BPS -gt 100 ]; then echo "   the minimum percent of confident bases must range from 0 to 100 (option -p)" ; exit ; fi
   ;;
  esac
done

#############################################
#  Verifying mandatory options              #
#############################################
if [ "$FASTQ_FWD" = "XXX" ]; then echo "   no FASTQ formatted input file (mandatory option -f)" ; exit ; fi
#############################################
#  Reading FASTQ cleaning steps to perform  #
#############################################
# step X <=> step CFQ
# step Y <=> step CQ
STEPS=$(echo "$STEPS" | sed -e 's/CFQ\|CQF\|FCQ\|FQC\|QCF\|QFC/X/g' | sed -e 's/CQ\|QC/Y/g')

OUT=""; STEP_NB=${#STEPS}; STEP_ID=0;
while [ $STEP_ID -lt $STEP_NB ]
do
  STEP=${STEPS:$STEP_ID:1};
  if [ "$FASTQ_REV" = "XXX" ]
  then
    if   [ "$STEP" = "Q" ] || [ "$STEP" = "F" ] || [ "$STEP" = "A" ] || [ "$STEP" = "C" ] || [ "$STEP" = "D" ]; then OUT="$OUT$STEP"; 
    elif [ "$STEP" = "X" ]; then OUT="$OUT""CFQ";
    elif [ "$STEP" = "Y" ]; then OUT="$OUT""CQ";
    fi
  else
    if [ "$STEP" = "Q" ] || [ "$STEP" = "F" ] || [ "$STEP" = "A" ] || [ "$STEP" = "C" ] || [ "$STEP" = "D" ] || [ "$STEP" = "d" ]; then OUT="$OUT$STEP"; 
    elif [ "$STEP" = "X" ]; then OUT="$OUT""CFQ";
    elif [ "$STEP" = "Y" ]; then OUT="$OUT""CQ";
    fi
  fi
  STEP_ID=$(( $STEP_ID + 1 ));
done
if [ "$STEPS" = "N" ]; then OUT="N"; fi
#############################################
#  Computing read lengths                   #
#############################################
seq=$(head -2 $FASTQ_FWD | tail -1); READ_LENGTH=${#seq}; 
#############################################
#  Starting timer                           #
#############################################
STIME=$SECONDS;


#################################################################################
# Launching FASTQ cleaning                                                      #
#################################################################################

###############################
#  Single-ends                #
###############################
if [ "$FASTQ_REV" = "XXX" ]
then
  NAME=${FASTQ_FWD%.*}; INFILE=$(randomfile $NAME); OUTFILE=$(randomfile $NAME) ;
  LOGFILE=$NAME.log.txt; if [ -e $LOGFILE ]; then rm $LOGFILE; fi
  
  echo -n "$(gettime $STIME)   Reading file $FASTQ_FWD ..." ;
  #fq2sanger  $FASTQ_FWD  $INFILE  $LOGFILE ;
  cp $FASTQ_FWD $INFILE
  if [ ! -e $INFILE ]; then echo -e "\033[31m[fail]\033[00m" ; exit ; fi
  echo " [ok]" ;
  
  # echo "  read lengths = $READ_LENGTH" >> $LOGFILE ;
  echo "Launching fqCleaner v.$VERSION:" >> $LOGFILE ;
  echo "  quality score threshold                 -q $QUALITY_THRESHOLD" >> $LOGFILE ;
  echo "  minimum read length                     -l $MIN_READ_LENGTH" >> $LOGFILE ;
  echo "  minimum percentage of confident bases   -p $PERCENT_CONFIDENT_BPS" >> $LOGFILE ;

  STEP_NB=${#STEPS}; STEP_ID=0;
  while [ $STEP_ID -lt $STEP_NB ]
  do
    STEP=${STEPS:$STEP_ID:1};
    case $STEP in
    N)
      echo -n "> number of reads:";
      card=$(fqsize $INFILE); echo " $card" ;
      echo "> number of reads: $card" >> $LOGFILE ;
     ;;
    Q)
      echo -n "$(gettime $STIME)   Trimming off low-quality bases ..." ;
      fqtrim_se  $INFILE  $QUALITY_THRESHOLD  $MIN_READ_LENGTH  $OUTFILE ; 
      if [ ! -e $OUTFILE ]; then echo -e "\033[31m[fail]\033[00m" ; echo "Trimming off low-quality bases ... [fail]" >> $LOGFILE ; rm $INFILE ; exit ; fi
      mv $OUTFILE $INFILE ; echo " [ok]" ;
      echo "Trimming off low-quality bases ... [ok]" >> $LOGFILE ;
     ;;
    F)
      echo -n "$(gettime $STIME)   Filtering out non-confident reads ..." ;
      fqfilter  $INFILE  $QUALITY_THRESHOLD  $PERCENT_CONFIDENT_BPS  $OUTFILE ; 
      if [ ! -e $OUTFILE ]; then echo -e "\033[31m[fail]\033[00m" ; echo "Filtering out non-confident reads ... [fail]" >> $LOGFILE ; rm $INFILE ; exit ; fi
      mv $OUTFILE $INFILE ; echo " [ok]" ;
      echo "Filtering out non-confident reads ... [ok]" >> $LOGFILE  ;
     ;;
    A)
      echo -n "$(gettime $STIME)   Filtering out artifactual reads ..." ;
      fqartfilter  $INFILE  $OUTFILE ; 
      if [ ! -e $OUTFILE ]; then echo -e "\033[31m[fail]\033[00m" ; echo "Filtering out artifactual reads ... [fail]" >> $LOGFILE ; rm $INFILE ; exit ; fi
      mv $OUTFILE $INFILE ; echo " [ok]" ;
      echo "Filtering out artifactual reads ... [ok]" >> $LOGFILE ;
     ;;
    C)
      echo -n "$(gettime $STIME)   Trimming off contaminant oligonucleotides ..." ;
      fqalienclipper_se  $INFILE  $ALIEN_SEQUENCES  0  0  $MIN_READ_LENGTH  $OUTFILE ; 
      if [ ! -e $OUTFILE ]; then echo -e "\033[31m[fail]\033[00m" ; echo "Trimming off contaminant oligonucleotides ... [fail]" >> $LOGFILE ; rm $INFILE ; exit ; fi
      mv $OUTFILE $INFILE ; echo " [ok]" ;
      echo "Trimming off contaminant oligonucleotides ... [ok]" >> $LOGFILE ;
     ;;
    X)
      echo -n "$(gettime $STIME)   Trimming off low-quality and contaminant residues, and filtering out non-confident reads ..." ;
      fqalienclipper_se  $INFILE  $ALIEN_SEQUENCES  $QUALITY_THRESHOLD  $PERCENT_CONFIDENT_BPS  $MIN_READ_LENGTH  $OUTFILE ; 
      if [ ! -e $OUTFILE ]; then echo -e "\033[31m[fail]\033[00m" ; echo "Trimming off low-quality and contaminant residues, and filtering out non-confident reads ... [fail]" >> $LOGFILE ; rm $INFILE ; exit ; fi
      mv $OUTFILE $INFILE ; echo " [ok]" ;
      echo "Trimming off low-quality and contaminant residues, and filtering out non-confident reads ... [ok]" >> $LOGFILE ;
     ;;
    Y)
      echo -n "$(gettime $STIME)  Trimming off low-quality bases and contaminant oligonucleotides  ..." ;
      fqalienclipper_se  $INFILE  $ALIEN_SEQUENCES  $QUALITY_THRESHOLD  0  $MIN_READ_LENGTH  $OUTFILE ; 
      if [ ! -e $OUTFILE ]; then echo -e "\033[31m[fail]\033[00m" ; echo "Trimming off low-quality and contaminant residues ... [fail]" >> $LOGFILE ; rm $INFILE ; exit ; fi
      mv $OUTFILE $INFILE ; echo " [ok]" ;
      echo "Trimming off low-quality and contaminant residues ... [ok]" >> $LOGFILE ;
     ;;
    D)
      echo -n "$(gettime $STIME)   Filtering out duplicated reads ..." ;
      fqnoduplicate_se  $INFILE  $OUTFILE ; 
      if [ ! -e $OUTFILE ]; then echo -e "\033[31m[fail]\033[00m" ; echo "Filtering out duplicated reads ... [fail]" >> $LOGFILE ; rm $INFILE ; exit ; fi
      mv $OUTFILE $INFILE ; echo " [ok]" ;
      echo "Filtering out duplicated reads ... [ok]" >> $LOGFILE ;
     ;;
    esac

    STEP_ID=$(( $STEP_ID + 1 ));
  done

  if [ "$OUTFASTQ_FWD" = "XXX" ]; then OUTFASTQ_FWD=$NAME.$OUT.fq; fi
  mv $INFILE $OUTFASTQ_FWD ;
  echo "$(gettime $STIME)   Remaining single-ends reads written into file $OUTFASTQ_FWD" ;
  chmod 775 $OUTFASTQ_FWD ;
  echo "Remaining single-ends reads written into file $OUTFASTQ_FWD" >> $LOGFILE ;
  echo "Total running time: $(gettime $STIME)" >> $LOGFILE ;
  chmod 775 $LOGFILE ;


###############################
#  Paired-ends                #
###############################
else
  NAME_FWD=${FASTQ_FWD%.*}; INFILE_FWD=$(randomfile $NAME_FWD); OUTFILE_FWD=$(randomfile $NAME_FWD) ;
  NAME_REV=${FASTQ_REV%.*}; INFILE_REV=$(randomfile $NAME_REV); OUTFILE_REV=$(randomfile $NAME_REV) ;
  NAME_SGL=$NAME_FWD; OUTFILE_SGL=$(randomfile $NAME_SGL) ;
  LOGFILE=$NAME_FWD.log.txt; if [ -e $LOGFILE ]; then rm $LOGFILE; fi
  
  echo -n "$(gettime $STIME)   Reading file $FASTQ_FWD ..." ;
  fq2sanger  $FASTQ_FWD  $INFILE_FWD  $LOGFILE ;
  if [ ! -e $INFILE_FWD ]; then echo -e "\033[31m[fail]\033[00m" ; exit ; fi
  echo " [ok]" ;
  echo -n "$(gettime $STIME)   Reading file $FASTQ_REV ..." ;
  fq2sanger  $FASTQ_REV  $INFILE_REV  $LOGFILE ;
  if [ ! -e $INFILE_REV ]; then echo -e "\033[31m[fail]\033[00m" ; exit ; fi
  echo " [ok]" ;

  # echo "  read lengths = $READ_LENGTH" >> $LOGFILE ;
  echo "Launching fqCleaner v.$VERSION:" >> $LOGFILE ;
  echo "  quality score threshold                 -q $QUALITY_THRESHOLD" >> $LOGFILE ;
  echo "  minimum read length                     -l $MIN_READ_LENGTH" >> $LOGFILE ;
  echo "  minimum percentage of confident bases   -p $PERCENT_CONFIDENT_BPS" >> $LOGFILE ;

  pe="true";

  STEP_NB=${#STEPS}; STEP_ID=0;
  while [ $STEP_ID -lt $STEP_NB ]
  do
    STEP=${STEPS:$STEP_ID:1};
    case $STEP in
    N)
      echo -n "> number of reads:"; 
      card_fwd=$(fqsize $INFILE_FWD);  echo -n "  fwd=$card_fwd";
      if [ $STEP_ID -eq 0 ] || [ "$pe" = "false" ]
      then card_rev=$(fqsize $INFILE_REV);  echo -n "  rev=$card_rev";
      else card_rev=$card_fwd;              echo -n "  rev=$card_rev";
      fi
      card_sgl=$(fqsize $OUTFILE_SGL); echo    "  sgl=$card_sgl";
      echo "> number of reads:  fwd=$card_fwd  rev=$card_rev  sgl=$card_sgl" >> $LOGFILE ;
      if [ $card_fwd -ne $card_rev ]; then pe="false"; fi
     ;;
    Q)
      if [ "$pe" = "false" ]
      then 
        echo -n "$(gettime $STIME)   Intersecting reads ..." ;
        fqintersect  $INFILE_FWD  $INFILE_REV  $OUTFILE_FWD  $OUTFILE_REV  $OUTFILE_SGL  A ; 
	mv $OUTFILE_FWD $INFILE_FWD ; mv $OUTFILE_REV $INFILE_REV ; 
        echo " [ok]" ;
        pe="true";
        echo "Intersecting reads ... [ok]" >> $LOGFILE ;
      fi
      echo -n "$(gettime $STIME)   Trimming off low-quality bases ..." ;
      tmpfile=$(randomfile $NAME_SGL) ;
      fqtrim_pe  $INFILE_FWD  $INFILE_REV  $QUALITY_THRESHOLD  $MIN_READ_LENGTH  $OUTFILE_FWD  $OUTFILE_REV  $tmpfile ;
      if [ ! -e $OUTFILE_FWD ]; then echo -e "\033[31m[fwd fail]\033[00m" ; echo "Trimming off low-quality bases ... [fwd fail]" >> $LOGFILE ; rm $INFILE_FWD $INFILE_REV $OUTFILE_SGL ; exit ; fi
      if [ ! -e $OUTFILE_REV ]; then echo -e "\033[31m[rev fail]\033[00m" ; echo "Trimming off low-quality bases ... [rev fail]" >> $LOGFILE ; rm $INFILE_FWD $INFILE_REV $OUTFILE_SGL ; exit ; fi
      if [ ! -e $tmpfile ];     then echo -e "\033[31m[sgl fail]\033[00m" ; echo "Trimming off low-quality bases ... [sgl fail]" >> $LOGFILE ; rm $INFILE_FWD $INFILE_REV $OUTFILE_SGL ; exit ; fi
      mv $OUTFILE_FWD $INFILE_FWD ; mv $OUTFILE_REV $INFILE_REV ; cat $tmpfile >> $OUTFILE_SGL ; rm $tmpfile ; echo " [ok]" ;
      echo "Trimming off low-quality bases ... [ok]" >> $LOGFILE ;
     ;;
    F)
      echo -n "$(gettime $STIME)   Filtering out non-confident reads ." ;
      fqfilter  $INFILE_FWD  $QUALITY_THRESHOLD  $PERCENT_CONFIDENT_BPS  $OUTFILE_FWD ; 
      if [ ! -e $OUTFILE_FWD ]; then echo -e ".. \033[31m[fwd fail]\033[00m" ; echo "Filtering out non-confident reads ... [fwd fail]" >> $LOGFILE ; rm $INFILE_FWD $INFILE_REV $OUTFILE_SGL ; exit ; fi
      mv $OUTFILE_FWD $INFILE_FWD ;
      echo -n "." ;
      fqfilter  $INFILE_REV  $QUALITY_THRESHOLD  $PERCENT_CONFIDENT_BPS  $OUTFILE_REV ; 
      if [ ! -e $OUTFILE_REV ]; then echo -e ". \033[31m[rev fail]\033[00m" ; echo "Filtering out non-confident reads ... [rev fail]" >> $LOGFILE ; rm $INFILE_FWD $INFILE_REV $OUTFILE_SGL ; exit ; fi
      mv $OUTFILE_REV $INFILE_REV ;
      echo ". [ok]" ;
      pe="false";
      echo "Filtering out non-confident reads ... [ok]" >> $LOGFILE ;
     ;;
    A)
      echo -n "$(gettime $STIME)   Filtering out artifactual reads ." ;
      fqartfilter  $INFILE_FWD  $OUTFILE_FWD ; 
      if [ ! -e $OUTFILE_FWD ]; then echo -e ".. \033[31m[fwd fail]\033[00m" ; echo "Filtering out artifactual reads ... [fwd fail]" >> $LOGFILE ; rm $INFILE_FWD $INFILE_REV $OUTFILE_SGL ; exit ; fi
      mv $OUTFILE_FWD $INFILE_FWD ;
      echo -n "." ;
      fqartfilter  $INFILE_REV  $OUTFILE_REV ; 
      if [ ! -e $OUTFILE_REV ]; then echo -e ". \033[31m[rev fail]\033[00m" ; echo "Filtering out artifactual reads ... [rev fail]" >> $LOGFILE ; rm $INFILE_FWD $INFILE_REV $OUTFILE_SGL ; exit ; fi
      mv $OUTFILE_REV $INFILE_REV ;
      echo ". [ok]" ;
      pe="false";
      echo "Filtering out artifactual reads ... [ok]" >> $LOGFILE ;
     ;;
    C)
      if [ "$pe" = "false" ]
      then 
        echo -n "$(gettime $STIME)   Intersecting reads ..." ;
        fqintersect  $INFILE_FWD  $INFILE_REV  $OUTFILE_FWD  $OUTFILE_REV  $OUTFILE_SGL  A ; 
	mv $OUTFILE_FWD $INFILE_FWD ; mv $OUTFILE_REV $INFILE_REV ; 
        echo " [ok]" ;
        pe="true";
        echo "Intersecting reads ... [ok]" >> $LOGFILE ;
      fi
      echo -n "$(gettime $STIME)   Trimming off contaminant oligonucleotides ..." ;
      tmpfile=$(randomfile $NAME_SGL) ;
      fqalienclipper_pe  $INFILE_FWD  $INFILE_REV  $ALIEN_SEQUENCES  0  0  $MIN_READ_LENGTH  $OUTFILE_FWD  $OUTFILE_REV  $tmpfile ;
      if [ ! -e $OUTFILE_FWD ]; then echo -e "\033[31m[fwd fail]\033[00m" ; echo "Trimming off contaminant oligonucleotides ... [fwd fail]" >> $LOGFILE ; rm $INFILE_FWD $INFILE_REV $OUTFILE_SGL ; exit ; fi
      if [ ! -e $OUTFILE_REV ]; then echo -e "\033[31m[rev fail]\033[00m" ; echo "Trimming off contaminant oligonucleotides ... [rev fail]" >> $LOGFILE ; rm $INFILE_FWD $INFILE_REV $OUTFILE_SGL ; exit ; fi
      if [ ! -e $tmpfile ];     then echo -e "\033[31m[sgl fail]\033[00m" ; echo "Trimming off contaminant oligonucleotides ... [sgl fail]" >> $LOGFILE ; rm $INFILE_FWD $INFILE_REV $OUTFILE_SGL ; exit ; fi
      mv $OUTFILE_FWD $INFILE_FWD ; mv $OUTFILE_REV $INFILE_REV ; cat $tmpfile >> $OUTFILE_SGL ; rm $tmpfile ;
      echo " [ok]" ;
      echo "Trimming off contaminant oligonucleotides ... [ok]" >> $LOGFILE ;
     ;;
    X)
      if [ "$pe" = "false" ]
      then 
        echo -n "$(gettime $STIME)   Intersecting reads ..." ;
        fqintersect  $INFILE_FWD  $INFILE_REV  $OUTFILE_FWD  $OUTFILE_REV  $OUTFILE_SGL  A ; 
	mv $OUTFILE_FWD $INFILE_FWD ; mv $OUTFILE_REV $INFILE_REV ; 
        echo " [ok]" ;
        pe="true";
        echo "Intersecting reads ... [ok]" >> $LOGFILE ;
      fi
      echo -n "$(gettime $STIME)   Trimming off low-quality and contaminant residues, and filtering out non-confident reads ..." ;
      tmpfile=$(randomfile $NAME_SGL) ;
      fqalienclipper_pe  $INFILE_FWD  $INFILE_REV  $ALIEN_SEQUENCES  $QUALITY_THRESHOLD  $PERCENT_CONFIDENT_BPS  $MIN_READ_LENGTH  $OUTFILE_FWD  $OUTFILE_REV  $tmpfile ;
      if [ ! -e $OUTFILE_FWD ]; then echo -e "\033[31m[fwd fail]\033[00m" ; echo "Trimming off low-quality and contaminant residues, and filtering out non-confident reads ... [fwd fail]" >> $LOGFILE ; rm $INFILE_FWD $INFILE_REV $OUTFILE_SGL ; exit ; fi
      if [ ! -e $OUTFILE_REV ]; then echo -e "\033[31m[rev fail]\033[00m" ; echo "Trimming off low-quality and contaminant residues, and filtering out non-confident reads ... [rev fail]" >> $LOGFILE ; rm $INFILE_FWD $INFILE_REV $OUTFILE_SGL ; exit ; fi
      if [ ! -e $tmpfile ];     then echo -e "\033[31m[sgl fail]\033[00m" ; echo "Trimming off low-quality and contaminant residues, and filtering out non-confident reads ... [sgl fail]" >> $LOGFILE ; rm $INFILE_FWD $INFILE_REV $OUTFILE_SGL ; exit ; fi
      mv $OUTFILE_FWD $INFILE_FWD ; mv $OUTFILE_REV $INFILE_REV ; cat $tmpfile >> $OUTFILE_SGL ; rm $tmpfile ;
      echo " [ok]" ;
      echo "Trimming off low-quality and contaminant residues, and filtering out non-confident reads ... [ok]" >> $LOGFILE ;
     ;;
    Y)
      if [ "$pe" = "false" ]
      then 
        echo -n "$(gettime $STIME)   Intersecting reads ..." ;
        fqintersect  $INFILE_FWD  $INFILE_REV  $OUTFILE_FWD  $OUTFILE_REV  $OUTFILE_SGL  A ; 
	mv $OUTFILE_FWD $INFILE_FWD ; mv $OUTFILE_REV $INFILE_REV ; 
        echo " [ok]" ;
        pe="true";
        echo "Intersecting reads ... [ok]" >> $LOGFILE ;
      fi
      echo -n "$(gettime $STIME)   Trimming off low-quality and contaminant  ..." ;
      tmpfile=$(randomfile $NAME_SGL) ;
      fqalienclipper_pe  $INFILE_FWD  $INFILE_REV  $ALIEN_SEQUENCES  $QUALITY_THRESHOLD  0  $MIN_READ_LENGTH  $OUTFILE_FWD  $OUTFILE_REV  $tmpfile ;
      if [ ! -e $OUTFILE_FWD ]; then echo -e "\033[31m[fwd fail]\033[00m" ; echo "Trimming off low-quality and contaminant reads ... [fwd fail]" >> $LOGFILE ; rm $INFILE_FWD $INFILE_REV $OUTFILE_SGL ; exit ; fi
      if [ ! -e $OUTFILE_REV ]; then echo -e "\033[31m[rev fail]\033[00m" ; echo "Trimming off low-quality and contaminant reads ... [rev fail]" >> $LOGFILE ; rm $INFILE_FWD $INFILE_REV $OUTFILE_SGL ; exit ; fi
      if [ ! -e $tmpfile ];     then echo -e "\033[31m[sgl fail]\033[00m" ; echo "Trimming off low-quality and contaminant reads ... [sgl fail]" >> $LOGFILE ; rm $INFILE_FWD $INFILE_REV $OUTFILE_SGL ; exit ; fi
      mv $OUTFILE_FWD $INFILE_FWD ; mv $OUTFILE_REV $INFILE_REV ; cat $tmpfile >> $OUTFILE_SGL ; rm $tmpfile ;
      echo " [ok]" ;
      echo "Trimming off low-quality and contaminant reads ... [ok]" >> $LOGFILE ;
     ;;
    D)
      if [ "$pe" = "false" ]
      then 
        echo -n "$(gettime $STIME)   Intersecting reads ..." ;
        fqintersect  $INFILE_FWD  $INFILE_REV  $OUTFILE_FWD  $OUTFILE_REV  $OUTFILE_SGL  A ; 
	mv $OUTFILE_FWD $INFILE_FWD ; mv $OUTFILE_REV $INFILE_REV ; 
        echo " [ok]" ;
        pe="true";
        echo "Intersecting reads ... [ok]" >> $LOGFILE ;
      fi
      echo -n "$(gettime $STIME)   Filtering out duplicated reads ..." ;
      fqnoduplicate_pe  $INFILE_FWD  $INFILE_REV  $OUTFILE_FWD  $OUTFILE_REV ;
      if [ ! -e $OUTFILE_FWD ]; then echo -e "\033[31m[fwd fail]\033[00m" ; echo "Filtering out duplicated reads ... [fwd fail]" >> $LOGFILE ; rm $INFILE_FWD $INFILE_REV $OUTFILE_SGL ; exit ; fi
      if [ ! -e $OUTFILE_REV ]; then echo -e "\033[31m[rev fail]\033[00m" ; echo "Filtering out duplicated reads ... [rev fail]" >> $LOGFILE ; rm $INFILE_FWD $INFILE_REV $OUTFILE_SGL ; exit ; fi
      mv $OUTFILE_FWD $INFILE_FWD ; mv $OUTFILE_REV $INFILE_REV ; 
      echo " [ok]" ;
      echo "Filtering out duplicated reads ... [ok]" >> $LOGFILE ;
     ;;
    d)
      echo -n "$(gettime $STIME)   Filtering out duplicated reads ." ;
      fqnoduplicate_se  $INFILE_FWD  $OUTFILE_FWD ; 
      if [ ! -e $OUTFILE_FWD ]; then echo -e ".. \033[31m[fwd fail]\033[00m" ; echo "Filtering out duplicated reads ... [fwd fail]" >> $LOGFILE ; rm $INFILE_FWD $INFILE_REV $OUTFILE_SGL ; exit ; fi
      mv $OUTFILE_FWD $INFILE_FWD ; 
      echo -n "." ;
      fqnoduplicate_se  $INFILE_REV  $OUTFILE_REV ; 
      if [ ! -e $OUTFILE_REV ]; then echo -e ". \033[31m[rev fail]\033[00m" ; echo "Filtering out duplicated reads ... [rev fail]" >> $LOGFILE ; rm $INFILE_FWD $INFILE_REV $OUTFILE_SGL ; exit ; fi
      mv $OUTFILE_REV $INFILE_REV ; 
      echo ". [ok]" ;
      echo "Filtering out duplicated reads ... [ok]" >> $LOGFILE ;
      pe="false";
     ;;
    esac

    STEP_ID=$(( $STEP_ID + 1 ));
  done

  if [ "$pe" = "false" ]
  then 
    echo -n "$(gettime $STIME)   Intersecting reads ..." ;
    fqintersect  $INFILE_FWD  $INFILE_REV  $OUTFILE_FWD  $OUTFILE_REV  $OUTFILE_SGL  A ; 
    mv $OUTFILE_FWD $INFILE_FWD ; mv $OUTFILE_REV $INFILE_REV ; 
    echo " [ok]" ;
    echo "Intersecting reads ... [ok]" >> $LOGFILE ;

    STEP_ID=$(( $STEP_ID - 1 )); STEP=${STEPS:$STEP_ID:1};
    if [ "$STEP" = "N" ]
    then 
      echo -n "> number of reads:"; 
      card_fwd=$(fqsize $INFILE_FWD);  echo -n "  fwd=$card_fwd";
      card_rev=$(fqsize $INFILE_REV);  echo -n "  rev=$card_rev";
      card_sgl=$(fqsize $OUTFILE_SGL); echo    "  sgl=$card_sgl";
      echo "> number of reads:  fwd=$card_fwd  rev=$card_rev  sgl=$card_sgl" >> $LOGFILE ;
    fi
  fi

  if [ "$OUTFASTQ_FWD" = "XXX" ]; then OUTFASTQ_FWD=$NAME_FWD.$OUT.fq; fi
  if [ "$OUTFASTQ_REV" = "XXX" ]; then OUTFASTQ_REV=$NAME_REV.$OUT.fq; fi
  if [ "$OUTFASTQ_SGL" = "XXX" ]; then OUTFASTQ_SGL=$NAME_SGL.$OUT.sgl.fq; fi
  mv $INFILE_FWD $OUTFASTQ_FWD ;
  mv $INFILE_REV $OUTFASTQ_REV ;
  if [ ! -e $OUTFILE_SGL ]; then echo > $OUTFILE_SGL; fi
  mv $OUTFILE_SGL $OUTFASTQ_SGL ;

  echo "$(gettime $STIME)   Remaining paired-ends reads written into files $OUTFASTQ_FWD and $OUTFASTQ_REV" ;
  chmod 775 $OUTFASTQ_FWD ; chmod 775 $OUTFASTQ_REV ;
  echo "$(gettime $STIME)   Single reads written into file $OUTFASTQ_SGL" ;
  chmod 775 $OUTFASTQ_SGL ;
  echo "Remaining paired-ends reads written into files $OUTFASTQ_FWD and $OUTFASTQ_REV" >> $LOGFILE ;
  echo "Single reads written into file $OUTFASTQ_SGL" >> $LOGFILE ;
  echo "Total running time: $(gettime $STIME)" >> $LOGFILE ;
  chmod 775 $LOGFILE ;

fi

