#!/bin/bash

usage()
{
cat << EOF
usage: $0 options

This script run the SNP and Indels analysis for one genome using fqCleaner for reads pre-processing (optional), bwa mem for reads mapping, GATK2 for variant calling, and snpEff for variant annotation (optional)

OPTIONS:
   -h      Show this message
   -i      Input reads (.fastq format)
   -j	   Input mate reads (.fastq format)
   -r      Reference genome (.fasta format)
   -p      Output files prefix
   -n      No quality control performed
   -a      Annotation of variants with SnpEff (see -c and -d options)
   -c	   SnpEff config file (mandatory if -a)
   -d      SnpEff database ID of the species  (mandatory if -a)
   -v      Verbose

EXAMPLE:
./script.sh -i reads1.fq -j reads2.fq -r ref.fasta -p prefix -n -a -c SNPeffconfig1.txt -d spombe
EOF
}

INPUT=
INPUTMATE=
REFERENCE=
VERBOSE=
NOQC=
ANNOT=
CONF=
ID=
while getopts “hi:j:r:o:p:c:d:vna” OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         i)
             INPUT=$OPTARG
             ;;
	 j)
	     INPUTMATE=$OPTARG
	     ;;
         r)
             REFERENCE=$OPTARG
             ;;
         p)
             PREFIX=$OPTARG
             ;;
         c)
             CONF=$OPTARG
             ;;
         d)
             ID=$OPTARG
             ;;
         v)
             VERBOSE=1
             ;;
         n)
             NOQC=1
             ;;
         a)
             ANNOT=1
             ;;
         ?)
             usage
             exit
             ;;
     esac
done

echo "[info] parsing command-line"

if [[ -z $INPUT ]] || [[ -z $REFERENCE ]]
then
     echo "ERROR : Please supply a .fasta reference genome (-r) and .fastq reads to compare (-i)"
     usage
     exit 1
fi

if [[ -z $INPUTMATE ]]
then
	echo "[info] no mate reads file supplied, input is single reads"
fi

if [[ ! -z $ANNOT ]] && ( [[ -z $CONF ]] || [[ -z $ID ]] )
then
	echo "ERROR : SnpEff config file (-c) and SnpEff database ID (-d) are mandatory to perform SnpEff annotation (-a)"
	usage
	exit 1
fi

#basenames
refname=$(basename "$REFERENCE")
extension="${refname##*.}"
refname="${refname%.*}"
readgroup="@RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\PU:unit1"

#prefix of outputs
if [[ -z $PREFIX ]]
then
      dir=$(dirname $INPUT)
      PREFIX=${name%.*}"vs"$refname 
      echo "[info] no prefix supplied, output files will begin with $PREFIX"
else
      dir=$(dirname $PREFIX)
      PREFIX=$(basename $PREFIX)
      echo "[info] output files will begin with $PREFIX"
fi

TMP=$dir"/tmp"
RESULTS=$dir"/results"
echo "[info] working in $dir/tmp - results stored in $dir/results"
if [[ ! -s $RESULTS ]]
     then
          mkdir $RESULTS
fi

if [[ ! -s $TMP ]]
     then
          mkdir $TMP
fi

outputsam=$TMP"/"$PREFIX".sam"
outputbam=$TMP"/"$PREFIX".bam"
dedupbam=$TMP"/dedup_"$PREFIX".bam"
dictname=$RESULTS"/"$refname".dict"
metrics=$TMP"/"$refname".metrics"
realignedbam=$TMP"/realigned_"$PREFIX".bam"
rawvariants=$TMP"/"$PREFIX"_raw.vcf"

logfile=$TMP"/"$PREFIX".log"
REFNAME=$(basename $REFERENCE})
REFNAME="${REFNAME%.*}"
refgenome=$RESULTS"/"$REFNAME".fa"

if [ ! -f $refgenome ]
then 
      cp $REFERENCE $refgenome
fi

#Versions
echo "[info] versions : bwa/0.7.5a GenomeAnalysisTK/2.7-2 snpEff/3.5"
module load bwa/0.7.5a
module load GenomeAnalysisTK/2.7-2
module load snpEff/3.5

#Index reference genome
fai=$TMP"/"$REFERENCE".fai"
echo "[info] indexing reference genome"
if [ ! -f $fai ]
then
	samtools faidx $refgenome 
fi

if [ ! -f $RESULTS"/"$REFNAME".bwt" ]
then
        bwa index -a is -p $RESULTS"/"$REFNAME $refgenome 
fi

if [ ! -f $dictname ]
then
        CreateSequenceDictionary REFERENCE=$refgenome OUTPUT=$dictname 
fi

###SNPs
##filenames and filters
rawsnps=$TMP"/"$PREFIX"_snps_raw.vcf"
filtexpr="DP<$res || FS > 60.0 || MQ < 40.0"
filtname="cov"$res
filteredsnps=$TMP"/"$PREFIX"snps_cov"$res".vcf"
filteredcovsnps=$TMP"/"$PREFIX"snps_cov"$res"selected.vcf"
filteredhetsnps=$RESULTS"/"$PREFIX"snps_filtered.vcf"

rawindels=$TMP"/"$PREFIX"_indels_raw.vcf"
filteredindels=$RESULTS"/"$PREFIX"indels_filtered.vcf"
filtindel="QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"


#Annotation using snpEff
if [[ ! -z $ANNOT ]]
then
	configfile=$CONF
	refID=$ID
	echo "[info] annotation using snpEff : config file $CONF and ref ID $ID"
	#snps
	snpeffannotation=$TMP"/"$PREFIX"snps_effects.vcf"
	snpeffannotated=${filteredhetsnps%.*}"_annotated.vcf"
	snpstats=$TMP"/"$PREFIX"snpannot.csv"
	echo "[cmd] snpEff -c $configfile -v -s $snpstats -csvStats -o gatk $refID $filteredhetsnps > $snpeffannotation"
	snpEff -c $configfile -v -s $snpstats -csvStats -o gatk $refID $filteredhetsnps > $snpeffannotation
	GenomeAnalysisTK -T VariantAnnotator -R $refgenome -A SnpEff --variant $filteredhetsnps --snpEffFile $snpeffannotation  -o $snpeffannotated
	#indels
	indeleffannotation=$TMP"/"$PREFIX"indels_effects.vcf"
	indeleffannotated=${filteredindels%.*}"_annotated.vcf"
	indelstats=$TMP"/"$PREFIX"indelsannot.csv"
	echo "[cmd] snpEff -c $configfile -v -s $indelstats -o gatk $refID $filteredindels > $indeleffannotation"
	snpEff -c $configfile -v -s $indelstats -csvStats -o gatk $refID $filteredindels > $indeleffannotation
	GenomeAnalysisTK -T VariantAnnotator -R $refgenome -A SnpEff --variant $filteredindels --snpEffFile $indeleffannotation  -o $indeleffannotated
	
	allvariants=$RESULTS"/"$PREFIX"allvariants_annotated.vcf"
	echo "[info] merging snps and indels in the same .vcf file"
	cat $snpeffannotated > $allvariants
	grep -v '#' $indeleffannotated >> $allvariants
fi

