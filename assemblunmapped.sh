unmapped=${1%%.*}_unmapped.sam
fastq=${1%%.*}_unmapped.fastq

head -n 2 $1 > $unmapped
samtools view -S -f 4 $1 >> $unmapped

SamToFastq INPUT=$unmapped FASTQ=$fastq

velveth-101 assemblyunmapped 99 -short -fastq $fastq
velvetg-101 assemblyunmapped/ -cov_cutoff 4 -min_contig_lgth 200

