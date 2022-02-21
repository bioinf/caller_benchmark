#!/bin/bash

REFERENCE='/media/array/callers_proj/raw_data/GRCh37.primary_assembly.genome.fa.bwt2_idx'
DATENOW=$( date )
echo "Alignment started at: $DATENOW"

#sleep 21600

for FQ in fastqs/*R1.fastq.gz
do
	TMPTAG=${FQ%%.R1.fastq.gz}
	TAG=${TMPTAG##fastqs/}
	SMP=$( echo $TAG | grep -oP '.*OME' )
	echo "Aligning $TAG"
	bowtie2 -x ${REFERENCE} -1 ${TMPTAG}.R1.fastq.gz -2 ${TMPTAG}.R2.fastq.gz -p 40 --local --rg-id $SMP --rg SM:$SMP --rg LB:LIBRARY --rg PL:ILLUMINA 2> ./logs/${SMP}_bowtie2_local.log | samtools view -bS - > ./bams/${TAG}_BOWTIE2LOC.bam 
done

DATENOW=$( date )
echo "Alignment finished at: $DATENOW"

