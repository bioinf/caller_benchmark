#!/bin/bash

REFERENCE='/media/array/callers_proj/raw_data/GRCh37.primary_assembly.genome.fa'
DATENOW=$( date )
echo "Alignment started at: $DATENOW"

#sleep 36000

for FQ in fastqs/HG00{6,7}*R1.fastq.gz
do
	TMPTAG=${FQ%%.R1.fastq.gz}
	TAG=${TMPTAG##fastqs/}
	SMP=$( echo $TAG | grep -oP '.*OME' )
	echo "Aligning $TAG"
	bwa mem -M -t 36 -R "@RG\tID:$TAG\tSM:S${SMP}\tLB:1\tPL:illumina" \
		$REFERENCE ${TMPTAG}.R1.fastq.gz ${TMPTAG}.R2.fastq.gz \
		2> ./logs/${TAG}.bwa.log | samtools view -bS - > ./bams/${TAG}_BWA.bam 
done

DATENOW=$( date )
echo "Alignment finished at: $DATENOW"

