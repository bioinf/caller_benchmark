#!/bin/bash

REFERENCE='/media/array/callers_proj/raw_data/GRCh37.primary_assembly.genome.fa'
DATENOW=$( date )
echo "Alignment started at: $DATENOW"

sleep 25200

for FQDIR in fastqs/HG001_GENOME
do
	cd $FQDIR
	TAG=${FQDIR##fastqs/}
	echo "Aligning $FQDIR"
	for R1 in $( find * | grep -P '_R1_.*gz' )
	do
		R2=${R1/R1/R2}
		NODIR=$( basename $R1 )
		RGID=${NODIR%%_R1*}
		echo $R1 $R2 $RGID
		bwa mem -M -t 40 -R "@RG\tID:$RGID\tSM:S${TAG}\tLB:1\tPL:illumina" \
			$REFERENCE $R1 $R2 \
			2> ../../logs/${TAG}.bwa.log | samtools view -bS - > ../../bams/WGS/${TAG}_${NODIR%%.fastq.gz}_BWA.bam 
	done
	cd ../../
done

DATENOW=$( date )
echo "Alignment finished at: $DATENOW"

