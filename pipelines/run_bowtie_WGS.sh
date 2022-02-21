#!/bin/bash

REFERENCE='/media/array/callers_proj/raw_data/GRCh37.primary_assembly.genome.fa.bwt2_idx'
DATENOW=$( date )
echo "Alignment started at: $DATENOW"

for FQDIR in fastqs/HG00{6,7}*_GENOME
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
		bowtie2 -x ${REFERENCE} -1 $R1 -2 $R2 -p 39 --rg-id $RGID --rg SM:$RGID --rg LB:LIBRARY --rg PL:ILLUMINA 2> ../../logs/${RGID}_bowtie2.log | samtools view -bS - > ../../bams/WGS/${TAG}_${NODIR%%.fastq.gz}_BOWTIE2.bam
#		bwa mem -M -t 40 -R "@RG\tID:$RGID\tSM:S${TAG}\tLB:1\tPL:illumina" \
#			$REFERENCE $R1 $R2 \
#			2> ../../logs/${TAG}.bwa.log | samtools view -bS - > ../../bams/WGS/${TAG}_${NODIR%%.fastq.gz}_BWA.bam 
	done
	cd ../../
done

DATENOW=$( date )
echo "Alignment finished at: $DATENOW"

