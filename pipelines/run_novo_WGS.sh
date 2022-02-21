#!/bin/bash

NOVOALIGN='/media/array/callers_proj/raw_data/novocraft/novoalign'
REFERENCE='/media/array/callers_proj/raw_data/GRCh37.primary_assembly.genome.fa.nvlgn'
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
		"${NOVOALIGN}" -d "${REFERENCE}" -c 40 -f $R1 $R2 -o SAM "@RG\tID:${RGID}\tSM:S${TAG}\tLB:1\tPL:illumina" 2> ../../logs/${TAG}.novo.log | samtools view -bS - > ../../bams/WGS/${TAG}_${NODIR%%.fastq.gz}_NOVOALIGN.bam 
	done
	cd ../../
done

