#!/bin/bash

NOVOALIGN='/media/array/callers_proj/raw_data/novocraft/novoalign'
REFERENCE='/media/array/callers_proj/raw_data/GRCh37.primary_assembly.genome.fa.nvlgn'
DATENOW=$( date )
echo "Alignment started at: $DATENOW"

for FQ in fastqs/HG00{6,7}*R1.fastq.gz
do
	TMPTAG=${FQ%%.R1.fastq*}
	TAG=${TMPTAG##fastqs/}
	SMP=$( echo $TAG | grep -oP '.*OME' )
	echo "Aligning $TAG"
	"${NOVOALIGN}" -d "${REFERENCE}" -c 40 -f $FQ ${FQ%%R1*}R2.fastq.gz -o SAM "@RG\tID:${TAG}\tSM:S${SMP}\tLB:1\tPL:illumina" 2> ./logs/${TAG}.novo.log | samtools view -bS  -> ./bams/${TAG}_NOVOALIGN.bam 
done

DATENOW=$( date )
echo "Alignment finished at: $DATENOW"

