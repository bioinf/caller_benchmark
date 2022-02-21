#!/bin/bash

REFERENCE='/media/array/callers_proj/raw_data/GRCh37.primary_assembly.genome.fa'
DATENOW=$( date )
echo "Alignment started at: $DATENOW"

for FQ in fastqs/HG00{6,7}*R1.fastq.gz
do
	TMPTAG=${FQ%%.R1.fastq.gz}
	TAG=${TMPTAG##fastqs/}
	SMP=$( echo $TAG | grep -oP '.*OME' )
	mkdir bams/${SMP}_ISAAC/ fastqs/${SMP}_ISAAC_input/
	ln -s ${PWD}/${TMPTAG}.R1.fastq.gz fastqs/${SMP}_ISAAC_input/lane1_read1.fastq.gz 
	ln -s ${PWD}/${TMPTAG}.R2.fastq.gz fastqs/${SMP}_ISAAC_input/lane1_read2.fastq.gz
	echo "Aligning $TAG"
	isaac-align -r $REFERENCE -b fastqs/${SMP}_ISAAC_input/ -f fastq-gz --output-directory bams/${SMP}_ISAAC -j 16 -m 80 --mark-duplicates 0  --verbosity 4 --variable-read-length 1 &> ./logs/${SMP}_isaac.log
done

DATENOW=$( date )
echo "Alignment finished at: $DATENOW"

