#!/bin/bash

NOVOALIGN='/media/array/callers_proj/raw_data/novocraft/novoalign'
REFERENCE='/media/array/callers_proj/raw_data/GRCh37.primary_assembly.genome.fa'
DATENOW=$( date )
echo "Alignment started at: $DATENOW"

for FQDIR in fastqs/HG00{1,2,3,4,5,6,7}*_GENOME
do
	cd $FQDIR
	TAG=${FQDIR##fastqs/}
	echo "Aligning $FQDIR"
	for R1 in $( find * | grep -P '_R1_.*gz' )
	do
		R2=${R1/R1/R2}
		NODIR=$( basename $R1 )
		RGID=${NODIR%%_R1*}
		mkdir ../../bams/WGS/${TAG}_${NODIR}_ISAAC/ ../../fastqs/${TAG}/${NODIR}_ISAAC_input/
		ln -s ${PWD}/$R1 ../../fastqs/${TAG}/${NODIR}_ISAAC_input/lane1_read1.fastq.gz
		ln -s ${PWD}/$R2 ../../fastqs/${TAG}/${NODIR}_ISAAC_input/lane1_read2.fastq.gz
		isaac-align -r $REFERENCE -b ../../fastqs/${TAG}/${NODIR}_ISAAC_input/ -f fastq-gz --output-directory ../../bams/WGS/${TAG}_${NODIR}_ISAAC/ -j 40 -m 80 --mark-duplicates 0  --verbosity 4 --variable-read-length 1 &> ../../logs/${NODIR}_isaac.log
	done
	cd ../../
done

#TMPTAG=${FQ%%.R1.fastq.gz}
#	TAG=${TMPTAG##fastqs/}
#	SMP=$( echo $TAG | grep -oP '.*OME' )
#	mkdir bams/${SMP}_ISAAC/ fastqs/${SMP}_ISAAC_input/
#	ln -s ${PWD}/${TMPTAG}.R1.fastq.gz fastqs/${SMP}_ISAAC_input/lane1_read1.fastq.gz
#	ln -s ${PWD}/${TMPTAG}.R2.fastq.gz fastqs/${SMP}_ISAAC_input/lane1_read2.fastq.gz
#	echo "Aligning $TAG"
#	isaac-align -r $REFERENCE -b fastqs/${SMP}_ISAAC_input/ -f fastq-gz --output-directory bams/${SMP}_ISAAC -j 16 -m 80 --mark-duplicates 0  --verbosity 4 --variable-read-length 1 &> ./logs/${SMP}_isaac.log

