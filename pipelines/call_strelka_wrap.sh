#!/bin/bash

#sleep 18000

for i in HG003_GENOME_{BWA,ISAAC,NOVOALIGN}*.bam HG004_GENOME_{BOWTIE2,BOWTIE2LOC,BWA,NOVOALIGN}*.bam HG005_GENOME_{BOWTIE2,BOWTIE2LOC,BWA}*.bam
do
	sleep 2
	while [ $( ps -AF | grep 'call_strelka.sh' | wc -l ) -ge 13 ] ; do sleep 1 ; done
	/media/array/callers_proj/raw_data/call_strelka.sh $i &
done
