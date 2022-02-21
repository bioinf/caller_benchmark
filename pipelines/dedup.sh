#!/bin/bash


for i in *.sorted.bam
do
	sleep 20
	while [ $( docker ps | grep broadinst | wc -l ) -ge 4 ] ; do sleep 1 ; done
	docker run -v ${PWD}:/data broadinstitute/gatk gatk --java-options "-Xmx24g" MarkDuplicates -I /data/$i -O /data/${i%%.sorted.bam}.dedup.bam -M /data/${i%%.sorted.bam}.MD.metrics & 
done
