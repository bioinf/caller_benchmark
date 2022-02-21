#!/bin/bash

sleep 14400

for i in *.bam
do
	samtools sort -@ 40 -o ${i%%.bam}.sorted.bam $i
done
