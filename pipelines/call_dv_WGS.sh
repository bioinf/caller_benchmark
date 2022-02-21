#!/bin/bash

REFERENCE='GRCh37.primary_assembly.genome.fa'


for j in BWA BOWTIE2 BOWTIE2LOC ISAAC NOVOALIGN
do
	for i in *${j}.dedup.bam
	do
		docker run -v /media/array/callers_proj/raw_data/:/data google/deepvariant:latest /opt/deepvariant/bin/run_deepvariant --model_type=WGS --ref=/data/$REFERENCE --reads=/data/bams/WGS/$i --regions=/data/gencode.v19.revamp.bed --output_vcf=/data/${i%%.dedup.bam}/${i%%.dedup.bam}_DV_STANDART.vcf --output_gvcf=/data/${i%%.dedup.bam}/${i%%.dedup.bam}.DV.g.vcf --num_shards=6 &
	done
	wait
done
