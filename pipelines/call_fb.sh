#!/bin/bash

REFERENCE='GRCh37.primary_assembly.genome.fa'
CALLINGREGIONS='gencode.v19.revamp.bed'

i=$1

# FreeBayes
freebayes --standard-filters -t /media/array/callers_proj/raw_data/$CALLINGREGIONS -f /media/array/callers_proj/raw_data/$REFERENCE $i > /media/array/callers_proj/raw_data/${i%%.dedup.bam}/${i%%.dedup.bam}.FB.raw.vcf

docker run -v /media/array/callers_proj/raw_data:/data broadinstitute/gatk:latest gatk LeftAlignAndTrimVariants \
	-V /data/${i%%.dedup.bam}/${i%%.dedup.bam}.FB.raw.vcf \
	-R /data/$REFERENCE --split-multi-allelics \
	-O /data/${i%%.dedup.bam}/${i%%.dedup.bam}.FB.split.vcf

docker run -v /media/array/callers_proj/raw_data:/data broadinstitute/gatk:latest gatk VariantFiltration \
	-V /data/${i%%.dedup.bam}/${i%%.dedup.bam}.FB.split.vcf -filter "RPR < 1" --filter-name "RPR1" \
	-filter "RPL < 1" --filter-name "RPL1" -filter "SAF < 1" --filter-name "SAF1" \
	-filter "SAR < 1" --filter-name "SAR1" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "QUAL / AO < 10.0" \
	--filter-name "QUAQLbyAO10" -O /data/${i%%.dedup.bam}/${i%%.dedup.bam}_FB_STANDART.vcf


