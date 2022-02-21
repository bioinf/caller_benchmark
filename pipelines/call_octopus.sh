#!/bin/bash

REFERENCE='GRCh37.primary_assembly.genome.fa'

i=$1

docker run -v /media/array/callers_proj/raw_data/:/data octopus octopus -R /data/$REFERENCE -I /data/bams/WGS/$i -t /data/gencode.v19.revamp.octopus.list -o /data/${i%%.dedup.bam}/${i%%.dedup.bam}_OCTOPUS_STANDART.vcf --threads 4

docker run -v /media/array/callers_proj/raw_data/:/data octopus octopus -R /data/$REFERENCE -I /data/bams/WGS/$i --filter-vcf /data/${i%%.dedup.bam}/${i%%.dedup.bam}_OCTOPUS_STANDART.vcf --forest /opt/octopus/resources/forests/germline.v0.7.4.forest.gz -o /data/${i%%.dedup.bam}/${i%%.dedup.bam}_OCTOPUS_FOREST.vcf --threads 4
