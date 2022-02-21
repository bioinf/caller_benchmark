#!/bin/bash

REFERENCE='GRCh37.primary_assembly.genome.fa'
CALLINGREGIONS='gencode.v19.revamp.bed.gz'
PREFIX='/media/array/callers_proj/raw_data'

i=$1

mkdir ${PREFIX}/${i%%.dedup.bam}/${i%%.dedup.bam}_clair3/

#wget http://www.bio8.cs.hku.hk/clair_models/illumina/12345.tar
#tar -xf 12345.tar

docker run -v ${PREFIX}:/data hkubal/clair3:latest \
	/opt/bin/run_clair3.sh --bam_fn=/data/bams/$i \
	--ref_fn=/data/$REFERENCE --bed_fn=/data/$CALLINGREGIONS \
	--threads=40 --platform="ilmn" --model_path="/opt/models/ilmn" \
	--output=/data/${i%%.dedup.bam}/${i%%.dedup.bam}_clair3/


