#!/bin/bash

REFERENCE='GRCh37.primary_assembly.genome.fa'
CALLINGREGIONS='gencode.v19.revamp.bed.gz'

i=$1

# Manta+Strelka
mkdir /media/array/callers_proj/raw_data/${i%%.dedup.bam}/{mantawd,strelkawd}
configManta.py --bam $i --referenceFasta /media/array/callers_proj/raw_data/$REFERENCE --runDir /media/array/callers_proj/raw_data/${i%%.dedup.bam}/mantawd --callRegions /media/array/callers_proj/raw_data/$CALLINGREGIONS --exome
/media/array/callers_proj/raw_data/${i%%.dedup.bam}/mantawd/runWorkflow.py -j 4

configureStrelkaGermlineWorkflow.py --bam $i --referenceFasta /media/array/callers_proj/raw_data/$REFERENCE \
	--runDir /media/array/callers_proj/raw_data/${i%%.dedup.bam}/strelkawd \
	--callRegions /media/array/callers_proj/raw_data/$CALLINGREGIONS \
	--indelCandidates /media/array/callers_proj/raw_data/${i%%.dedup.bam}/mantawd/results/variants/candidateSmallIndels.vcf.gz \
	--exome
/media/array/callers_proj/raw_data/${i%%.dedup.bam}/strelkawd/runWorkflow.py -m local -j 4

cp /media/array/callers_proj/raw_data/${i%%.dedup.bam}/strelkawd/results/variants/variants.vcf.gz /media/array/callers_proj/raw_data/${i%%.dedup.bam}/${i%%.dedup.bam}_STRELKA_STANDART.vcf.gz
