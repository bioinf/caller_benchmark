#!/bin/bash

for i in *.vcf.gz
do 
	smp=$( echo $i | grep -oP 'HG00\d' )
	/media/array/callers_proj/rtg-tools-3.12/rtg vcfeval -b /media/array/callers_proj/concordance_analysis/${smp}*gz --bed-regions /media/array/callers_proj/concordance_analysis/GRCh37_WES_CDS.bed -c $i -e /media/array/callers_proj/giab/giab42/master.bed -o ${i%%.vcf.gz}.rtg/ -m combine -t /media/array/callers_proj/GRCh37.primay_assembly.genome.rtg/
done
