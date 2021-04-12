#!/bin/bash

cat $( find * | grep extended | head -n1 ) | head -n1 | perl -pe 's|^|Sample,ExperimentType,Aligner,CallerFilter,|' | perl -pe 's|,|\t|g'

for i in *eval_data
do
	SAMPLE=${i%%_*}
	OME=$( echo $i | grep -oP '[A-Z]+OME' )
	ALIGNER=$( echo $i | grep -oP 'NOVOALIGN|ISAAC4|BOWTIE2|BWA' )
	CALLER_FILTER=$(echo $i | grep -oP '(DV|HC|FB|ST)_[A-Z0-9]+' )
#	echo $SAMPLE $OME $ALIGNER $CALLER_FILTER
	tail -n +2 ${i}/report.extended.csv | perl -pe 's|^|'"$SAMPLE,$OME,$ALIGNER,$CALLER_FILTER"',|' | perl -pe 's|,|\t|g' 
done
