#!/bin/bash

for i in *.dedup.bam
do
#	sleep 2
#	while [ $( docker ps | grep 'clair' | wc -l ) -ge 10 ] ; do sleep 1 ; done
	/media/array/callers_proj/raw_data/call_clair3.sh $i
done
