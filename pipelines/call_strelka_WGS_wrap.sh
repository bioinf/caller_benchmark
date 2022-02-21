#!/bin/bash

#sleep 90

for i in *.dedup.bam
do
	sleep 2
	while [ $( ps -AF | grep 'call_strelka.sh' | wc -l ) -ge 13 ] ; do sleep 1 ; done
	/media/array/callers_proj/raw_data/call_strelka_WGS.sh $i &
done
