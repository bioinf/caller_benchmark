#!/bin/bash

for i in *.dedup.bam
do
	sleep 1
	while [ $( ps -AF | grep call_octopus | wc -l ) -ge 11 ] ; do sleep 1 ; done
	/media/array/callers_proj/raw_data/call_octopus.sh $i &
done
