#!/bin/bash

#sleep 18000

for i in *.dedup.bam
do
	sleep 2
	while [ $( ps -AF | grep 'call_fb.sh' | wc -l ) -ge 19 ] ; do sleep 1 ; done
	/media/array/callers_proj/raw_data/call_fb.sh $i &
done
