#!/bin/bash

for i in *.dedup.bam
do
	sleep 2
	while [ $( ps -AF | grep 'call_hc_WGS.sh' | wc -l ) -ge 21 ] ; do sleep 1 ; done
	../../call_hc_WGS.sh $i &
done
