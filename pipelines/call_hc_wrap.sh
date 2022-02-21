#!/bin/bash

for i in *.dedup.bam
do
	sleep 2
	while [ $( ps -AF | grep 'call_hc.sh' | wc -l ) -ge 19 ] ; do sleep 1 ; done
	../call_hc.sh $i &
done
