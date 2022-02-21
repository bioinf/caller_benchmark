#!/bin/bash

for i in $( ls *.gz | grep -oP 'HG00\d_[EXGEN]+OME_BWA' | sort -u ) ; do for j in ${i}*gz ; do for k in ${i}*gz ; do ./eval_concord.py $j $k ; done ; done ; done
