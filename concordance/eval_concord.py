#!/usr/bin/env python

import sys
import gzip
import re

first_vcf, second_vcf = sys.argv[1:]
case = re.findall('HG00\d_[EXGEN]+OME', first_vcf)[0]
first = re.findall('((DV|CLAIR3|STRELKA|OCTOPUS|FB|HC)_[A-Z0-9]+)', first_vcf)[0][0]
second = re.findall('((DV|CLAIR3|STRELKA|OCTOPUS|FB|HC)_[A-Z0-9]+)', second_vcf)[0][0]

variants = {}
with gzip.open(second_vcf, 'rt') as vcf_file:
    for line in vcf_file:
        if line.startswith('#'):
            continue
        content = line.strip().split('\t')
        if content[6] != 'PASS':
            continue
        variants[':'.join(content[:2])] = None

total_cnt = 0
uq_discordant = 0
with gzip.open(first_vcf, 'rt') as vcf_file:
    for line in vcf_file:
        if line.startswith('#'):
            continue
        content = line.strip().split('\t')
        if content[6] != 'PASS':
            continue
        total_cnt += 1
        if ':'.join(content[:2]) not in variants:
            uq_discordant += 1

print(f'{case}\t{first}\t{second}\t{total_cnt}\t{uq_discordant}\t{uq_discordant/total_cnt}')
