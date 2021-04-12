#!/usr/bin/env bash

#MAIN_DIR=$1

function make_comparison()
{
  PARAMS=("$@") 
  CURR_VCF="${PARAMS[0]}"

  TMP=$(basename "${CURR_VCF}")

  SAMPLE=$(cut -d '_' -f 1 <<< "${TMP}")
  TYPE=$(cut -d '_' -f 2 <<< "${TMP}")
  CALLER=$(cut -d '_' -f 5 <<< "${TMP}")
  ALIGNER=$(cut -d '_' -f 4 <<< "${TMP}")

  echo $SAMPLE $TYPE $ALIGNER $CALLER 

  if [[ "${SAMPLE}" == "HG001" ]]; then
      ref_bed=/media/array/callers_proj/giab/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed
      ref_vcf=/media/array/callers_proj/giab/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz
      if [[ "${TYPE}" == "GENOME" ]]; then
        bed=/media/array/callers_proj/bed/gencode.r33.hg19.transcripts_padded_150.bed
      else
        bed=/media/array/callers_proj/bed/HG001_exome_regions_padded_150.bed
      fi
  fi
  if [[ "${SAMPLE}" == "HG002" ]]; then
      ref_bed=/media/array/callers_proj/giab/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed
      ref_vcf=/media/array/callers_proj/giab/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz
      if [[ "${TYPE}" == "GENOME" ]]; then
        bed=/media/array/callers_proj/bed/gencode.r33.hg19.transcripts_padded_150.bed
      else
        bed=/media/array/callers_proj/bed/HG002_exome_regions_padded_150.bed
      fi
  fi
  if [[ "${SAMPLE}" == "HG003" ]]; then
      ref_bed=/media/array/callers_proj/giab/HG003_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed
      ref_vcf=/media/array/callers_proj/giab/HG003_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v.3.3.2_highconf.vcf.gz
      if [[ "${TYPE}" == "GENOME" ]]; then
        bed=/media/array/callers_proj/bed/gencode.r33.hg19.transcripts_padded_150.bed
      else
        bed=/media/array/callers_proj/bed/HG003_exome_regions_padded_150.bed
      fi
  fi
  if [[ "${SAMPLE}" == "HG004" ]]; then
      ref_bed=/media/array/callers_proj/giab/HG004_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed
      ref_vcf=/media/array/callers_proj/giab/HG004_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v.3.3.2_highconf.vcf.gz
      if [[ "${TYPE}" == "GENOME" ]]; then
        bed=/media/array/callers_proj/bed/gencode.r33.hg19.transcripts_padded_150.bed
      else
        bed=/media/array/callers_proj/bed/HG004_exome_regions_padded_150.bed
      fi
  fi
  if [[ "${SAMPLE}" == "HG005" ]]; then
      ref_bed=/media/array/callers_proj/giab/HG005_GRCh37_highconf_CG-IllFB-IllGATKHC-Ion-SOLID_CHROM1-22_v.3.3.2_highconf_noMetaSV.bed
      ref_vcf=/media/array/callers_proj/giab/HG005_GRCh37_highconf_CG-IllFB-IllGATKHC-Ion-SOLID_CHROM1-22_v.3.3.2_highconf.vcf.gz
      if [[ "${TYPE}" == "GENOME" ]]; then
        bed=/media/array/callers_proj/bed/gencode.r33.hg19.transcripts_padded_150.bed
      else
        bed=/media/array/callers_proj/bed/HG005_exome_regions_padded_150.bed
      fi
  fi

  HAPPY=/media/array/callers_proj/happy/bin/hap.py
  REF=/media/array/callers_proj/GRCh37.p13.genome.fa
  REF_RTG=/media/array/callers_proj/GRCh37.p13.genome.rtg/
  RTG_SH=/media/array/callers_proj/rtg-tools-3.12/rtg
  EVAL_DIR=/media/array/callers_proj/eval_stratified/${TMP%%.vcf.gz}_eval_data/
#  rm -rf $EVAL_DIR
  mkdir $EVAL_DIR

  python $HAPPY $ref_vcf $CURR_VCF -r $REF -f $ref_bed --threads 4 --engine vcfeval -T $bed -o "${EVAL_DIR}/report" --engine-vcfeval-template $REF_RTG --engine-vcfeval-path $RTG_SH --stratification /media/array/callers_proj/GRCh37-giab-stratifications/v2.0-GRCh37-stratifications_plusCustom.tsv
  echo $SAMPLE $TYPE $ALIGNER $CALLER 
  echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
}
export -f make_comparison

/usr/bin/parallel -j 10 make_comparison ::: /media/array/callers_proj/bgzipped/*.vcf.gz
wait

for i in $( tree /media/array/callers_proj/eval_stratified | grep -B 1 report.runinf | grep -oP 'HG.*' | sed 's/_eval_data/.vcf.gz/' | perl -pe 's|^|/media/array/callers_proj/bgzipped/|' ) ; do make_comparison $i ; done


