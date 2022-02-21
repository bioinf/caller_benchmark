#!/usr/bin/env bash

#MAIN_DIR=$1

function make_comparison()
{
  PARAMS=("$@") 
  CURR_VCF="${PARAMS[0]}"

  TMP=$(basename "${CURR_VCF}")

  SAMPLE=$(cut -d '_' -f 1 <<< "${TMP}")
  TYPE=$(cut -d '_' -f 2 <<< "${TMP}")
  CALLER=$(cut -d '_' -f 4 <<< "${TMP}")
  ALIGNER=$(cut -d '_' -f 3 <<< "${TMP}")

  echo $SAMPLE $TYPE $ALIGNER $CALLER 

  if [[ "${SAMPLE}" == "HG001" ]]; then
      ref_bed=/media/array/callers_proj/giab/master.bed
      ref_vcf=/media/array/callers_proj/giab/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz
      if [[ "${TYPE}" == "GENOME" ]]; then
        bed=/media/array/callers_proj/raw_data/gencode.v19.revamp.bed
      else
        bed=/media/array/callers_proj/bed/exome_sureselect.bed
      fi
  fi
  if [[ "${SAMPLE}" == "HG002" ]]; then
      ref_bed=/media/array/callers_proj/giab/master.bed
      ref_vcf=/media/array/callers_proj/giab/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz
      if [[ "${TYPE}" == "GENOME" ]]; then
        bed=/media/array/callers_proj/raw_data/gencode.v19.revamp.bed
      else
        bed=/media/array/callers_proj/bed/exome_sureselect.bed
      fi
  fi
  if [[ "${SAMPLE}" == "HG003" ]]; then
      ref_bed=/media/array/callers_proj/giab/master.bed
      ref_vcf=/media/array/callers_proj/giab/HG003_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v.3.3.2_highconf.vcf.gz
      if [[ "${TYPE}" == "GENOME" ]]; then
        bed=/media/array/callers_proj/raw_data/gencode.v19.revamp.bed
      else
        bed=/media/array/callers_proj/bed/exome_sureselect.bed
      fi
  fi
  if [[ "${SAMPLE}" == "HG004" ]]; then
      ref_bed=/media/array/callers_proj/giab/master.bed
      ref_vcf=/media/array/callers_proj/giab/HG004_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v.3.3.2_highconf.vcf.gz
      if [[ "${TYPE}" == "GENOME" ]]; then
        bed=/media/array/callers_proj/raw_data/gencode.v19.revamp.bed
      else
        bed=/media/array/callers_proj/bed/exome_sureselect.bed
      fi
  fi
  if [[ "${SAMPLE}" == "HG005" ]]; then
      ref_bed=/media/array/callers_proj/giab/master.bed
      ref_vcf=/media/array/callers_proj/giab/HG005_GRCh37_highconf_CG-IllFB-IllGATKHC-Ion-SOLID_CHROM1-22_v.3.3.2_highconf.vcf.gz
      if [[ "${TYPE}" == "GENOME" ]]; then
        bed=/media/array/callers_proj/raw_data/gencode.v19.revamp.bed
      else
        bed=/media/array/callers_proj/bed/exome_sureselect.bed
      fi
  fi
  if [[ "${SAMPLE}" == "HG006" ]]; then
      ref_bed=/media/array/callers_proj/giab/master.bed
      ref_vcf=/media/array/callers_proj/giab/HG006_GIAB_GRCh37_highconf_CG-IllFB-IllSNT-10X_CHROM1-22_v.3.3.2_highconf.vcf.gz
      if [[ "${TYPE}" == "GENOME" ]]; then
        bed=/media/array/callers_proj/raw_data/gencode.v19.revamp.bed
      else
        bed=/media/array/callers_proj/bed/exome_sureselect.bed
      fi
  fi
  if [[ "${SAMPLE}" == "HG007" ]]; then
      ref_bed=/media/array/callers_proj/giab/master.bed
      ref_vcf=/media/array/callers_proj/giab/HG007_GIAB_GRCh37_highconf_CG-IllFB-IllSNT-10X_CHROM1-22_v.3.3.2_highconf.vcf.gz
      if [[ "${TYPE}" == "GENOME" ]]; then
        bed=/media/array/callers_proj/raw_data/gencode.v19.revamp.bed
      else
        bed=/media/array/callers_proj/bed/exome_sureselect.bed
      fi
  fi

  HAPPY=/media/array/callers_proj/happy/bin/hap.py
  REF=/media/array/callers_proj/raw_data/GRCh37.primary_assembly.genome.fa
  REF_RTG=/media/array/callers_proj/GRCh37.primay_assembly.genome.rtg/
  RTG_SH=/media/array/callers_proj/rtg-tools-3.12/rtg
  EVAL_DIR=/media/array/callers_proj/eval_revised/${TMP%%.vcf.gz}_eval_data/
#  rm -rf $EVAL_DIR
  mkdir $EVAL_DIR

  python $HAPPY $ref_vcf $CURR_VCF -r $REF -f $ref_bed --threads 4 --engine vcfeval -T $bed -o "${EVAL_DIR}/report" --engine-vcfeval-template $REF_RTG --engine-vcfeval-path $RTG_SH --stratification /media/array/callers_proj/GRCh37-giab-stratifications/v2.0-GRCh37-stratifications_plusCustom.tsv
  echo $SAMPLE $TYPE $ALIGNER $CALLER 
  echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
}
export -f make_comparison

/usr/bin/parallel -j 14 make_comparison ::: /media/array/callers_proj/bgzipped/*GENOME*STRELKA_STANDART*.vcf.gz
wait

#for i in $( tree /media/array/callers_proj/eval_stratified | grep -B 1 report.runinf | grep -oP 'HG.*' | sed 's/_eval_data/.vcf.gz/' | perl -pe 's|^|/media/array/callers_proj/bgzipped/|' ) ; do make_comparison $i ; done


