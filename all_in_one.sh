#!/usr/bin/env bash

eval "$(conda shell.bash hook)"

REF=/home/lmbs02/bio/databases/referenses/hg19_37/gencode/GRCh37_p13_33_release/GRCh37.primary_assembly.genome.fa
BOWTIE2_REF=/home/lmbs02/bio/databases/referenses/hg19_37/gencode/GRCh37_p13_33_release/GRCh37.primary_assembly.genome
NOVOALIGN_REF=/home/lmbs02/bio/databases/referenses/hg19_37/gencode/GRCh37_p13_33_release/GRCh37.primary_assembly.genome
MILS_INDELS=/home/lmbs02/bio/databases/gatk_bundle/hg19/hg19_v0_Mills_and_1000G_gold_standard.indels.b37.vcf.gz
GENOMES=/home/lmbs02/bio/databases/gatk_bundle/hg19/hg19_v0_Homo_sapiens_assembly19.known_indels_20120518.vcf.gz
DBSNP=/home/lmbs02/bio/databases/gatk_bundle/hg19/hg19_v0_Homo_sapiens_assembly19.dbsnp138.vcf.gz
GENOME_FILE=/home/lmbs02/bio/databases/referenses/hg19_37/gencode/GRCh37_p13_32_release/for_bedtools.txt
REF_DICT=/home/lmbs02/bio/databases/referenses/hg19_37/gencode/GRCh37_p13_33_release/GRCh37.primary_assembly.genome.dict

BWA=/home/lmbs02/bio/biosoft/bwa/bwa-0.7.17/bwa
BOWTIE2=/home/lmbs02/bio/biosoft/bowtie2/bowtie2-2.4.1-linux-x86_64/bowtie2
NOVOALIGN=/home/lmbs02/bio/biosoft/novocraft/novocraft/novoalign
ISAAC4=/home/lmbs02/bio/biosoft/isaac4/Isaac4-Isaac-04.18.11.09/isaac-bild/bin/isaac-align
SAMTOOLS=/home/lmbs02/bio/biosoft/samtools/samtools-1.10/samtools
BCFTOOLS=/home/lmbs02/bio/biosoft/bcftools/bcftools-1.10.2/bcftools
PICARD=/home/lmbs02/bio/biosoft/picard/picard-2.22.0/picard.jar
GATK4=/home/lmbs02/bio/biosoft/gatk/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar
BEDTOOLS=/home/lmbs02/bio/biosoft/bedtools/bedtools2-2.29.0/bin/bedtools
DEEPVARIANT=/home/lmbs02/bio/biosoft/deepvariant/deepvariant/binaries/DeepVariant/0.10.0/DeepVariant-0.10.0
DEEPVARIANT_WES_MODEL=/home/lmbs02/bio/biosoft/deepvariant/deepvariant/models/DeepVariant/0.10.0/DeepVariant-inception_v3-0.10.0+data-wes_standard/model.ckpt
DEEPVARIANT_WGS_MODEL=/home/lmbs02/bio/biosoft/deepvariant/deepvariant/models/DeepVariant/0.10.0/DeepVariant-inception_v3-0.10.0+data-wes_standard/model.ckpt
STRELKA=/home/lmbs02/bio/biosoft/strelka/strelka-2.9.10.centos6_x86_64/bin/
MANTA=/home/lmbs02/bio/biosoft/manta/manta-1.6.0.centos6_x86_64/bin/
FREEBAYES=/home/lmbs02/bio/biosoft/freebayes/freebayes-1.3.1/bin/freebayes

function align_with_bwa()
{
  PARAMS=("$@") 
  CURR_DIR="${PARAMS[0]}"
  BWA="${PARAMS[1]}"
  SAMTOOLS="${PARAMS[2]}"
  REF="${PARAMS[3]}"

  SAMPLE=$(basename "${CURR_DIR}")
  ANALYSIS_DIR=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/analysis/
  ANALYSIS_SAMPLE_DIR="${ANALYSIS_DIR}/${SAMPLE}"

  mkdir -p "${ANALYSIS_SAMPLE_DIR}"
  sleep 1
  
  for FASTQ in $(find "${CURR_DIR}" -type f -maxdepth 1 -name *R1*.fastq.gz); do
    tmp=$(basename "${FASTQ}")
    lane="${tmp/_R1_001.fastq.gz/}"
    R1="${FASTQ}"
    R2="${R1/R1/R2}"  
    LANE_BAM="${ANALYSIS_SAMPLE_DIR}/${lane}.bam"
    "${BWA}" mem -R "@RG\tID:${lane}\tPL:ILLUMINA\tLB:LIBRARY\tSM:${SAMPLE}" "${REF}" -K 100000000 -M "$R1" "$R2" -t 24 | "${SAMTOOLS}" view -bT "${REF}" -> "${LANE_BAM}"
  done
}
export -f align_with_bwa

function align_with_bowtie2()
{
  PARAMS=("$@") 
  CURR_DIR="${PARAMS[0]}"
  BOWTIE2="${PARAMS[1]}"
  SAMTOOLS="${PARAMS[2]}"
  BOWTIE2_REF="${PARAMS[3]}"
  REF="${PARAMS[4]}"

  SAMPLE=$(basename "${CURR_DIR}")
  ANALYSIS_DIR=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/analysis/
  ANALYSIS_SAMPLE_DIR="${ANALYSIS_DIR}/${SAMPLE}"

  mkdir -p "${ANALYSIS_SAMPLE_DIR}"
  sleep 1
  
  for FASTQ in $(find "${CURR_DIR}" -type f -maxdepth 1 -name *R1*.fastq.gz); do
    tmp=$(basename "${FASTQ}")
    lane="${tmp/_R1_001.fastq.gz/}"
    R1="${FASTQ}"
    R2="${R1/R1/R2}"  
    LANE_BAM="${ANALYSIS_SAMPLE_DIR}/${lane}.bam"
    "${BOWTIE2}" -x "${BOWTIE2_REF}" -1 "$R1" -2 "$R2" -p 16 --rg-id ${lane} --rg SM:${SAMPLE} --rg LB:LIBRARY --rg PL:ILLUMINA  | "${SAMTOOLS}" view -bT "${REF}" -> "${LANE_BAM}"
  done
}
export -f align_with_bowtie2

function align_with_isaac4()
{
  PARAMS=("$@") 
  CURR_DIR="${PARAMS[0]}"
  ISAAC4="${PARAMS[1]}"
  SAMTOOLS="${PARAMS[2]}"
  REF="${PARAMS[3]}"

  SAMPLE=$(basename "${CURR_DIR}")
  ANALYSIS_DIR=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/analysis/
  ANALYSIS_SAMPLE_DIR="${ANALYSIS_DIR}/${SAMPLE}"

  mkdir -p "${ANALYSIS_SAMPLE_DIR}"
 
  for FASTQ in $(find "${CURR_DIR}" -type f -maxdepth 1 -name *R1*.fastq.gz); do
    name=$(basename $FASTQ)
    tmp="${name%_R1_001.fastq.gz}"
    lane="${tmp##*_}"
    R1="${FASTQ}"
    R2="${R1/R1/R2}"
    numb="${lane#L00}"
    ln -f $R1 "${CURR_DIR}/lane${numb}_read1.fastq.gz" 
    ln -f $R2 "${CURR_DIR}/lane${numb}_read2.fastq.gz"
  done

  "${ISAAC4}" -r $REF -b $CURR_DIR -f fastq-gz --output-directory $ANALYSIS_SAMPLE_DIR -m 55 --mark-duplicates 0  --verbosity 4 --variable-read-length 1
  mv "${ANALYSIS_SAMPLE_DIR}/Projects/default/default/sorted.bam" "${ANALYSIS_SAMPLE_DIR}/${SAMPLE}_L001.bam"
  rm -rf "${ANALYSIS_SAMPLE_DIR}/Projects/"
  rm -rf "${ANALYSIS_SAMPLE_DIR}/Reports/"
  rm -rf "${ANALYSIS_SAMPLE_DIR}/Stats/"
}
export -f align_with_isaac4

function align_with_novoalign()
{
  PARAMS=("$@") 
  CURR_DIR="${PARAMS[0]}"
  NOVOALIGN="${PARAMS[1]}"
  SAMTOOLS="${PARAMS[2]}"
  NOVOALIGN_REF="${PARAMS[3]}"
  REF="${PARAMS[4]}"

  SAMPLE=$(basename "${CURR_DIR}")
  ANALYSIS_DIR=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/analysis/
  ANALYSIS_SAMPLE_DIR="${ANALYSIS_DIR}/${SAMPLE}"

  mkdir -p "${ANALYSIS_SAMPLE_DIR}"
  sleep 1
  
  for FASTQ in $(find "${CURR_DIR}" -type f -maxdepth 1 -name *R1*.fastq.gz); do
    tmp=$(basename "${FASTQ}")
    lane="${tmp/_R1_001.fastq.gz/}"
    R1="${FASTQ}"
    R2="${R1/R1/R2}"  
    LANE_BAM="${ANALYSIS_SAMPLE_DIR}/${lane}.bam"
    "${NOVOALIGN}" -d "${NOVOALIGN_REF}" -f "$R1" "$R2" -o SAM  "@RG\tID:${lane}\tSM:${SAMPLE}\tPL:illumina" 2> /dev/null | "${SAMTOOLS}" view -bT "${REF}"  -> "${LANE_BAM}"
  done

}
export -f align_with_novoalign

function merge_lane_bams()
{
  PARAMS=("$@") 
  CURR_DIR="${PARAMS[0]}"
  SAMTOOLS="${PARAMS[1]}"

  SAMPLE=$(basename "${CURR_DIR}")
  ANALYSIS_DIR=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/analysis/
  ANALYSIS_SAMPLE_DIR="${ANALYSIS_DIR}/${SAMPLE}"
  RAW_BAM="${ANALYSIS_SAMPLE_DIR}/${SAMPLE}_raw.bam"
  BAMS_COUNT=$(find "${ANALYSIS_SAMPLE_DIR}" -maxdepth 1 -type f -name '*.bam' -printf x | wc -c)

  if [ "$BAMS_COUNT" -gt 1 ];
  then

    for LANE_BAM in $(find "${ANALYSIS_SAMPLE_DIR}" -maxdepth 1 -type f -name *.bam); do
      input_str+=" "
      input_str+="${LANE_BAM}"
    done

    eval "${SAMTOOLS} merge -f -l 9 -@ 16 ${RAW_BAM} ${input_str}"

  else

    LANE_BAM=$(find "${ANALYSIS_SAMPLE_DIR}" -maxdepth 1 -type f -name '*.bam')
    mv "${LANE_BAM}" "${RAW_BAM}"

  fi
}
export -f merge_lane_bams

function sort_after_align()
{
  PARAMS=("$@") 
  CURR_DIR="${PARAMS[0]}"
  SAMTOOLS="${PARAMS[1]}"

  SAMPLE=$(basename "${CURR_DIR}")
  ANALYSIS_DIR=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/analysis/
  ANALYSIS_SAMPLE_DIR="${ANALYSIS_DIR}/${SAMPLE}"
  RAW_BAM="${ANALYSIS_SAMPLE_DIR}/${SAMPLE}_raw.bam"
  SORTED_BAM="${ANALYSIS_SAMPLE_DIR}/${SAMPLE}_sorted.bam"

  "${SAMTOOLS}" sort -o "${SORTED_BAM}" -l 9 -@ 16 "${RAW_BAM}"
  "${SAMTOOLS}" index -@ 16 "${SORTED_BAM}"

  for LANE_BAM in $(find "${ANALYSIS_SAMPLE_DIR}" -type f -maxdepth 1 -not -name *sorted*); do
  	echo ${LANE_BAM}
    rm ${LANE_BAM}
  done
}
export -f sort_after_align

function mark_dublicates()
{
  PARAMS=("$@") 
  CURR_DIR="${PARAMS[0]}"
  PICARD="${PARAMS[1]}"
  SAMTOOLS="${PARAMS[2]}"

  SAMPLE=$(basename "${CURR_DIR}")
  ANALYSIS_DIR=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/analysis/
  ANALYSIS_SAMPLE_DIR="${ANALYSIS_DIR}/${SAMPLE}"
  SORTED_BAM="${ANALYSIS_SAMPLE_DIR}/${SAMPLE}_sorted.bam"
  SORTED_BAI="${ANALYSIS_SAMPLE_DIR}/${SAMPLE}_sorted.bam.bai"
  DUB_MARKED_BAM=${ANALYSIS_SAMPLE_DIR}/${SAMPLE}_dub_marked.bam
  DUB_MARK_METRICS=${ANALYSIS_SAMPLE_DIR}/${SAMPLE}_dub_mark_mertics.txt

  java -jar -Xmx11g "${PICARD}" MarkDuplicates I="${SORTED_BAM}" O="${DUB_MARKED_BAM}" M="${DUB_MARK_METRICS}" SORTING_COLLECTION_SIZE_RATIO=0.1
  "${SAMTOOLS}" index -@ 16 "${DUB_MARKED_BAM}"

  rm "${SORTED_BAM}" "${SORTED_BAI}"
}
export -f mark_dublicates

function base_recalibration()
{
  PARAMS=("$@") 
  CURR_DIR="${PARAMS[0]}"
  GATK4="${PARAMS[1]}"
  REF="${PARAMS[2]}"
  MILS_INDELS="${PARAMS[3]}"
  GENOMES="${PARAMS[4]}"
  DBSNP="${PARAMS[5]}"

  SAMPLE=$(basename "${CURR_DIR}")
  ANALYSIS_DIR=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/analysis/
  ANALYSIS_SAMPLE_DIR="${ANALYSIS_DIR}/${SAMPLE}"
  DUB_MARKED_BAM=${ANALYSIS_SAMPLE_DIR}/${SAMPLE}_dub_marked.bam
  RECALIB_REPORT=${ANALYSIS_SAMPLE_DIR}/${SAMPLE}_recal_data.table
  PADDED_BED=${CURR_DIR}/gencode.r33.hg19.transcripts_padded_150.bed

  java -jar -jar -Xmx11g "${GATK4}" BaseRecalibrator -R "${REF}" -I "${DUB_MARKED_BAM}" --known-sites "${DBSNP}" --known-sites "${MILS_INDELS}" --known-sites "${GENOMES}" -O "${RECALIB_REPORT}" -L "${PADDED_BED}"
}
export -f base_recalibration

function apply_recalibration()
{
  PARAMS=("$@") 
  CURR_DIR="${PARAMS[0]}"
  GATK4="${PARAMS[1]}"
  REF="${PARAMS[2]}"

  SAMPLE=$(basename "${CURR_DIR}")
  ANALYSIS_DIR=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/analysis/
  ANALYSIS_SAMPLE_DIR="${ANALYSIS_DIR}/${SAMPLE}"
  DUB_MARKED_BAM=${ANALYSIS_SAMPLE_DIR}/${SAMPLE}_dub_marked.bam
  RECALIB_REPORT=${ANALYSIS_SAMPLE_DIR}/${SAMPLE}_recal_data.table
  PADDED_BED=${CURR_DIR}/gencode.r33.hg19.transcripts_padded_150.bed
  RECALIBRATED_BAM=${ANALYSIS_SAMPLE_DIR}/${SAMPLE}_recalibrated.bam

  java -jar -jar -Xmx11g "${GATK4}" ApplyBQSR -R "${REF}" -I "${DUB_MARKED_BAM}" -O "${RECALIBRATED_BAM}" --bqsr-recal-file "${RECALIB_REPORT}"
  rm "${RECALIB_REPORT}"
}
export -f apply_recalibration

function gtk4_hc()
{
  PARAMS=("$@") 
  CURR_DIR="${PARAMS[0]}"
  GATK4="${PARAMS[1]}"
  REF="${PARAMS[2]}"

  SAMPLE=$(basename "${CURR_DIR}")
  ANALYSIS_DIR=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/analysis/
  ANALYSIS_SAMPLE_DIR="${ANALYSIS_DIR}/${SAMPLE}"
  PADDED_BED=${CURR_DIR}/gencode.r33.hg19.transcripts_padded_150.bed
  RECALIBRATED_BAM=${ANALYSIS_SAMPLE_DIR}/${SAMPLE}_recalibrated.bam
  BAMOUT_BAM=${ANALYSIS_SAMPLE_DIR}/${SAMPLE}_bamout.bam
  HC_VCF=${ANALYSIS_SAMPLE_DIR}/${SAMPLE}_GENOME_HG19_BWA_HC_HARDFILTER.vcf
  HC_GVCF=${ANALYSIS_SAMPLE_DIR}/${SAMPLE}_HC.g.vcf.gz

  echo "Directory : ${CURR_DIR}"

  java -jar "${GATK4}" HaplotypeCaller -R "${REF}" -I "${RECALIBRATED_BAM}" -bamout "${BAMOUT_BAM}" -O "${HC_VCF}" -L "${PADDED_BED}"
  java -jar "${GATK4}" CNNScoreVariants -R "${REF}" -V "${HC_VCF}" -O "${HC_VCF}.tmp.vcf" -tensor-type reference
  mv "${HC_VCF}.tmp.vcf" "${HC_VCF}"
  java -jar "${GATK4}" CNNScoreVariants -I "${BAMOUT_BAM}" -R "${REF}" -V "${HC_VCF}" -O "${HC_VCF}.tmp.vcf" -tensor-type read_tensor
  mv "${HC_VCF}.tmp.vcf" "${HC_VCF}"
  java -jar "${GATK4}" FilterVariantTranches -V "${HC_VCF}" --output "${HC_VCF}.tmp.vcf" \
  --info-key CNN_1D --snp-tranche 90 --indel-tranche 90 \
  --resource /home/lmbs02/bio/databases/gatk_bundle/hg19/hg19_v0_Mills_and_1000G_gold_standard.indels.b37.vcf.gz \
  --resource /home/lmbs02/bio/databases/gatk_bundle/hg19/hg19_v0_1000G_phase1.snps.high_confidence.b37.vcf.gz \
  --resource /home/lmbs02/bio/databases/gatk_bundle/hg19/hg19_v0_hapmap_3.3.b37.vcf.gz \
  mv "${HC_VCF}.tmp.vcf" "${HC_VCF}"
  java -jar "${GATK4}" FilterVariantTranches -V "${HC_VCF}" --output "${HC_VCF}.tmp.vcf" \
  --resource /home/lmbs02/bio/databases/gatk_bundle/hg19/hg19_v0_1000G_phase1.snps.high_confidence.b37.vcf.gz \
  --resource /home/lmbs02/bio/databases/gatk_bundle/hg19/hg19_v0_hapmap_3.3.b37.vcf.gz \
  --resource /home/lmbs02/bio/databases/gatk_bundle/hg19/hg19_v0_Mills_and_1000G_gold_standard.indels.b37.vcf.gz \
  --info-key CNN_2D --snp-tranche 99.9 --indel-tranche 99.5
  mv "${HC_VCF}.tmp.vcf" "${HC_VCF}"
  java -jar "${GATK4}" SelectVariants -V "${HC_VCF}" -select-type SNP -select-type MIXED -O "${HC_VCF}.snps.vcf"
  java -jar "${GATK4}" VariantFiltration -V "${HC_VCF}.snps.vcf" -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" \
  -filter "SOR > 3.0" --filter-name "SOR3" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" \
  --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O "${HC_VCF}.snps.tmp.vcf"
  mv "${HC_VCF}.snps.tmp.vcf" "${HC_VCF}.snps.vcf"
  java -jar "${GATK4}" SelectVariants -V "${HC_VCF}" -select-type INDEL -select-type MIXED -O "${HC_VCF}.indels.vcf"
  java -jar "${GATK4}" VariantFiltration -V "${HC_VCF}.indels.vcf" -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" \
  -filter "FS > 200.0" --filter-name "FS200" -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" -O "${HC_VCF}.indels.tmp.vcf"
  mv "${HC_VCF}.indels.tmp.vcf" "${HC_VCF}.indels.vcf"
  java -jar "${GATK4}" MergeVcfs -I "${HC_VCF}.snps.vcf" -I "${HC_VCF}.indels.vcf" -O "${HC_VCF}.tmp.vcf"
  mv "${HC_VCF}.tmp.vcf" "${HC_VCF}"
  rm "${HC_VCF}.snps.vcf" "${HC_VCF}.indels.vcf"
  java -jar "${GATK4}" CalculateGenotypePosteriors -V "${HC_VCF}" -supporting /home/lmbs02/bio/databases/gatk_bundle/hg19/b37_1000G_phase3_v4_20130502.sites.vcf.gz \
  -O "${HC_VCF}.tmp.vcf"
  mv "${HC_VCF}.tmp.vcf" "${HC_VCF}"
  java -jar "${GATK4}" VariantFiltration -V "${HC_VCF}" --genotype-filter-expression "GQ < 20" --genotype-filter-name "GQ20" -O "${HC_VCF}.tmp.vcf"
  mv "${HC_VCF}.tmp.vcf" "${HC_VCF}"
  bcftools filter -s GQ20 -e "FORMAT/FT[*]!=''" "${HC_VCF}" > "${HC_VCF}.tmp.vcf"
  mv "${HC_VCF}.tmp.vcf" "${HC_VCF}"

}
export -f gtk4_hc

function strelka_call()
{
  PARAMS=("$@") 
  CURR_DIR="${PARAMS[0]}"
  STRELKA="${PARAMS[1]}"
  MANTA="${PARAMS[2]}"
  REF="${PARAMS[3]}"

  SAMPLE=$(basename "${CURR_DIR}")
  ANALYSIS_DIR=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/analysis/
  ANALYSIS_SAMPLE_DIR="${ANALYSIS_DIR}/${SAMPLE}"
  PADDED_BED=${CURR_DIR}/gencode.r33.hg19.transcripts_padded_150.bed
  DUB_MARKED_BAM=${ANALYSIS_SAMPLE_DIR}/${SAMPLE}_dub_marked.bam
  ST2_VCF=${ANALYSIS_SAMPLE_DIR}/${SAMPLE}_GENOME_HG19_BWA_ST_STANDART.vcf
  STRELKA_DIR=${ANALYSIS_SAMPLE_DIR}/STRELKA/
  MANTA_DIR=${ANALYSIS_SAMPLE_DIR}/MANTA/
  CANDIDATE_INDELS="${MANTA_DIR}/results/variants/candidateSmallIndels.vcf.gz"

  rm -rf "${MANTA_DIR}"
  "${MANTA}configManta.py" --bam "${DUB_MARKED_BAM}" --referenceFasta "${REF}" --runDir "${MANTA_DIR}" --callRegions "${PADDED_BED}.gz" #--exome for exomes
  "${MANTA_DIR}runWorkflow.py" -j 12 
  rm -rf "${STRELKA_DIR}"
  "${STRELKA}/configureStrelkaGermlineWorkflow.py" --bam "${DUB_MARKED_BAM}" --referenceFasta "${REF}" --runDir "${STRELKA_DIR}" \
  --callRegions "${PADDED_BED}.gz" --indelCandidates "${CANDIDATE_INDELS}" #--exome for exomes
  "${STRELKA_DIR}/runWorkflow.py" -m local -j 12
  mv -f "${STRELKA_DIR}/results/variants/variants.vcf.gz" "${ST2_VCF}.gz"
  gzip -dkf "${ST2_VCF}.gz"
  rm "${ST2_VCF}.gz"

  rm -rf "${MANTA_DIR}"
  rm -rf "${STRELKA_DIR}"

  cp "${ST2_VCF}" /home/lmbs02/bio/work/process_panels/NIST_EXOMES/VCF_TO_URA/
}
export -f strelka_call

function deepvariant_call()
{
  PARAMS=("$@") 
  CURR_DIR="${PARAMS[0]}"
  DEEPVARIANT="${PARAMS[1]}"
  REF="${PARAMS[2]}"
  MODEL="${PARAMS[3]}"

  SAMPLE=$(basename "${CURR_DIR}")
  ANALYSIS_DIR=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/analysis/
  ANALYSIS_SAMPLE_DIR="${ANALYSIS_DIR}/${SAMPLE}"
  DV_VCF=${ANALYSIS_SAMPLE_DIR}/${SAMPLE}_GENOME_HG19_BWA_DV_STANDART.vcf
  DV_GVCF=${ANALYSIS_SAMPLE_DIR}/${SAMPLE}_DV.g.vcf.gz
  DUB_MARKED_BAM=${ANALYSIS_SAMPLE_DIR}/${SAMPLE}_dub_marked.bam
  PADDED_BED=${CURR_DIR}/tmp.bed #${CURR_DIR}/gencode.r33.hg19.transcripts_padded_150.bed
  MAKE_EXAMPLES="${DEEPVARIANT}/make_examples.zip"
  CALL_VARIANTS="${DEEPVARIANT}/call_variants.zip"
  POSTPROCESS_VARIANTS="${DEEPVARIANT}/postprocess_variants.zip"

  docker run --rm -v /home/lmbs02/:/home/lmbs02/ google/deepvariant:0.10.0 /opt/deepvariant/bin/run_deepvariant --model_type=WGS \
  --ref="${REF}" --reads="${DUB_MARKED_BAM}" --regions "${PADDED_BED}" --output_vcf="${DV_VCF}" --num_shards=24  #--model_type=WES for exomes

}
export -f deepvariant_call

function freebayes_call()
{
  PARAMS=("$@") 
  CURR_DIR="${PARAMS[0]}"
  FREEBAYES="${PARAMS[1]}"
  GATK4="${PARAMS[2]}"
  REF="${PARAMS[3]}"
  SAMPLE=$(basename "${CURR_DIR}")
  ANALYSIS_DIR=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/analysis/
  ANALYSIS_SAMPLE_DIR="${ANALYSIS_DIR}/${SAMPLE}"
  PADDED_BED=${CURR_DIR}/regions_padded_150.bed #gencode.r33.hg19.transcripts_padded_150.bed
  FB_VCF=${ANALYSIS_SAMPLE_DIR}/${SAMPLE}_EXOME_HG19_BOWTIE2_FB_STANDART.vcf
  DUB_MARKED_BAM=${ANALYSIS_SAMPLE_DIR}/${SAMPLE}_dub_marked.bam

  "${FREEBAYES}" --standard-filters -t "${PADDED_BED}" -f "${REF}" "${DUB_MARKED_BAM}" > "${FB_VCF}"
  java -jar "${GATK4}" VariantFiltration -V "${FB_VCF}" -filter "RPR < 1" --filter-name "RPR1" -filter "RPL < 1" --filter-name "RPL1" \
  -filter "SAF < 1" --filter-name "SAF1" -filter "SAR < 1" --filter-name "SAR1" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "QUAL / AO < 10.0" \
  --filter-name "QUAQLbyAO10" -O "${FB_VCF}.tmp.vcf"
  mv "${FB_VCF}.tmp.vcf" "${FB_VCF}"
}
export -f freebayes_call

function make_comparsion()
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
      ref_bed=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/refs/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed
      ref_vcf=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/refs/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz
      if [[ "${TYPE}" == "GENOME" ]]; then
        bed=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/fastq/GIABS/genome/HG001/gencode.r33.hg19.transcripts_padded_150.bed
      else
        bed=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/fastq/GIABS/exome/HG001/regions_padded_150.bed
      fi
  fi
  if [[ "${SAMPLE}" == "HG002" ]]; then
      ref_bed=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/refs/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed
      ref_vcf=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/refs/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz
      if [[ "${TYPE}" == "GENOME" ]]; then
        bed=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/fastq/GIABS/genome/HG001/gencode.r33.hg19.transcripts_padded_150.bed
      else
        bed=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/fastq/GIABS/exome/HG002/regions_padded_150.bed
      fi
  fi
  if [[ "${SAMPLE}" == "HG003" ]]; then
      ref_bed=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/refs/HG003_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed
      ref_vcf=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/refs/HG003_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v.3.3.2_highconf.vcf.gz
      if [[ "${TYPE}" == "GENOME" ]]; then
        bed=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/fastq/GIABS/genome/HG001/gencode.r33.hg19.transcripts_padded_150.bed
      else
        bed=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/fastq/GIABS/exome/HG003/regions_padded_150.bed
      fi
  fi
  if [[ "${SAMPLE}" == "HG004" ]]; then
      ref_bed=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/refs/HG004_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed
      ref_vcf=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/refs/HG004_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v.3.3.2_highconf.vcf.gz
      if [[ "${TYPE}" == "GENOME" ]]; then
        bed=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/fastq/GIABS/genome/HG001/gencode.r33.hg19.transcripts_padded_150.bed
      else
        bed=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/fastq/GIABS/exome/HG004/regions_padded_150.bed
      fi
  fi
  if [[ "${SAMPLE}" == "HG005" ]]; then
      ref_bed=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/refs/HG005_GRCh37_highconf_CG-IllFB-IllGATKHC-Ion-SOLID_CHROM1-22_v.3.3.2_highconf_noMetaSV.bed
      ref_vcf=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/refs/HG005_GRCh37_highconf_CG-IllFB-IllGATKHC-Ion-SOLID_CHROM1-22_v.3.3.2_highconf.vcf.gz
      if [[ "${TYPE}" == "GENOME" ]]; then
        bed=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/fastq/GIABS/genome/HG001/gencode.r33.hg19.transcripts_padded_150.bed
      else
        bed=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/fastq/GIABS/exome/HG005/regions_padded_150.bed
      fi
  fi

  DIR_PREFIX=$(basename "${CURR_VCF}" .vcf)
  HAPPY=/home/lmbs02/bio/biosoft/happy/hap.py-build/bin/hap.py
  REF=/home/lmbs02/bio/databases/referenses/hg19_37/gencode/GRCh37_p13_33_release/GRCh37.primary_assembly.genome.fa
  REF_RTG=/home/lmbs02/bio/databases/referenses/hg19_37/gencode/GRCh37_p13_33_release/GRCh37.primary_assembly.genome.rtg/
  RTG_SH=/home/lmbs02/bio/biosoft/rtg/rtg-core-non-commercial-3.10.1/rtg
  EVAL_DIR="/home/lmbs02/bio/work/process_panels/NIST_EXOMES/COMPARSION/${DIR_PREFIX}"
  rm -rf $EVAL_DIR
  mkdir $EVAL_DIR

  python2 $HAPPY $ref_vcf $CURR_VCF -r $REF -f $ref_bed --threads 26 --engine vcfeval -T $bed -o "${EVAL_DIR}/report" --engine-vcfeval-template $REF_RTG --engine-vcfeval-path $RTG_SH
  echo $SAMPLE $TYPE $ALIGNER $CALLER 
  echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
}
export -f make_comparsion

for DIR in $(find "${WORK_DIR}" -maxdepth 1 -mindepth 1 -type d -not -wholename "${WORK_DIR}")
do
  TMP=$(basename "${DIR}")
  SAMPLE=$(cut -d '_' -f 1 <<< "${TMP}")
  TYPE=$(cut -d '_' -f 2 <<< "${TMP}")
  echo $SAMPLE $TYPE
  BAM="${DIR}/${SAMPLE}_dub_marked.bam"
  METRICS="${DIR}/${SAMPLE}.HS.metrics"

  if [[ "${SAMPLE}" == "HG001" ]]; then
      if [[ "${TYPE}" == "GENOME" ]]; then
        bed=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/fastq/GIABS/genome/HG001/gencode.r33.hg19.transcripts_padded_150.bed
      else
        bed=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/fastq/GIABS/exome/HG001/regions_padded_150.bed
      fi
  fi
  if [[ "${SAMPLE}" == "HG002" ]]; then
      if [[ "${TYPE}" == "GENOME" ]]; then
        bed=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/fastq/GIABS/genome/HG001/gencode.r33.hg19.transcripts_padded_150.bed
      else
        bed=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/fastq/GIABS/exome/HG002/regions_padded_150.bed
      fi
  fi
  if [[ "${SAMPLE}" == "HG003" ]]; then
      if [[ "${TYPE}" == "GENOME" ]]; then
        bed=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/fastq/GIABS/genome/HG001/gencode.r33.hg19.transcripts_padded_150.bed
      else
        bed=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/fastq/GIABS/exome/HG003/regions_padded_150.bed
      fi
  fi
  if [[ "${SAMPLE}" == "HG004" ]]; then
      if [[ "${TYPE}" == "GENOME" ]]; then
        bed=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/fastq/GIABS/genome/HG001/gencode.r33.hg19.transcripts_padded_150.bed
      else
        bed=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/fastq/GIABS/exome/HG004/regions_padded_150.bed
      fi
  fi
  if [[ "${SAMPLE}" == "HG005" ]]; then
      if [[ "${TYPE}" == "GENOME" ]]; then
        bed=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/fastq/GIABS/genome/HG001/gencode.r33.hg19.transcripts_padded_150.bed
      else
        bed=/home/lmbs02/bio/work/process_panels/NIST_EXOMES/fastq/GIABS/exome/HG005/regions_padded_150.bed
      fi
  fi
  list=${bed/.bed/.interval_list}
  java -jar $PICARD BedToIntervalList I=$bed O=$list SD=$DICT
  java -jar -Xmx110g  $PICARD CollectHsMetrics R=$REF BAIT_INTERVALS=$list TARGET_INTERVALS=$list I=$BAM O=$METRICS NEAR_DISTANCE=0
done

#parallel --delay 10 --load 80% --noswap --memfree 20G --jobs 1 align_with_isaac4 ::: "${DIRS}" ::: "${ISAAC4}" ::: "${SAMTOOLS}" ::: "${REF}"

#parallel --delay 2.5 --load 80% --noswap --memfree 20G --jobs 1 merge_lane_bams ::: "${DIRS}" ::: "${SAMTOOLS}"

#parallel --delay 2.5 --load 80% --noswap --memfree 20G --jobs 1 sort_after_align ::: "${DIRS}" ::: "${SAMTOOLS}"

#parallel --delay 2.5 --load 80% --noswap --memfree 20G --jobs 1 mark_dublicates ::: "${DIRS}" ::: "${PICARD}" ::: "${SAMTOOLS}"

#parallel --delay 2.5 --load 80% --noswap --memfree 20G --jobs 10 base_recalibration ::: "${DIRS}" ::: "${GATK4}" ::: "${REF}" ::: "${MILS_INDELS}" ::: "${GENOMES}" ::: "${DBSNP}"

#parallel --delay 2.5 --load 80% --noswap --memfree 20G --jobs 10 apply_recalibration ::: "${DIRS}" ::: "${GATK4}" ::: "${REF}"

#conda activate gatk4141
#time  parallel --delay 1 --load 80% --noswap --memfree 20G --jobs 5 gtk4_hc ::: "${DIRS}" ::: "${GATK4}" ::: "${REF}"
#conda deactivate

#parallel --delay 1 --load 90% --noswap --memfree 10G --jobs 5 strelka_call ::: "${DIRS}" ::: "${STRELKA}" ::: "${MANTA}" ::: "${REF}"

#parallel --delay 5 --load 90% --noswap --memfree 10G --jobs 5 freebayes_call ::: "${DIRS}" ::: "${FREEBAYES}" ::: "${GATK4}" ::: "${REF}"

#parallel --delay 1 --load 80% --noswap --memfree 20G --jobs 1 deepvariant_call ::: "${DIRS}" ::: "${DEEPVARIANT}" ::: "${REF}" ::: "${DEEPVARIANT_WES_MODEL}"