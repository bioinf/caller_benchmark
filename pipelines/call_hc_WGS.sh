#!/bin/bash

REFERENCE='GRCh37.primary_assembly.genome.fa'
CALLINGREGIONS='gencode.v19.revamp.bed'
DBSNP='hg19_v0_Homo_sapiens_assembly19.dbsnp138.hg19.vcf.gz'
MILLS_INDELS='hg19_v0_Mills_and_1000G_gold_standard.indels.b37.hg19.vcf.gz'
GENOMES='hg19_v0_Homo_sapiens_assembly19.known_indels_20120518.hg19.vcf.gz'
KG_HICONF_INDELS='hg19_v0_1000G_phase1.snps.high_confidence.b37.hg19.vcf.gz'
HAPMAP='hg19_v0_hapmap_3.3.b37.hg19.vcf.gz'
KG_SITES='b37_1000G_phase3_v4_20130502.sites.hg19.vcf.gz'

i=$1

docker run -v /media/array/callers_proj/raw_data:/data broadinstitute/gatk:latest gatk --java-options "-Xmx8g" BaseRecalibrator -R /data/$REFERENCE -I /data/bams/WGS/$i --known-sites /data/bundle/$DBSNP --known-sites /data/bundle/$MILLS_INDELS --known-sites /data/bundle/$GENOMES -O /data/${i%%.dedup.bam}/${i%%.dedup.bam}.recal.table -L /data/$CALLINGREGIONS
docker run -v /media/array/callers_proj/raw_data:/data broadinstitute/gatk:latest gatk --java-options "-Xmx8g" ApplyBQSR -R /data/$REFERENCE -I /data/bams/WGS/$i -bqsr /data/${i%%.dedup.bam}/${i%%.dedup.bam}.recal.table -O /data/${i%%.dedup.bam}/${i%%.dedup.bam}.recal.bam -L /data/$CALLINGREGIONS
	
# Raw HaplotypeCaller
docker run -v /media/array/callers_proj/raw_data:/data broadinstitute/gatk:latest gatk --java-options "-Xmx8g"  HaplotypeCaller -R /data/$REFERENCE -I /data/${i%%.dedup.bam}/${i%%.dedup.bam}.recal.bam -bamout /data/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.bam -O /data/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.RAW.vcf -L /data/$CALLINGREGIONS
        
# CNN scoring and filtration
docker run -v /media/array/callers_proj/raw_data:/data broadinstitute/gatk:latest gatk --java-options "-Xmx8g" CNNScoreVariants -R /data/$REFERENCE -V /data/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.RAW.vcf -O /data/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.CNN1.vcf -tensor-type reference -L /data/$CALLINGREGIONS -ip 20
docker run -v /media/array/callers_proj/raw_data:/data broadinstitute/gatk:latest gatk --java-options "-Xmx8g" CNNScoreVariants -I /data/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.bam -R /data/$REFERENCE -V /data/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.CNN1.vcf -O /data/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.CNN2.vcf -tensor-type read_tensor -L /data/$CALLINGREGIONS -ip 20

docker run -v /media/array/callers_proj/raw_data:/data broadinstitute/gatk:latest gatk --java-options "-Xmx8g" FilterVariantTranches -V /data/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.CNN1.vcf --output /data/${i%%.dedup.bam}/${i%%.dedup.bam}_HC_1DCNN.vcf \
  --info-key CNN_1D --snp-tranche 99.9 --indel-tranche 99.0 \
  --resource /data/bundle/$MILLS_INDELS \
  --resource /data/bundle/$KG_HICONF_INDELS \
  --resource /data/bundle/$HAPMAP 
docker run -v /media/array/callers_proj/raw_data:/data broadinstitute/gatk:latest gatk --java-options "-Xmx8g" FilterVariantTranches -V /data/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.CNN2.vcf --output /data/${i%%.dedup.bam}/${i%%.dedup.bam}_HC_2DCNN.vcf \
  --resource /data/bundle/$MILLS_INDELS \
  --resource /data/bundle/$KG_HICONF_INDELS \
  --resource /data/bundle/$HAPMAP \
  --info-key CNN_2D --snp-tranche 99.9 --indel-tranche 99.5
	
# Hard filters
docker run -v /media/array/callers_proj/raw_data:/data broadinstitute/gatk:latest gatk --java-options "-Xmx8g" SelectVariants -V /data/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.RAW.vcf -select-type SNP -select-type MIXED -O /data/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.snps.vcf
docker run -v /media/array/callers_proj/raw_data:/data broadinstitute/gatk:latest gatk --java-options "-Xmx8g" VariantFiltration -V /data/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.snps.vcf -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" \
  -filter "SOR > 3.0" --filter-name "SOR3" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" \
  --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O /data/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.snps.flt.vcf
  	
docker run -v /media/array/callers_proj/raw_data:/data broadinstitute/gatk:latest gatk --java-options "-Xmx8g" SelectVariants -V /data/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.RAW.vcf -select-type INDEL -select-type MIXED -O /data/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.indels.vcf
docker run -v /media/array/callers_proj/raw_data:/data broadinstitute/gatk:latest gatk --java-options "-Xmx8g" VariantFiltration -V /data/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.indels.vcf -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" \
  -filter "FS > 200.0" --filter-name "FS200" -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" -O /data/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.indels.flt.vcf
	
docker run -v /media/array/callers_proj/raw_data:/data broadinstitute/gatk:latest gatk --java-options "-Xmx8g" MergeVcfs -I /data/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.snps.flt.vcf -I /data/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.indels.flt.vcf -O /data/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.hardfilter_tmp.vcf
	
docker run -v /media/array/callers_proj/raw_data:/data broadinstitute/gatk:latest gatk --java-options "-Xmx8g" CalculateGenotypePosteriors -V /data/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.hardfilter_tmp.vcf -supporting /data/bundle/$KG_SITES \
  -O /data/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.hardfilter_tmp2.vcf
	
docker run -v /media/array/callers_proj/raw_data:/data broadinstitute/gatk:latest gatk --java-options "-Xmx8g" VariantFiltration -V /data/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.hardfilter_tmp2.vcf --genotype-filter-expression "GQ < 20" --genotype-filter-name "GQ20" -O /data/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.hardfilter_tmp3.vcf

bcftools filter -s GQ20 -e "FORMAT/FT[*]!=''" /media/array/callers_proj/raw_data/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.hardfilter_tmp3.vcf > /media/array/callers_proj/raw_data/${i%%.dedup.bam}/${i%%.dedup.bam}_HC_HARDFILTER.vcf
