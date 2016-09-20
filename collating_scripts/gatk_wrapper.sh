#!/bin/sh
BUILD_ROOT=/dataset/AG_1000_sheep/scratch/gatk_run_032016
export BUILD_ROOT
BUILD_BIN=/dataset/AG_1000_sheep/active/bin/resequencing_pipeline/Resequencing-workflows
export BUILD_BIN


# wrapper for all-bams run of GATK , to enable setting 
# optimal environment for a run 
# takes two args - intervals file, and output file
# 
INTERVALS_FILE=$1
OUTPUT_FILE=$2
LOG_FILE=$3
mx=$4
if [ -z "$mx" ]; then
   mx=10G
fi

#cd $BUILD_ROOT

# sanity check we are somewhere sensible
#if [ ! -d ./link_farm ]; then
#   echo "link_farm missing !"
#   exit 1
#fi

bam_args=""
for bamfile in $BUILD_ROOT/link_farm/*.bam ; do
   bam_args="$bam_args -I $bamfile "
done
interval=`cat $INTERVALS_FILE`

echo "java -Xmx$mx -jar /dataset/AG_1000_sheep/active/bin/GATK-3.4-46/GenomeAnalysisTK.jar -R /dataset/datacache/scratch/ensembl/oar3.1/indexes/Ovis_aries.Oar_v3.1.dna_sm.toplevel.fa -T UnifiedGenotyper $bam_args -o $OUTPUT_FILE -out_mode EMIT_VARIANTS_ONLY -stand_call_conf 20 -stand_emit_conf 20 -A FisherStrand -A StrandOddsRatio -A StrandBiasBySample -rf BadCigar --genotype_likelihoods_model BOTH -L $interval -nt 2 -nct 2 -log $LOG_FILE"

java -Xmx$mx -jar /dataset/AG_1000_sheep/active/bin/GATK-3.4-46/GenomeAnalysisTK.jar -R /dataset/datacache/scratch/ensembl/oar3.1/indexes/Ovis_aries.Oar_v3.1.dna_sm.toplevel.fa -T UnifiedGenotyper $bam_args -o $OUTPUT_FILE -out_mode EMIT_VARIANTS_ONLY -stand_call_conf 20 -stand_emit_conf 20 -A FisherStrand -A StrandOddsRatio -A StrandBiasBySample -rf BadCigar --genotype_likelihoods_model BOTH -L $interval -nt 2 -nct 2 -log $LOG_FILE


exit






#For usage see https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_genotyper_UnifiedGenotyper.php

#Here's what we want

java -Xmx10g -jar /dataset/AG_1000_sheep/active/bin/GATK-3.4-46/GenomeAnalysisTK.jar
-R /dataset/datacache/scratch/ensembl/oar3.1/indexes/Ovis_aries.Oar_v3.1.dna_sm.toplevel.fa
-T UnifiedGenotyper
-I bam_file1 -I bam_file2 -I bam_file3
-o my_slice.vcf
-out_mode EMIT_VARIANTS_ONLY
-stand_call_conf 20
-stand_emit_conf 20
-A FisherStrand
-A StrandOddsRatio
-A StrandBiasBySample
-rf BadCigar
--genotype_likelihoods_model BOTH
-L chrom:start-stop
-nt 2
-nct 2
-log my_lice.log

#-I needs to feature all SGD bam files

#here's an example for just 2 bam files
java -Xmx10g -jar /dataset/AG_1000_sheep/active/bin/GATK-3.4-46/GenomeAnalysisTK.jar -R /dataset/datacache/scratch/ensembl/oar3.1/indexes/Ovis_aries.Oar_v3.1.dna_sm.toplevel.fa -T UnifiedGenotyper -I /dataset/AG_1000_sheep/scratch/general_processing_062015/NZCMPF100018040066/NZCMPF100018040066_realigned.bam -I /dataset/AG_1000_sheep/scratch/general_processing_062015/NZCMPM100018040065/NZCMPM100018040065_realigned.bam -o 1_36M_test2.vcf -out_mode EMIT_VARIANTS_ONLY -stand_call_conf 20 -stand_emit_conf 20 -A FisherStrand -A StrandOddsRatio -A StrandBiasBySample -rf BadCigar --genotype_likelihoods_model BOTH -L 1:1000000-36000000 -nt 2 -nct 2 -log UG2.log
#runtime: 60min (last time it took 30min)
#RAM usage: 3 GB

#there is an interval file available with 35Mbp slices, 35Mbp_intervals
/dataset/AG_1000_sheep/scratch/gatk_run/35Mbp_intervals
#that's 90 35Mbp slices and some 5600 smaller slices (chrUn).
#this can be submitted as one long list or more sensibly as 1000s of small jobs for the farm.




