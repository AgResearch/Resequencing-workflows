#!/bin/sh

BUILD_ROOT=/dataset/AG_1000_sheep/scratch/gatk_run_032016
export BUILD_ROOT
BUILD_BIN=/dataset/AG_1000_sheep/active/bin/resequencing_pipeline/Resequencing-workflows
export BUILD_BIN

mkdir -p $BUILD_ROOT/link_farm

function set_up_link_farm_ag_XXXXOLDXXXX() {
   mkdir -p $BUILD_ROOT/link_farm

   echo "links in $BUILD_ROOT/link_farm"

   # add AG bams
   for animal in `cat AG_for_gatk_run.txt`; do
      found=0
      for source_folder in /dataset/AG_1000_sheep/archive/2015_processing_results /dataset/AG_1000_sheep_overflow/archive/2015_processing_results /dataset/AG_1000_sheep/scratch/general_processing_062015  /dataset/AG_1000_sheep_overflow/scratch/overflow_processing_112015; do
      #echo $source_folder
         #echo "checking $source_folder/$animal/${animal}_realigned.bam"
         if [ -f $source_folder/$animal/${animal}_realigned.bam ]; then
            ln -fs $source_folder/$animal/${animal}_realigned.bam $BUILD_ROOT/link_farm/${animal}_realigned.bam
            ln -fs $source_folder/$animal/${animal}_realigned.bam.bai $BUILD_ROOT/link_farm/${animal}_realigned.bam.bai
            found=1
         fi
      done
      if [ $found == "0" ]; then
         echo $animal not found
      fi
   done  
}

function set_up_link_farm_ag() {
   # the lis of animals is compiled from
   # ../Resequencing-workflows/export_animals.py -n -d /tmp `cat animals_list.txt`
   # - this uses the list of all animals provided by count_animals.py, and
   # generates a statement to copy to disk - we can parse that
   # to get the paths to all of the animals
   # ../Resequencing-workflows/export_animals.py -n -d /tmp `cat ../data_transfer_032016/animals_list.txt`
   # this is then checked, and used to create the link-fram
   mkdir -p $BUILD_ROOT/link_farm

   for path in `cat ag_files.txt`; do
      linkname=`basename $path`
      
      # check existence
      ls  $path > /dev/null 2>&1
      if [ $? != 0 ]; then
         echo "WARNING  $path does not exist"
      fi
      ls  ${path}.bai > /dev/null 2>&1
      if [ $? != 0 ]; then
         echo "WARNING  ${path}.bai does not exist"
      fi
      ln -s $path $BUILD_ROOT/link_farm/$linkname
      ln -s ${path}.bai $BUILD_ROOT/link_farm/${linkname}.bai
   done
}

function set_up_link_farm_csiro() {
   mkdir -p $BUILD_ROOT/link_farm

   #for bam in /dataset/AG_1000_sheep/scratch/csiro_sheep/*_renamed.bam ; do
   for bam in /mnt/sheep_1000_temp/csiro/isgc/bam/*_renamed.bam ; do
      base=`basename $bam`
      ln -s $bam $BUILD_ROOT/link_farm/$base
      ln -s ${bam}.bai $BUILD_ROOT/link_farm/${base}.bai
   done
}

function set_up_link_farm_melbourne() {
   mkdir -p $BUILD_ROOT/link_farm

   #for bam in /dataset/AG_1000_sheep/scratch/melbourne_sheep/*.bam ; do
   for bam in /mnt/sheep_1000_temp/melbourne/bam/*.bam ; do
      base=`basename $bam`
      ln -s $bam $BUILD_ROOT/link_farm/$base
      ln -s ${bam}.bai $BUILD_ROOT/link_farm/${base}.bai
   done
}

function set_up_link_farm_exclusionsXXXXXOLDXXXXX() {
   for animal in `cat melbourne_exclusions.txt`; do
     ls $BUILD_ROOT/link_farm/${animal}*.bam* > /dev/null 2>&1 
     if [ $? == 0 ]; then
        echo "excluding $BUILD_ROOT/link_farm/${animal}*.bam*"
        rm $BUILD_ROOT/link_farm/${animal}*.bam*
     else
        echo "oops cant exclude $animal !?"
     fi 
     
   done
}

function set_up_link_farm_exclusions() {
   echo "no exclusions this time"
}

function run_wrapper() {
   INTERVALS_FILE=/dataset/AG_1000_sheep/active/bin/resequencing_pipeline/vcf_build_032016/35Mbp_intervals_with_10kpadding
   #tardis.py -k -w -c 1 -d $BUILD_ROOT -hpctype local $BUILD_BIN/gatk_wrapper.sh _condition_text_input_/dataset/AG_1000_sheep/scratch/gatk_run/35Mbp_intervals1 _condition_text_output_$BUILD_ROOT/not_1000_sheep1.vcf _condition_text_output_$BUILD_ROOT/not_1000_sheep1.log
   #tardis.py -k -w -c 1 -t $BUILD_BIN/gatk_condor_template.txt -d $BUILD_ROOT -hpctype condor $BUILD_BIN/gatk_wrapper.sh _condition_text_input_/dataset/AG_1000_sheep/scratch/gatk_run/35Mbp_intervals2 _condition_text_output_$BUILD_ROOT/not_1000_sheep2.vcf _condition_text_output_$BUILD_ROOT/not_1000_sheep2.log 8G
   #tardis.py -k -w -c 1 -t $BUILD_BIN/gatk_condor_template.txt -d $BUILD_ROOT -hpctype condor $BUILD_BIN/gatk_wrapper.sh _condition_text_input_/dataset/AG_1000_sheep/scratch/gatk_run/35Mbp_intervals1.1 _condition_text_output_$BUILD_ROOT/not_1000_sheep1.1.vcf _condition_text_output_$BUILD_ROOT/not_1000_sheep1.1.log 8G

   # 3/2016 runs
   # dry run
   #tardis.py -k -w -c 1 -v -t $BUILD_BIN/gatk_condor_template.txt -d $BUILD_ROOT -hpctype condor $BUILD_BIN/gatk_wrapper.sh _condition_text_input_$BUILD_BIN/35Mbp_intervals1.1 _condition_text_output_$BUILD_ROOT/1000_sheep.vcf.cat _condition_text_output_$BUILD_ROOT/1000_sheep.logi.cat.log 8G
   # run with crook padding
   #tardis.py -k -w -c 1 -t $BUILD_BIN/gatk_condor_template.txt -d $BUILD_ROOT -hpctype condor $BUILD_BIN/gatk_wrapper.sh _condition_text_input_$BUILD_BIN/35Mbp_intervals  _condition_text_output_$BUILD_ROOT/1000_sheep.vcf.cat _condition_text_output_$BUILD_ROOT/1000_sheep.log.cat.log 8G
   tardis.py -k -w -c 1 -t $BUILD_BIN/gatk_condor_template.txt -d $BUILD_ROOT -hpctype condor $BUILD_BIN/gatk_wrapper.sh _condition_text_input_$INTERVALS_FILE  _condition_text_output_$BUILD_ROOT/1000_sheep.vcf _condition_text_output_$BUILD_ROOT/1000_sheep_vcfbuild.log 8G
}

#set_up_link_farm_ag
#set_up_link_farm_csiro
#set_up_link_farm_melbourne
#set_up_link_farm_exclusions
run_wrapper

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
