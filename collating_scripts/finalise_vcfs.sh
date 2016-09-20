#!/bin/sh

INTERVALS_FILE_WITHOUT_PADDING=/dataset/AG_1000_sheep/active/bin/resequencing_pipeline/vcf_build_032016/35Mbp_intervals_with_nopadding
# (obtained as ./infoseq_to_slices.py /dataset/AG_1000_sheep/scratch/gatk_run/oarv3.1.infoseq > 35Mbp_intervals_with_nopadding)
BUILD_ROOT=/dataset/AG_1000_sheep/scratch/gatk_run_032016
BUILD_BIN=/dataset/AG_1000_sheep/active/bin/resequencing_pipeline/Resequencing-workflows
VCF_FOLDER=/dataset/AG_1000_sheep/scratch/gatk_run_032016/tardis_gd2_ss

function generate_unpadding_script() {
#
# generate a series of commands to "un-pad" sub-region vcf files using commands like 
# vcftools --vcf /dataset/AG_1000_sheep/scratch/gatk_run_032016/tardis_gd2_ss/1000_sheep.00001.vcf --chr 1 --from-bp 1 --to-bp 35000000 --recode --recode-INFO-all --out /dataset/AG_1000_sheep/scratch/gatk_run_032016/unpadded/1000_sheep.00001.vcf
# Ths is used to process a series of sub-region VCF files that cover overlapping sub-regions of a chromosome (which in turn result from 
# a "paralellised" run of gatk on a compute cluster, with each node analysing variation across N animals 
#, for a specific genomic region. Overlapping regions used to avoid edge effects)
# The listing should be directed to a command file, which can then be run like this : 
# tardis.py -w -c 1 -t $BUILD_BIN/gatk_condor_template.txt -d $BUILD_ROOT -hpctype condor source _condition_text_input_unpad_commands.txt 1\> _condition_text_output_unpadding.stdout 2\>_condition_text_output_unpadding.stderr
#
chunk_number=1
for record in `cat $INTERVALS_FILE_WITHOUT_PADDING`; do
   # construct a vcftools command to list the vcf file without padding 
   # this one doesn't do a very good job .....
   #vcf_file=`echo $chunk_number | awk '{printf("1000_sheep.%05d.vcf",$1)}' -`
   #command="vcf-query $VCF_FOLDER/$vcf_file -f '%LINE' -r $record > $BUILD_ROOT/unpadded/$vcf_file"
   #echo $command

   # this one does a better job
   chr=`echo $record | awk -F\: '{print $1}' -`
   from=`echo $record | awk -F\: '{print $2}' - | awk -F- '{print $1}' -`
   to=`echo $record | awk -F\: '{print $2}' - | awk -F- '{print $2}' -`
   echo $chunk_number,$chr,$from,$to | awk -F, '{printf("vcftools --vcf /dataset/AG_1000_sheep/scratch/gatk_run_032016/tardis_gd2_ss/1000_sheep.%05d.vcf --chr %s --from-bp %s --to-bp %s --recode --recode-INFO-all --out /dataset/AG_1000_sheep/scratch/gatk_run_032016/unpadded/1000_sheep.%05d.vcf\n",$1,$2,$3,$4,$1)}' -

   let chunk_number=chunk_number+1
done
}

function generate_concat_listfiles() {
# generate a listfile, containing all of the (non-overlapping) sub-region VCF files for each chromosome. 
# Such a listfile can be used as an argument to vcf-concat
concat_type=$1
if [ $concat_type == "unpadded" ]; then
   chunk_number=1
   for record in `cat $INTERVALS_FILE_WITHOUT_PADDING`; do
      # create a seperate concat listfile for each chromosome

      chr=`echo $record | awk -F\: '{print $1}' -`
      listfile="other_vcf.list"
      case $chr in (1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|23|24|25|26|X|MT)
         listfile="${chr}_vcf.list"
      esac

      echo $chunk_number | awk '{printf("/dataset/AG_1000_sheep/scratch/gatk_run_032016/unpadded/1000_sheep.%05d.vcf.recode.vcf\n",$1)}' - >> $listfile

      let chunk_number=chunk_number+1
   done
elif [ $concat_type == "filtered" ]; then
   chunk_number=1
   for record in `cat $INTERVALS_FILE_WITHOUT_PADDING`; do
      # create a seperate concat listfile for each chromosome

      chr=`echo $record | awk -F\: '{print $1}' -`
      listfile="other_vcf_filtered.list"
      case $chr in (1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|23|24|25|26|X|MT)
         listfile="${chr}_vcf_filtered.list"
      esac

      echo $chunk_number | awk '{printf("/dataset/AG_1000_sheep/scratch/gatk_run_032016/filtered/1000_sheep.%05d.recode.filtered.vcf\n",$1)}' - >> $listfile

      let chunk_number=chunk_number+1
   done
fi
}

function run_concats() {
# this method used to concatenate VCF files from non-overlapping sub-regions of a chromosome, using
# a listfile obtained as above
   # this doesn't actually run them - just produces a command-file
   concat_type=$1
   if [ $concat_type == "unpadded" ]; then
      for file in *_vcf.list; do
         chr=`basename $file _vcf.list`
         echo "nohup vcf-concat -f $file > $BUILD_ROOT/${chr}.vcf &"
      done
   elif [ $concat_type == "filtered" ]; then
      for file in *_vcf_filtered.list; do
         chr=`basename $file _vcf_filtered.list`
         echo "nohup vcf-concat -f $file > $BUILD_ROOT/${chr}.filtered.vcf &"
      done
   fi
}

repair_zero_length_vcf() {
# this method a one-off to repair some VCFs that were incomplete due to 
# cluster outage
# manually run these :
#/dataset/AG_1000_sheep/scratch/gatk_run_032016/tardis_gd2_ss/1000_sheep.04178.vcf
#/dataset/AG_1000_sheep/scratch/gatk_run_032016/tardis_gd2_ss/1000_sheep.01829.vcf
#/dataset/AG_1000_sheep/scratch/gatk_run_032016/tardis_gd2_ss/1000_sheep.05219.vcf
#/dataset/AG_1000_sheep/scratch/gatk_run_032016/tardis_gd2_ss/1000_sheep.03389.vcf
#/dataset/AG_1000_sheep/scratch/gatk_run_032016/tardis_gd2_ss/1000_sheep.03862.vcf
#/dataset/AG_1000_sheep/scratch/gatk_run_032016/tardis_gd2_ss/1000_sheep.00760.vcf
#/dataset/AG_1000_sheep/scratch/gatk_run_032016/tardis_gd2_ss/1000_sheep.02766.vcf
#/dataset/AG_1000_sheep/scratch/gatk_run_032016/tardis_gd2_ss/1000_sheep.04091.vcf
#/dataset/AG_1000_sheep/scratch/gatk_run_032016/tardis_gd2_ss/1000_sheep.01809.vcf

rm unpad_commands_subset.txt
for chunk in 4178 1829 5219 3389 3862 760 2766 4091 1809; do
   #echo "running $VCF_FOLDER/run${chunk}.shxxxx"
   #$VCF_FOLDER/run${chunk}.shxxxx
   name_frag=`echo $chunk | awk '{printf("1000_sheep.%05d.vcf\n",$1)}' -`
   grep $name_frag unpad_commands.txt >> unpad_commands_subset.txt
   source unpad_commands_subset.txt
done

}

generate_filter_commands() {
# generate commands which will run all of the (overlapping) sub-region VCF files through a chain of filters
# - redirect output to a command file and then run it with e.g. 
# tardis.py -w -c 1 -t $BUILD_BIN/gatk_condor_template.txt -d $BUILD_ROOT -hpctype condor source _condition_text_input_filter_commands.txt 1\> _condition_text_output_filtering.stdout 2\>_condition_text_output_filtering.stderr

   FILTERS_BIN=/dataset/AG_1000_sheep/active/bin/resequencing_pipeline/vcf_build_032016/filters
   BUILD_DIR=/dataset/AG_1000_sheep/scratch/gatk_run_032016/filtered
   for file in /dataset/AG_1000_sheep/scratch/gatk_run_032016/unpadded/*.vcf; do
      BASE=`basename $file .vcf.recode.vcf`
      outfile=$BUILD_DIR/${BASE}.recode.filtered.vcf 
      echo "$FILTERS_BIN/ag_run_filters.sh $file $outfile"
   done
}

fun_filtering() {
   BUILD_BIN=/dataset/AG_1000_sheep/active/bin/resequencing_pipeline/Resequencing-workflows
   BUILD_ROOT=/dataset/AG_1000_sheep/scratch/gatk_run_032016
   #tardis.py -w -c 1 -t $BUILD_BIN/gatk_condor_template.txt -d $BUILD_ROOT -hpctype condor source _condition_text_input_filter_commands.txt 1\> _condition_text_output_filtering.stdout 2\>_condition_text_output_filtering.stderr

   # this version of the command splices a header of the chunk number into the collated stdout and stderr files to help with tracking down
   # which chunk generated which error message
   tardis.py -w -c 1 -k -t $BUILD_BIN/gatk_condor_template.txt -d $BUILD_ROOT -hpctype condor source _condition_text_input_filter1_commands.txt 1\> _condition_text_output_filtering1.stdout 2\> \>\(echo _condition_output_interval \> _condition_text_output_filtering1.stderr \; cat \>\> _condition_text_output_filtering1.stderr \)
}


######################################################################################################
# uncomment and run the above to finalise the vcf's 
######################################################################################################
#generate_unpadding_script
#generate_concat_listfiles unpadded
run_concats unpadded
#repair_zero_length_vcf
#generate_filter_commands
#
# noet no run_filter method - filtering run manually by run_filtering.sh  
#generate_concat_listfiles filtered
#run_concats filtered
