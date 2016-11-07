#!/bin/sh

#echo "get_fastq.sh is processing $@"

# examples :   
#./get_fastq.sh link_files /dataset/farmIQ_2012/archive/Nov_2012/C15BA/Sample_907834/907834_NoIndex_L002_R1_001.fastq.gz  /dataset/AG_1000_sheep_2017/scratch/processing/NZROMM100017134122/907834_NoIndex_L002_R1_001.fastq.gz
#./get_fastq.sh sra /dataset/AG_1000_sheep_2017/scratch/SRS1709518/SRR4290367_1.fastq.gz /dataset/AG_1000_sheep_2017/scratch/processing/SRS1709518/SRR4290367_1.fastq.gz

FASTQ_DUMP=/dataset/bioinformatics_dev/active/sra-tools/sratoolkit.2.4.5-2-centos_linux64/bin/fastq-dump
method=$1

set -x
if [ $method == "link_files" ]; then
   target=$2
   link_name=$3
   ln -s $target $link_name
elif [ $method == "sra" ]; then 
   target=$2
   output_name=$3

   # from the target, figure out the sra file name
   sra_dir=`dirname $target`
   sra_base=`basename $target _1.fastq.gz`
   sra_file=$sra_dir/${sra_base}.sra/${sra_base}.sra

   # from the output_name get the output folder
   output_folder=`dirname $output_name`

   # get the fastq files (- this command will get both pairs)
   $FASTQ_DUMP --split-files --gzip --outdir $output_folder $sra_file
fi
