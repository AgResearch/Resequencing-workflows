#!/bin/sh

#echo "get_fastq.sh is processing $@"
# example : ./get_fastq.sh link_files /dataset/farmIQ_2012/archive/Nov_2012/C15BA/Sample_907834 /dataset/AG_1000_sheep_2017/scratch/processing/NZROMM100017134122/907834_NoIndex_L002_R1_001.fastq.gz

method=$1

set -x
if [ $method == "link_files" ]; then
   target=$2
   link_name=$3
   ln -s $target $link_name
fi
set +x


