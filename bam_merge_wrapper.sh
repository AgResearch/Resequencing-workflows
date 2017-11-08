#!/bin/bash

# this wrapper is used to make the bamtools merge interface work 
# like the sambamba interface - i.e. so can just pass a list of input 
# files , without having to prefix each one with -in
out_bam=$1
in_bams=""
shift 
while [ ! -z "$1" ]; do
   in_bams="$in_bams -in $1"
   shift
done
command="bamtools merge $in_bams -out $out_bam"
$command
