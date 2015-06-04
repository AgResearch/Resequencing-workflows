#!/bin/sh

# wrap a call to quadtrim, so that we run it in a sub-dir thats passed in 
QUADTRIM=/dataset/hiseq/active/bin/quadtrim
R1=$1
R2=$2
quad_dir=$3

mkdir $quad_dir
cd $quad_dir

$QUADTRIM -d sheep $R1 $R2

# output examples : 
#1016456575_GTGAAA_L001_R1_001-pass.fq  1016456575_GTGAAA_L001_RX_001-discard.fq    quad_wrapper.sh
#1016456575_GTGAAA_L001_R2_001-pass.fq  1016456575_GTGAAA_L001_RX_001-singleton.fq
