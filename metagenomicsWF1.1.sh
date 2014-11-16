#!/bin/sh

RUN1DIR=/dataset/metagdata/archive/nzgl/140627_M02810_0023_000000000-A856J/processed_trimmed
RUN2DIR=/dataset/metagdata/archive/nzgl/140714_M02810_0028_000000000-AAJE2/processed_trimmed
BUILDDIR=/dataset/metagdata/scratch/nzgl/buildR2

make -f metagenomicsWF1.1.mk -d --no-builtin-rules -j 8 run1dir=$RUN1DIR run2dir=$RUN2DIR builddir=$BUILDDIR  R2tax_assignments  > buildR2filtered.log 2>&1 
module load qiime
make -f metagenomicsWF1.1.mk -d --no-builtin-rules -j 8 run1dir=$RUN1DIR run2dir=$RUN2DIR builddir=$BUILDDIR  R2tax_assignment_summaries  >> buildR2filtered.log 2>&1
