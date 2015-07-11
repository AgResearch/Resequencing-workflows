#!/bin/sh

DATADIR=/dataset/AG_1000_sheep/active/bin/resequencing_pipeline/Resequencing-workflows/ResequencingWF1.5.testdata
OUTDIR=/dataset/AG_1000_sheep/active/bin/resequencing_pipeline/Resequencing-workflows/ResequencingWF1.5.testout
TARDIS=~/galaxy/hpc/dev/tardis.py


$TARDIS -w -c 1000 -hpctype local -d $OUTDIR /dataset/hiseq/active/bin/quadtrim _condition_paired_fastq_input_$DATADIR/C6KHFANXX-1418-37-05-01_ATGTCA_L001_R1_001.fastq.gz  _condition_paired_fastq_input_$DATADIR/C6KHFANXX-1418-37-05-01_ATGTCA_L001_R2_001.fastq.gz  '_condition_fastq_product_\S*?_R1_\S*?\d{5}-pass.fq,/dataset/AG_1000_sheep/active/bin/resequencing_pipeline/Resequencing-workflows/ResequencingWF1.5.testout/C6KHFANXX-1418-37-05-01_ATGTCA_L001_R1_001.fastq.quadtrim.gz'
