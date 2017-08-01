#!/bin/sh

# wrapper for quadtrim to give control over file-naming 
# This script takes 3 or more arguments
# example 
#quadtrim.sh HCTJGALXX_6_161117_FR07923227_Other__R_160818_SHACLA_DNA_M002_R1.00001.fastq HFVNHALXX_5_161117_FR07923226_Other__R_160818_SHACLA_DNA_M002_R2.00001.fastq quadtrim -d sheep

R1=$1
R2=$2

shift 2
quadtrim_call=$@

# when quadtrim is run it does this : 

#HFVNHALXX_5_161117_FR07923226_Other__R_160818_SHACLA_DNA_M002_R1.00001.fastq
#HFVNHALXX_5_161117_FR07923226_Other__R_160818_SHACLA_DNA_M002_R2.00001.fastq
#
#yields
#
#HFVNHALXX_5_161117_FR07923226_Other__R_160818_SHACLA_DNA_M002_R1.00001-pass.fq
#HFVNHALXX_5_161117_FR07923226_Other__R_160818_SHACLA_DNA_M002_R2.00001-pass.fq
#HFVNHALXX_5_161117_FR07923226_Other__R_160818_SHACLA_DNA_M002_RX.00001-discard.fq
#HFVNHALXX_5_161117_FR07923226_Other__R_160818_SHACLA_DNA_M002_RX.00001-singleton.fq
#
# which this script renames as
#
#HFVNHALXX_5_161117_FR07923226_Other__R_160818_SHACLA_DNA_M002_R1.00001.pass
#HFVNHALXX_5_161117_FR07923226_Other__R_160818_SHACLA_DNA_M002_R2.00001.pass
#HFVNHALXX_5_161117_FR07923226_Other__R_160818_SHACLA_DNA_M002_RX.00001.pass - i.e. singletons
#HFVNHALXX_5_161117_FR07923226_Other__R_160818_SHACLA_DNA_M002_RX.00001.discard - i.e. singletons
# 
# it also yields  (e.g.)
#HFVNHALXX_5_161117_FR07923226_Other__R_160818_SHACLA_DNA_M002_RX.00001.quadtrim.stdout
#HFVNHALXX_5_161117_FR07923226_Other__R_160818_SHACLA_DNA_M002_RX.00001.quadtrim.stderr
#
# 
R1_base=`basename $R1 .fastq`
R2_base=`basename $R2 .fastq`
dir=`dirname $R1`
chunk=`echo $R1_base | awk -F\. '{print $NF}' -`

#run quadtrim
$quadtrim_call $R1 $R2 > $R1.quadtrim.stdout 2>$R1.quadtrim.stderr


mv $dir/${R1_base}-pass.fq $dir/${R1_base}.pass
mv $dir/${R2_base}-pass.fq $dir/${R2_base}.pass
singles_file=`ls $dir/*.${chunk}-singleton.fq`
singles_base=`basename $singles_file -singleton.fq`
mv $singles_file $dir/${singles_base}.pass
