#!/bin/sh
# this version uses makefile 1.5 which uses quadtrim rather than flexbar

function get_opts() {

DRY_RUN=no
p1=""
p2=""
prestr=""
midstr=""
poststr=""
R1_FILE_PATTERN="*_L*_R1_*.fastq.gz"

help_text="
 examples : \n
 ./ResequencingWF1.5.sh  -n -S Sample_932593 -D /dataset/hiseq/active/140516_D00390_0041_BC4KPCACXX/bcl2fastq/Unaligned/Project_SalmonSheepWGS/Sample_932593 -G /dataset/Salmo_salar_2014/active/ICSASG_v1/Ssa_ASM_3.6.fa  -B /dataset/Salmo_salar_2014/scratch  -p 9 -q 1 -r 2 -s _R -t _0\n
 (the sample folders contains paired files like  \n
  932593_GTCCGC_L003_R1_005.fastq.gz, 932593_GTCCGC_L003_R1_007.fastq.gz, 932593_GTCCGC_L004_R1_002.fastq.gz\n
  932593_GTCCGC_L003_R2_005.fastq.gz, 932593_GTCCGC_L003_R2_007.fastq.gz, 932593_GTCCGC_L004_R2_002.fastq.gz\n
  etc \n

 ./ResequencingWF1.5.sh  -n -S NZCPWF000001391796 -D /dataset/hiseq/scratch/postprocessing_dev/150508_D00390_0226_AC6H4RANXX.processed/bcl2fastq/Project_WGS_NICHETRAITS/Sample_1016456575  -G /dataset/OARv3.0/active/current_version/sheep.v3.0.14th.final.fa -B /dataset/AG_1000_sheep/scratch/PHUAS_processing_052015 -T /dataset/AG_1000_sheep/scratch/PHUAS_processing_052015 -p 1016456575_GTGAAA -q 1 -r 2 -s _R -t _0\n

 ./ResequencingWF1.5.sh  -n -S NZCPWF000001391796_test -D /dataset/AG_1000_sheep/scratch/PHUAS_processing_052015/testdata -G /dataset/OARv3.0/active/current_version/sheep.v3.0.14th.final.fa -B /dataset/AG_1000_sheep/scratch/PHUAS_processing_052015 -T /dataset/AG_1000_sheep/scratch/PHUAS_processing_052015 -p 1016456575_GTGAAA -q 1 -r 2 -s _R -t _0\n

 ./ResequencingWF1.5.sh  -n -S NZCPWF000001391796_oldparms -D /dataset/hiseq/scratch/postprocessing_dev/150508_D00390_0226_AC6H4RANXX.processed/bcl2fastq/Project_WGS_NICHETRAITS/Sample_1016456575 -G /dataset/OARv3.0/active/current_version/sheep.v3.0.14th.final.fa -B /dataset/AG_1000_sheep/scratch/PHUAS_processing_052015 -T /dataset/AG_1000_sheep/scratch/PHUAS_processing_052015 -p 1016456575_GTGAAA -q 1 -r 2 -s _R -t _0\n

 ./ResequencingWF1.5.sh  -n -S NZCPWF000001391796_test -D /dataset/AG_1000_sheep/scratch/PHUAS_processing_052015/testdata -G /dataset/datacache/scratch/ensembl/oar3.1/indexes/Ovis_aries.Oar_v3.1.dna_sm.toplevel.fa -B /dataset/AG_1000_sheep/scratch/PHUAS_processing_052015 -T /dataset/AG_1000_sheep/scratch/PHUAS_processing_052015 -p 1016456575_GTGAAA -q 1 -r 2 -s _R -t _0\n

 ./ResequencingWF1.5.sh  -S NZCPWF000001391796 -D /dataset/hiseq/scratch/postprocessing_dev/150508_D00390_0226_AC6H4RANXX.processed/bcl2fastq/Project_WGS_NICHETRAITS/Sample_1016456575 -G /dataset/datacache/scratch/ensembl/oar3.1/indexes/Ovis_aries.Oar_v3.1.dna_sm.toplevel.fa -B /dataset/AG_1000_sheep/scratch/PHUAS_processing_052015 -T /dataset/AG_1000_sheep/scratch/PHUAS_processing_052015 -p 1016456575_GTGAAA -q 1 -r 2 -s _R -t _0\n

./ResequencingWF1.5.sh -n -S NZCPWF100017865294 -X "*TTAGGC_L*_R1_*.fastq.gz" -D /dataset/BLGsheep/archive/NZGL01418_1/Raw  -G /dataset/datacache/scratch/ensembl/oar3.1/indexes/Ovis_aries.Oar_v3.1.dna_sm.toplevel.fa -B /dataset/AG_1000_sheep/scratch/general_processing_062015 -T /dataset/AG_1000_sheep/scratch/general_processing_062015 -p C6KFHANXX-1418-2-05-01_TTAGGC -q 1 -r 2 -s _R -t _0
(using the -X argument to specify a pattern to match when marshalling files from the 
data folder - e.g. when there are a number of different barcoded outputs in a single folder)

"

while getopts ":nhS:B:D:G:X:T:p:q:r:s:t:" opt; do
  case $opt in
    n)
      DRY_RUN=yes
      ;;
    h)
      echo -e $help_text
      exit 0
      ;;
    T)
      TEMP_ROOT=$OPTARG
      ;;
    S)
      SAMPLE=$OPTARG
      ;;
    B)
      BUILD_ROOT=$OPTARG
      ;;
    D)
      DATA_DIR=$OPTARG
      ;;
    G)
      REF_GENOME=$OPTARG
      ;;
    X)
      R1_FILE_PATTERN=$OPTARG
      ;;
    p)
      prestr=$OPTARG
      ;;
    q)
      p1=$OPTARG
      ;;
    r)
      p2=$OPTARG
      ;;
    s)
      midstr=$OPTARG
      ;;
    t)
      poststr=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done
}

function check_opts() {
  if [ ! -d "$BUILD_ROOT" ]; then
    echo "BUILD_ROOT $BUILD_ROOT not found "
    exit 1
  fi
  if [ ! -d "$DATA_DIR" ]; then
    echo "DATA_DIR $DATA_DIR not found "
    exit 1
  fi
  if [ ! -d "$TEMP_ROOT" ]; then
    echo "TEMP_ROOT $TEMP_ROOT not found "
    exit 1
  fi
  if [ -z "$SAMPLE" ]; then
    echo "must specify sample name using -S option"
    exit 1
  fi
}

function echo_opts() {
  echo SAMPLE=$SAMPLE
  echo DATA_DIR=$DATA_DIR
  echo BUILD_ROOT=$BUILD_ROOT
  echo LANE_TARGETS=$LANE_TARGETS
  echo DRY_RUN=$DRY_RUN
}

function load_modules() {
  module load samtools/0.1.19
}


get_opts $@

check_opts

echo_opts

load_modules

##### from here inline code to run the processing
TEMP_DIR=$TEMP_ROOT/${SAMPLE}tmp
BUILD_DIR=$BUILD_ROOT/${SAMPLE}
if [ ! -d $TEMP_DIR ]; then
   mkdir $TEMP_DIR
fi

if [ ! -d $BUILD_DIR ]; then
   mkdir $BUILD_DIR
fi



# copy or init tardis config file 
if [ -f .tardishrc ]; then
   cp  .tardishrc $TEMP_DIR
else
echo "
[tardish]

[tardis_engine]
job_template_name=condor_send_env_job
shell_template_name=condor_shell
" > $TEMP_DIR/.tardishrc
fi

# get the lane monikers string needed by the makefile 
moniker_string=""
moniker_list_file=/tmp/${$}moniker_list.tmp
echo  > $moniker_list_file 
R1_PATH_PATTERN=$DATA_DIR/$R1_FILE_PATTERN
#for R1file in $DATA_DIR/*_L*_R1_*.fastq.gz; do 
for R1file in $R1_PATH_PATTERN; do 
   moniker=`basename $R1file .fastq.gz` 
   echo ${BUILD_DIR}/${moniker}.lanemergedbam | sed 's/_R1_/_/g' - >> $moniker_list_file 
done

for moniker in `sort -u $moniker_list_file`; do
   moniker_string="$moniker_string $moniker"
done


# get the Readgroup prefix to be passed in 
RGPREFIX='@RG\\tSM:'${SAMPLE}'\\tPL:illumina\\tLB:'${SAMPLE}'\\tCN:AgResearch'


cp ResequencingWF1.5.mk $TEMP_DIR
cd $TEMP_DIR
if [ $DRY_RUN != "no" ]; then
   echo "***** dry run only *****"
   set -x
   make -f ResequencingWF1.5.mk -d --no-builtin-rules -j 8 -n BWA_reference=$REF_GENOME dd=$DATA_DIR p1=$p1 p2=$p2 prestr=$prestr midstr=$midstr poststr=$poststr rgprefix=${RGPREFIX} mytmp=$TEMP_DIR builddir=$BUILD_DIR removeSubjectDuplicates=y removeLaneDuplicates=n  lanemergedBAMIncludeList="${moniker_string}" ${BUILD_DIR}/${SAMPLE}.all > ${SAMPLE}.log 2>&1
else
   set -x
   make -f ResequencingWF1.5.mk -d --no-builtin-rules -j 8 BWA_reference=$REF_GENOME dd=$DATA_DIR p1=$p1 p2=$p2 prestr=$prestr midstr=$midstr poststr=$poststr rgprefix=${RGPREFIX} mytmp=$TEMP_DIR builddir=$BUILD_DIR removeSubjectDuplicates=y removeLaneDuplicates=n  lanemergedBAMIncludeList="${moniker_string}" ${BUILD_DIR}/${SAMPLE}.all > ${SAMPLE}.log 2>&1
fi
set +x

if [ $? == 0 ]; then
   echo "cleaning tardis/quadtrim tempdata"
   set -x
   find ${BUILD_DIR} -name "tardis*" -type d -exec rm -rf {} \;
   set +x
else
   echo "(build appears to have failed so skipping the tardis/quadtrim clean)"
fi



# make a precis of the log file
make -i -f ResequencingWF1.5.mk ${SAMPLE}.logprecis 
# prepare a listing of software versions
make -i -f ResequencingWF1.5.mk ${SAMPLE}.versions 


cat ${SAMPLE}.logprecis

