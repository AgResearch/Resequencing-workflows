#!/bin/sh
# this version uses makefile 1.5 which uses quadtrim rather than flexbar

function get_opts() {

DRY_RUN=no
DEBUG=no
R1_FILE_PATTERN="*_L*_R1*.fastq.gz"
HPCTYPE=local
THREADS=8
INPUT_METHOD=link_files
p1=""
p2=""
prestr=""
midstr=""
poststr=""
quadtrim_option_set=sheep_set

help_text="
 examples : \n
\n
 ./ResequencingWF1.5.sh  -S NZCPWF000001391796 -D /dataset/hiseq/scratch/postprocessing_dev/150508_D00390_0226_AC6H4RANXX.processed/bcl2fastq/Project_WGS_NICHETRAITS/Sample_1016456575 -G /dataset/datacache/scratch/ensembl/oar3.1/indexes/Ovis_aries.Oar_v3.1.dna_sm.toplevel.fa -B /dataset/AG_1000_sheep/scratch/PHUAS_processing_052015 -T /dataset/AG_1000_sheep/scratch/PHUAS_processing_052015 -p 1016456575_GTGAAA -q 1 -r 2 -s _R -t _0\n
\n
./ResequencingWF1.5.sh -n -S NZCPWF100017865294 -X "*TTAGGC_L*_R1_*.fastq.gz" -D /dataset/BLGsheep/archive/NZGL01418_1/Raw  -G /dataset/datacache/scratch/ensembl/oar3.1/indexes/Ovis_aries.Oar_v3.1.dna_sm.toplevel.fa -B /dataset/AG_1000_sheep/scratch/general_processing_062015 -T /dataset/AG_1000_sheep/scratch/general_processing_062015 -p C6KFHANXX-1418-2-05-01_TTAGGC -q 1 -r 2 -s _R -t _0\n
(using the -X argument to specify a pattern to match when marshalling files from the \n
data folder - e.g. when there are a number of different barcoded outputs in a single folder) \n
\n
for more exampples see ResequencingWF1.5.examples\n

"

while getopts ":ndhe:S:B:D:G:X:T:J:C:p:q:r:s:t:I:" opt; do
  case $opt in
    n)
      DRY_RUN=yes
      ;;
    d)
      DEBUG=yes
      ;;
    h)
      echo -e $help_text
      exit 0
      ;;
    T)
      TEMP_ROOT=$OPTARG
      ;;
    J)
      THREADS=$OPTARG
      ;;
    C)
      HPCTYPE=$OPTARG
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
    e)
      quadtrim_option_set=$OPTARG
      ;;
    X)
      R1_FILE_PATTERN=$OPTARG
      ;;
    I)
      INPUT_METHOD=$OPTARG
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
  if [[ "$HPCTYPE" != "condor" && "$HPCTYPE" != "local" ]]; then
    echo "hpctype must be condor or local"
    exit 1
  fi
  if [ $INPUT_METHOD == "sra" ]; then
     if [ $R1_FILE_PATTERN == "*_L*_R1*.fastq.gz" ]; then
        R1_FILE_PATTERN="*.sra/*.sra"
     fi
  fi


}

function echo_opts() {
  echo SAMPLE=$SAMPLE
  echo DATA_DIR=$DATA_DIR
  echo BUILD_ROOT=$BUILD_ROOT
  echo LANE_TARGETS=$LANE_TARGETS
  echo DRY_RUN=$DRY_RUN
  echo DEBUG=$DEBUG
  echo HPCTYPE=$HPCTYPE
  echo THREADS=$THREADS
  echo quadtrim_option_set=$quadtrim_option_set
  echo INPUT_METHOD=$INPUT_METHOD
}

#
# edit this method to set required environment (or set up
# before running this script)
#
function configure_env() {
  if [ -d /usr/lib64/samtools0119/bin ]; then
     module load samtools/0.1.19
  else
     echo "assuming environment has path to samtools 0.1.19 or equivalent"
  fi

  if [ -f /dataset/AG_1000_bulls/active/bin/sambamba ]; then
     PATH=$PATH:/dataset/AG_1000_bulls/active/bin
  else
     echo "assuming environment has path to sambamba"
  fi
     
  if [ -f /dataset/hiseq/active/bin/quadtrim ]; then
     PATH=$PATH:/dataset/hiseq/active/bin
  else
     echo "assuming environment has path to quadtrim"
  fi

  if [ -f /dataset/AFC_dairy_cows/archive/GenomeAnalysisTKLite-2.3-9-gdcdccbb/GenomeAnalysisTKLite.jar  ]; then
     GATK_LITE_JAR=/dataset/AFC_dairy_cows/archive/GenomeAnalysisTKLite-2.3-9-gdcdccbb/GenomeAnalysisTKLite.jar
     export GATK_LITE_JAR
  else
     echo "assuming environment has GATK_LITE_JAR containing path to  GenomeAnalysisTKLite.jar"
  fi

  echo "configured the following environment variables :"
  echo "PATH=$PATH"
  echo "GATK_LITE_JAR=$GATK_LITE_JAR"

  TEMP_DIR=$TEMP_ROOT/${SAMPLE}tmp
  BUILD_DIR=$BUILD_ROOT/${SAMPLE}
  if [ ! -d $TEMP_DIR ]; then
     mkdir $TEMP_DIR
  fi

  if [ ! -d $BUILD_DIR ]; then
     mkdir $BUILD_DIR
  fi

}


function get_lane_targets() {
   if [ $INPUT_METHOD == "link_files" ]; then 
      # get the lane monikers string needed by the makefile. Basically this looks at all the input files 
      # in the data folder, and works out what targets they imply. (The make process will then work
      # backwards from those targets to in the input files in the data folder !  ) 
      set -x
      moniker_string=""
      moniker_list_file=/tmp/${$}moniker_list.tmp
      echo  > $moniker_list_file 
      R1_PATH_PATTERN=$DATA_DIR/$R1_FILE_PATTERN
      for R1file in $R1_PATH_PATTERN; do 
         moniker=`basename $R1file .fastq.gz` 
         if [ -z $poststr ]; then
            echo ${BUILD_DIR}/${moniker}.lanemergedbam | sed 's/_R1//g' - >> $moniker_list_file
         else
            echo ${BUILD_DIR}/${moniker}.lanemergedbam | sed 's/_R1_/_/g' - >> $moniker_list_file 
         fi
      done

      for moniker in `sort -u $moniker_list_file`; do
         moniker_string="$moniker_string $moniker"
      done
      set +x
   elif [ $INPUT_METHOD == "sra" ]; then
      # the input sra file(s) is/are under the data folder - each sra file contains a run of the 
      # sample , there may be one or many. For example ERS1171396 may contain one or more folders such as SRR4291010 each of these
      # will contain a file like SRR4291010 from which we can extract SRR4291010_1.fastq.gz  SRR4291010_2.fastq.gz
      # Hence - enurate all *.fastq.* which need to be extracted from a sample, and pack these into 
      # the moiker_string which becomes part of the lane target
      # 
      set -x
      moniker_string=""
      moniker_list_file=/tmp/${$}moniker_list.tmp
      R1_PATH_PATTERN=$DATA_DIR/$R1_FILE_PATTERN
      for R1file in $R1_PATH_PATTERN; do 
         moniker=`basename $R1file .sra` 
         echo ${BUILD_DIR}/${moniker}.lanemergedbam >> $moniker_list_file
      done

      for moniker in `sort -u $moniker_list_file`; do
         moniker_string="$moniker_string $moniker"
      done
      set +x
   else
      echo "unknown input method $INPUT_METHOD"
      exit 1
   fi
}


get_opts $@

check_opts

echo_opts

configure_env

get_lane_targets

##### from here inline code to run the processing


# copy or init tardis config file 
if [ -f .tardishrc ]; then
   echo "found .tardishrc in current directory will use this (copying to $TEMP_DIR)"
   cp  .tardishrc $TEMP_DIR
else
   if [ $HPCTYPE == "condor" ]; then
      echo "setting up the following tardis startup file in $TEMP_DIR :"
      echo "
[tardish]

[tardis_engine]
job_template_name=condor_send_env_job
shell_template_name=condor_shell
"
      echo "
[tardish]

[tardis_engine]
job_template_name=condor_send_env_job
shell_template_name=condor_shell
" > $TEMP_DIR/.tardishrc
   else
      echo "setting up the following tardis startup file in $TEMP_DIR :"
      echo "
[tardish]

[tardis_engine]
shell_template_name=local_shell
max_processes=$THREADS
hpctype=local
"
      echo "
[tardish]

[tardis_engine]
shell_template_name=local_shell
max_processes=$THREADS
hpctype=local
" > $TEMP_DIR/.tardishrc
   fi
fi


# get the Readgroup prefix to be passed in 
RGPREFIX='@RG\\tSM:'${SAMPLE}'\\tPL:illumina\\tLB:'${SAMPLE}'\\tCN:AgResearch'


cp ResequencingWF1.5.mk $TEMP_DIR
cp get_fastq.sh $TEMP_DIR
cd $TEMP_DIR
if [ $DRY_RUN != "no" ]; then
   echo "***** dry run only *****"
   set -x
   make -f ResequencingWF1.5.mk -d --no-builtin-rules -j $THREADS -n input_method=$INPUT_METHOD BWA_reference=$REF_GENOME quadtrim_option_set=$quadtrim_option_set dd=$DATA_DIR p1=$p1 p2=$p2 prestr=$prestr midstr=$midstr poststr=$poststr rgprefix=${RGPREFIX} mytmp=$TEMP_DIR builddir=$BUILD_DIR removeSubjectDuplicates=y removeLaneDuplicates=n  lanemergedBAMIncludeList="${moniker_string}" ${BUILD_DIR}/${SAMPLE}.all > ${SAMPLE}.log 2>&1
else
   set -x
   make -f ResequencingWF1.5.mk -d --no-builtin-rules -j $THREADS input_method=$INPUT_METHOD BWA_reference=$REF_GENOME quadtrim_option_set=$quadtrim_option_set dd=$DATA_DIR p1=$p1 p2=$p2 prestr=$prestr midstr=$midstr poststr=$poststr rgprefix=${RGPREFIX} mytmp=$TEMP_DIR builddir=$BUILD_DIR removeSubjectDuplicates=y removeLaneDuplicates=n  lanemergedBAMIncludeList="${moniker_string}" ${BUILD_DIR}/${SAMPLE}.all > ${SAMPLE}.log 2>&1
fi
set +x

if [[ ( $DEBUG != "yes" )  && ( $DRY_RUN != "yes" )  ]]; then
   if [ $? == 0 ]; then
      echo "cleaning tardis/quadtrim tempdata"
      set -x
      find ${BUILD_DIR} -name "tardis*" -type d -exec rm -rf {} \;
      set +x
   else
      echo "(build appears to have failed so skipping the tardis/quadtrim clean)"
   fi
fi

# make a precis of the log file
make -i -f ResequencingWF1.5.mk ${SAMPLE}.logprecis 
# prepare a listing of software versions
make -i -f ResequencingWF1.5.mk ${SAMPLE}.versions 


cat ${SAMPLE}.logprecis

# suggest an archiving step.
echo "suggest you check the run then archive it. Example command to archive : 
./ArchiveBuild1.5.sh   -S $SAMPLE -A /dataset/AG_1000_sheep/archive/2015_processing_results -B /dataset/AG_1000_sheep/scratch/general_processing_062015  -T /dataset/AG_1000_sheep/scratch/general_processing_062015
"

