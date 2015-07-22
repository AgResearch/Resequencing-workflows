#!/bin/sh
# driver script for an associated makefile, used to move a completed project to
# an archive filesystem and clean up intermediate files 

function get_opts() {

DRY_RUN=no
help_text="
 examples : \n
 ./ArchiveBuild1.5.sh  -n -S NZCPWF000001391796 -A /dataset/AG_1000_sheep/archive/2015_processing_results -B /dataset/AG_1000_sheep/scratch/PHUAS_processing_052015 -T /dataset/AG_1000_sheep/scratch/PHUAS_processing_052015 \n
"

while getopts ":nhS:B:D:G:T:A:" opt; do
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
    A)
      ARCHIVE_ROOT=$OPTARG
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
  echo BUILD_ROOT=$BUILD_ROOT
  echo DRY_RUN=$DRY_RUN
}


get_opts $@

check_opts

echo_opts

##### from here inline code to run the processing
TEMP_DIR=$TEMP_ROOT/${SAMPLE}tmp
BUILD_DIR=$BUILD_ROOT/${SAMPLE}
ARCHIVE_DIR=$ARCHIVE_ROOT/${SAMPLE}
if [ ! -d $TEMP_DIR ]; then
   mkdir $TEMP_DIR
fi

if [ ! -d $BUILD_DIR ]; then
   mkdir $BUILD_DIR
fi

if [ ! -d $ARCHIVE_DIR ]; then
   mkdir $ARCHIVE_DIR
fi

cp ArchiveBuild1.5.mk $TEMP_DIR
cd $TEMP_DIR
echo "(logging run to $TEMP_DIR/${SAMPLE}.archivelog)" 
if [ $DRY_RUN != "no" ]; then
   echo "***** dry run only *****"
   set -x
   make -f ArchiveBuild1.5.mk -d --no-builtin-rules -j 8 -n builddir=$BUILD_DIR archivedir=${ARCHIVE_DIR} tempdir=$TEMP_DIR ${ARCHIVE_DIR}/${SAMPLE}.all > ${SAMPLE}.archivelog 2>&1
else
   set -x
   make -f ArchiveBuild1.5.mk -d --no-builtin-rules -j 8 builddir=$BUILD_DIR archivedir=${ARCHIVE_DIR} tempdir=$TEMP_DIR ${ARCHIVE_DIR}/${SAMPLE}.all > ${SAMPLE}.archivelog 2>&1
fi
set +x

# make a precis of the log file
make -i -f ArchiveBuild1.5.mk ${SAMPLE}.archivelogprecis > /dev/null  2>&1

cat ${SAMPLE}.archivelogprecis

