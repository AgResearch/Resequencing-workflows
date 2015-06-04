#!/bin/sh
#
# build a resequencing project - i.e. align lanes of sequencing and call snps
#
# !!!!!!!!!!!!! this was nto actually tested - abandoned and went to 1.4 !!!!!!!!!!!!!!

function get_opts() {
DRY_RUN=no
help_text="
 examples : \n
 ./ResequencingWF1.3.sh  -n -A ferdinand -D /dataset/bulls/archive/nzgl/Raw -B /dataset/bulls/scratch -L  \"H9FHAADXX-1263-05-5-1_TAGCTT_L001_001,H9FHAADXX-1263-05-5-1_TAGCTT_L002_001,H9WBCADXX-1263-05-5-1_TAGCTT_L002_001,H9WUGADXX-1263-05-5-1_TAGCTT_L001_001,H9WUGADXX-1263-05-5-1_TAGCTT_L002_001\" \n
"

while getopts ":nhA:B:D:L:" opt; do
  case $opt in
    n)
      DRY_RUN=yes
      ;;
    h)
      echo -e $help_text
      exit 0
      ;;
    A)
      ANIMAL=$OPTARG
      ;;
    B)
      BUILD_ROOT=$OPTARG
      ;;
    D)
      DATA_DIR=$OPTARG
      ;;
    L)
      LANE_TARGETS=$OPTARG
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
  if [ -z "$ANIMAL" ]; then
    echo "must specify animal using -A option"
    exit 1
  fi
  if [ -z "$LANE_TARGETS" ]; then
    echo "must specify lane targets using -L option"
    exit 1
  fi

}

function echo_opts() {
  echo ANIMAL=$ANIMAL
  echo DATA_DIR=$DATA_DIR
  echo BUILD_ROOT=$BUILD_ROOT
  echo LANE_TARGETS=$LANE_TARGETS
  echo DRY_RUN=$DRY_RUN
}


get_opts $@

check_opts

echo_opts

# from here in-line code to do run
BUILD_DIR=$BUILD_ROOT/${ANIMAL}
RGPREFIX='@RG\tSM:'${ANIMAL}'\tPL:illumina\tLB:'${ANIMAL}'\tCN:AgResearch'

# patch lane targets  - they are passed in comma seperated - change to
# blank seperated and prefix with build root
LANE_TARGETS1=`echo $LANE_TARGETS | sed 's/,/ /g' -`
LANE_TARGETS=""
for target in $LANE_TARGETS1; do
   LANE_TARGETS="$BUILD_ROOT/$target $LANE_TARGETS"
done


if [ $DRY_RUN == "no" ]; then
  echo make -f ResequencingWF1.3.mk -d --no-builtin-rules -j 4 dd=$DATA_DIR p1=1 p2=2 prestr=H9 midstr=_R poststr=_001 rgprefix=${RGPREFIX} mytmp=/dataset/bulls/scratch/tmp removeSampleDuplicates=y removeLaneDuplicates=n  lanemergedBAMIncludeList="$LANE_TARGETS"  ${BUILD_DIR}/${ANIMAL}.all  > ${ANIMAL}.log 2>&1
else
  echo make -f ResequencingWF1.3.mk -d --no-builtin-rules -j 4 -n dd=$DATA_DIR p1=1 p2=2 prestr=H9 midstr=_R poststr=_001 rgprefix=${RGPREFIX} mytmp=/dataset/bulls/scratch/tmp removeSampleDuplicates=y removeLaneDuplicates=n  lanemergedBAMIncludeList="$LANE_TARGETS"  ${BUILD_DIR}/${ANIMAL}.all  > ${ANIMAL}.log 2>&1
fi
