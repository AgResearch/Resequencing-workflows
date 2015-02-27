#!/bin/sh
# example of building the target re-sequencing output files using ResequencingWF1.3.mk

ANIMAL=ferdinand
LANE_MONIKER=

BUILD_DIR=/dataset/bulls/scratch/${ANIMAL}
DATA_DIR=/dataset/bulls/archive/nzgl/Raw 
RGPREFIX='@RG\tSM:'${ANIMAL}'\tPL:illumina\tLB:'${ANIMAL}'\tCN:AgResearch'

# do a dry run first then uncomment the real run 

# dry run : 
make -f ResequencingWF1.3.mk -d --no-builtin-rules -j 4 -n dd=$DATA_DIR p1=1 p2=2 prestr=H9 midstr=_R poststr=_001 rgprefix=${RGPREFIX} mytmp=/dataset/bulls/scratch/tmp removeSampleDuplicates=y removeLaneDuplicates=n  lanemergedBAMIncludeList="${BUILD_DIR}/H9FHAADXX-1263-05-5-1_TAGCTT_L001_001.lanemergedbam ${BUILD_DIR}/H9FHAADXX-1263-05-5-1_TAGCTT_L002_001.lanemergedbam  ${BUILD_DIR}/H9WBCADXX-1263-05-5-1_TAGCTT_L002_001.lanemergedbam ${BUILD_DIR}/H9WUGADXX-1263-05-5-1_TAGCTT_L001_001.lanemergedbam ${BUILD_DIR}/H9WUGADXX-1263-05-5-1_TAGCTT_L002_001.lanemergedbam"  ${BUILD_DIR}/${ANIMAL}.all  > ${ANIMAL}.log 2>&1 

# actual run 
#make -f ResequencingWF1.3.mk -d --no-builtin-rules -j 4 dd=$DATA_DIR p1=1 p2=2 prestr=H9 midstr=_R poststr=_001 rgprefix=${RGPREFIX} mytmp=/dataset/bulls/scratch/tmp removeSampleDuplicates=y removeLaneDuplicates=n  lanemergedBAMIncludeList="${BUILD_DIR}/H9FHAADXX-1263-05-5-1_TAGCTT_L001_001.lanemergedbam ${BUILD_DIR}/H9FHAADXX-1263-05-5-1_TAGCTT_L002_001.lanemergedbam  ${BUILD_DIR}/H9WBCADXX-1263-05-5-1_TAGCTT_L002_001.lanemergedbam ${BUILD_DIR}/H9WUGADXX-1263-05-5-1_TAGCTT_L001_001.lanemergedbam ${BUILD_DIR}/H9WUGADXX-1263-05-5-1_TAGCTT_L002_001.lanemergedbam"  ${BUILD_DIR}/${ANIMAL}.all  > ${ANIMAL}.log 2>&1 
