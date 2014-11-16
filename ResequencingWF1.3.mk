#
# Resequencing workflow version 3. 
#***************************************************************************************
# changes 
#***************************************************************************************
# 31/7/2014 flexbar run :  --min-read-length 35 should have been --min-read-length 50
# 31/7/2014 added a new fastqc target , doing fastqc on the flexbar output 
# 
#***************************************************************************************
# references:
#***************************************************************************************
# make: 
#     http://www.gnu.org/software/make/manual/make.html
# gatk:
#     http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html
#     http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_coverage_DepthOfCoverage.html
#     https://www.broadinstitute.org/gatk/guide/tagged?tag=bwa
# SAM
#     http://samtools.github.io/hts-specs/SAMv1.pdf
#
# Note - you should avoid making your targets in the same folder as your source data (in case you overwrite or clean up
# something important). Create an empty build folder, and then your target to make is empty_build_folder/animal_name.all    
#
#*****************************************************************************************
# examples  
#*****************************************************************************************
#*****************
# (this example uses test data and builds output for an animal JJLNZLM000000000394 sequenced across two lanes. It uses 
# hard coded "dummy" targetIntervals, for testing (real intervals are slow)
#*****************
#make -f ResequencingWF1.3.mk -d --no-builtin-rules dd=/dataset/AG_1000_bulls/active/testdata p1=1 p2=2 prestr=C3PACACXX-1143- midstr=_R poststr=_001 rgprefix='@RG\tSM:JJLNZLM000000000394\tPL:illumina\tLB:JJLNZLM000000000394\tCN:AgResearch' mytmp=/dataset/AG_1000_bulls/scratch/tmp removeSampleDuplicates=y removeLaneDuplicates=n targetIntervals='Chr1:100-200' lanemergedBAMIncludeList="/dataset/AG_1000_bulls/scratch/mktest3/C3PACACXX-1143-01-5-1_GCCAAT_L001_001.lanemergedbam /dataset/AG_1000_bulls/scratch/mktest3/C3PACACXX-1143-01-5-1_GCCAAT_L002_001.lanemergedbam"  /dataset/AG_1000_bulls/scratch/mktest3/JJLNZLM000000000394.all
#
#******************
# (this example is from an actual run and uses shell variable settings for the animalid and datafilename to 
# improve readability  - i.e. these commands would be executed in a batch file)
#*****************
#
#ANIMAL=AANNZLM017568091689
#LANE_MONIKER=C3PF2ACXX-1143-12-5-1_CAGATC_L008_001

#BUILD_DIR=/dataset/AG_1000_bulls/scratch/${ANIMAL}
#DATA_DIR=/dataset/AG_1000_bulls/archive/nzgl01143/NZGL01143_C3PF2ACXX/Raw
#RGPREFIX='@RG\tSM:'${ANIMAL}'\tPL:illumina\tLB:'${ANIMAL}'\tCN:AgResearch'

#make -f ResequencingWF1.3.mk -d --no-builtin-rules -j 4 dd=$DATA_DIR p1=1 p2=2 prestr=C3PF2ACXX-1143- midstr=_R poststr=_001 rgprefix=${RGPREFIX} mytmp=/dataset/AG_1000_bulls/scratch/tmp removeSampleDuplicates=y removeLaneDuplicates=n  lanemergedBAMIncludeList="${BUILD_DIR}/${LANE_MONIKER}.lanemergedbam"  ${BUILD_DIR}/${ANIMAL}.all  > ${ANIMAL}.log 2>&1
#
#
#*****************
# (this example is from an actual run for a bull that was run on multiple lanes)
# (also used shell variables, but the animal and rgprefix were hard-coded)
#*****************
#BUILD_DIR=/dataset/AG_1000_bulls/scratch/JJLNZLM000000000394
#DATA_DIR=/dataset/AG_1000_bulls/archive/nzgl01143/NZGL01143_C3PACACXX/Raw
#MYTMP=/dataset/AG_1000_bulls/active/tmp

#make -f ResequencingWF1.3.mk -d --no-builtin-rules -j 4 dd=$DATA_DIR p1=1 p2=2 prestr=C3PACACXX-1143- midstr=_R poststr=_001 rgprefix='@RG\tSM:JJLNZLM000000000394\tPL:illumina\tLB:JJLNZLM000000000394\tCN:AgResearch' mytmp=${MYTMP} removeSampleDuplicates=y removeLaneDuplicates=n  lanemergedBAMIncludeList="${BUILD_DIR}/C3PACACXX-1143-01-5-1_GCCAAT_L001_001.lanemergedbam ${BUILD_DIR}/C3PACACXX-1143-01-5-1_GCCAAT_L002_001.lanemergedbam ${BUILD_DIR}/C3PACACXX-1143-01-5-1_GCCAAT_L003_001.lanemergedbam ${BUILD_DIR}/C3PACACXX-1143-01-5-1_GCCAAT_L004_001.lanemergedbam ${BUILD_DIR}/C3PACACXX-1143-01-5-1_GCCAAT_L005_001.lanemergedbam ${BUILD_DIR}/C3PACACXX-1143-01-5-1_GCCAAT_L006_001.lanemergedbam"  ${BUILD_DIR}/JJLNZLM000000000394.all  > JJLNZLM000000000394.log 2>&1
#
#
# ******************************************************************************************
# initialise project specific variables - these are usually overridden by command-line settings as above
# ******************************************************************************************
rgvarname=rgprefix
targetIntervalsvarname=targetIntervals
dd=/does_not_exist
p1=1
p2=2
prestr=
midstr=
poststr=
mytmp=/tmp
adaptersFile=/dataset/AG_1000_bulls/active/adapters.txt

#
#*******************************************************************************************
# project filename construction support examples
#*******************************************************************************************
#  - i.e. the variables p1,p2,prestr, midstr, poststr : these allow the rules to construct and refer to the filenames 
# implied by the target name 
# example:
# the input paired-end filenames are C3PACACXX-1143-04-5-1_TTAGGC_L008_R1_001.fastq.gz and C3PACACXX-1143-04-5-1_TTAGGC_L008_R2_001.fastq.gz
#
# These names are parsed (for the purposes of this makefile) as
# "prestring" + "stem" + "midstring" + "pair identifier" + "poststring" + ".fastq.gz"
#
# e.g. in this case we decide to parse as 
#
# "C3PACACXX-1143-" + "04-5-1_TTAGGC_L008"  + "_R" +  1 + "_001"  + ".fastq.gz"
#
# - here you would make (e.g.) C3PACACXX-1143-04-5-1_TTAGGC_L008_001.lanemergedbam,  using 
#
# p1=1 p2=2 prestr=C3PACACXX-1143- midstr=_R poststr=_001
#
# so that % will then match  a stem 04-5-1_TTAGGC_L008 in this case. The source and target filenames are
# then built back up from these.
# (You could use other values of these variables to build the same target - however the above choices mean
# we can use the same settings to build all files in a folder that have the pattern
# C3PACACXX-1143-*_R[1|2]_001). Whereas if we used (say) 
# p1=1 p2=2 prestr=C3PACACXX-1143-04-5-1_TTAGGC_ midstr=_R poststr=_001
# - so that the stem matches the lane number - this will not work for files that 
# have a different barcode or date string. (The stem that is matched is only used 
# internally - so it doesn't matter particularly what is matched)
# 
# 

# ******************************************************************************************
# other variables (not project specific)
# ******************************************************************************************
#RUN_TARDIS=echo tardis.py # this is for a dry run. Change to next lines for the real run
#RUN_FASTQC=echo fastqc # this is for a dry run. Change to next lines for the real run
#
RUN_TARDIS=tardis.py
RUN_FASTQC=fastqc
RUN_SAMBAMBA=/dataset/AG_1000_bulls/active/bin/sambamba
RUN_JAVA=/usr/bin/java
GATK=/dataset/AFC_dairy_cows/archive/GenomeAnalysisTKLite-2.3-9-gdcdccbb/GenomeAnalysisTKLite.jar

# variables for tardis and other apps
TARDIS_chunksize=300000
#TARDIS_chunksize=1000 for testing small files
TARDIS_workdir=$(mytmp)
BWA_reference=/dataset/AFC_dairy_cows/active/1000_bulls/umd_3_1_reference_1000_bull_genomes.fa


#*****************************************************************************************
# from here are the targets and rules
#*****************************************************************************************


###############################################
# top level phony test targets  (used for various testing / debugging)
###############################################
#.PHONY : $(prestr)%$(poststr).test
#$(prestr)%$(poststr).test: $(prestr)%$(poststr).vcf $(prestr)%$(poststr).bamstats $(prestr)%$(poststr).coverage.sample_summary
#	echo "making test target"
.PHONY : %.test
%.test:  $(lanemergedBAMIncludeList)
	echo "making test target..."
	echo $(lanemergedBAMIncludeList) 

.PHONY : %.testdependency 
$(prestr)%$(poststr).testdependency:
	echo making $*  
	# make the complete read group
        # for example the complete readgroup is '@RG\tID:C3PACACXX-1143-01-5-1_GCCAAT_L001\tSM:JJLNZLM000000000394\tPL:illumina\tLB:JJLNZLM000000000394\tPU:C3PACACXX-1143-01-5-1_GCCAAT_L001\tCN:AgResearch'
        # we are given part of this as the variable rgprefix :
        # '@RG\tSM:JJLNZLM000000000394\tPL:illumina\tLB:JJLNZLM000000000394\tCN:AgResearch' 
        # and we will append \tID:C3PACACXX-1143-01-5-1_GCCAAT_L001\tPU:C3PACACXX-1143-01-5-1_GCCAAT_L001'
        # 
	echo will use readgroup info
	echo '$(rgprefix)\tID:$(prestr)$(*F)$(poststr)\tPU:$(prestr)$(*F)$(poststr)'

###############################################
# top level target "logprecis" . This extracts from the log file the 
# relevant commands that were run, in a readable order 
# - run this after the actual build has been completed
###############################################
%.logprecis: %.log
	echo "fastqc" > $*.logprecis
	echo "------" >> $*.logprecis
	egrep "^$(RUN_FASTQC)" $*.log >> $*.logprecis
	echo "flexbar" >> $*.logprecis
	echo "-------" >> $*.logprecis
	egrep "^$(RUN_TARDIS) -w" $*.log  >> $*.logprecis
	echo "bwa" >> $*.logprecis
	echo "---" >> $*.logprecis
	egrep "^bwa" $*.log  >> $*.logprecis
	echo "bamstats" >> $*.logprecis
	echo "--------" >> $*.logprecis
	egrep "bamtools stats" $*.log  >> $*.logprecis
	echo "bam filtering" >> $*.logprecis
	echo "-------------" >> $*.logprecis
	egrep "^$(RUN_SAMBAMBA) view" $*.log >> $*.logprecis
	echo "bam sorting" >> $*.logprecis
	echo "-----------" >> $*.logprecis
	egrep "^$(RUN_SAMBAMBA) sort" $*.log >> $*.logprecis
	echo "bam deduplication" >> $*.logprecis
	echo "-----------------" >> $*.logprecis
	egrep "^$(RUN_SAMBAMBA) markdup" $*.log >> $*.logprecis
	echo "bam merging" >> $*.logprecis
	echo "-----------" >> $*.logprecis
	egrep "^$(RUN_SAMBAMBA) merge" $*.log >> $*.logprecis
	echo "gatk : intervals" >> $*.logprecis	
	echo "----------------" >> $*.logprecis
	egrep "^$(RUN_JAVA)" $*.log | grep RealignerTargetCreator >> $*.logprecis
	echo "gatk : realigner" >> $*.logprecis	
	echo "----------------" >> $*.logprecis
	egrep "^$(RUN_JAVA)" $*.log | grep IndelRealigner >> $*.logprecis
	echo "vcf" >> $*.logprecis	
	echo "---" >> $*.logprecis
	egrep "^samtools mpileup" $*.log >> $*.logprecis
	echo "gatk : coverage" >> $*.logprecis	
	echo "---------------" >> $*.logprecis
	egrep "^$(RUN_JAVA)" $*.log | grep DepthOfCoverage >> $*.logprecis


###############################################
# top level phony target "versions"  - output versions of all tools 
# - note , you need to tell make to ignore errors 
# for this target - for some tools , the only way to get a version
# string is to run them with options that provoke an error
# Where useful, this reports the package name as well 
# as getting the tool to report its version.
# (in some cases e.g. bamtools this is all the version
# info there is as the tool itself doesn't report any
# version info)
# (not all the tools that were used are installed 
# as packages currently )
###############################################
.PHONY : versions.log
versions.log:
	echo "Tool versions : " > versions.log
	echo "Java" >> versions.log
	echo "----"  >> versions.log
	echo $(RUN_JAVA) -version  >> versions.log
	$(RUN_JAVA) -version  >> versions.log 2>&1
	echo "fastqc"  >> versions.log
	echo "------"  >> versions.log
	echo $(RUN_FASTQC) -version  >> versions.log
	$(RUN_FASTQC) -version  >> versions.log  2>&1
	echo rpm -q fastqc >> versions.log
	rpm -q fastqc >> versions.log 2>&1
	echo "sambamba"  >> versions.log
	echo "--------"  >> versions.log
	echo $(RUN_SAMBAMBA) -h  >> versions.log
	$(RUN_SAMBAMBA) -h >> versions.log  2>&1
	echo "gatk"  >> versions.log
	echo "----"  >> versions.log
	echo $(RUN_JAVA) -Xmx4G -jar -Djava.io.tmpdir=$(mytmp) -jar $(GATK) >> versions.log
	$(RUN_JAVA) -Xmx4G -jar -Djava.io.tmpdir=$(mytmp) -jar $(GATK) -h | head -5 >> versions.log  2>&1
	echo "bwa"  >> versions.log
	echo "---"  >> versions.log
	echo bwa >> versions.log
	bwa >> versions.log  2>&1
	echo rpm -q bwa >> versions.log
	rpm -q bwa >> versions.log 2>&1
	echo "samtools"  >> versions.log
	echo "--------"  >> versions.log
	echo samtools >> versions.log
	samtools >> versions.log  2>&1
	echo rpm -q samtools >> versions.log
	rpm -q samtools >> versions.log 2>&1
	echo "bamtools"  >> versions.log
	echo "--------"  >> versions.log
	echo bamtools >> versions.log
	bamtools >> versions.log  2>&1
	echo rpm -q bamtools >> versions.log
	rpm -q bamtools >> versions.log 2>&1
	echo "flexbar"  >> versions.log
	echo "-------"  >> versions.log
	echo flexbar -version >> versions.log
	flexbar -version >> versions.log  2>&1
	echo rpm -q flexbar >> versions.log
	rpm -q flexbar >> versions.log 2>&1


###############################################
# top level phony target "all"  - this is the one thats usually used 
###############################################
.PHONY : %.all 
#%.all: %.vcf     # use this if you don't want to bother with coverage (- its slow) 
%.all: %.vcf %.coverage.sample_summary 
	echo "making all targets"


###############################################
# how to make coverage output files 
# (This makes a number of files like *.coverage.* 
# as well as the coverage.sample_summary file)
###############################################
%.coverage.sample_summary: %.realignedbam
	$(RUN_JAVA) -Xmx4G -jar -Djava.io.tmpdir=$(mytmp) -jar $(GATK) -T DepthOfCoverage -R $(BWA_reference) -I $*_realigned.bam --omitDepthOutputAtEachBase --logging_level ERROR --summaryCoverageThreshold 10 --summaryCoverageThreshold 20 --summaryCoverageThreshold 30 --summaryCoverageThreshold 40 --summaryCoverageThreshold 50 --summaryCoverageThreshold 80 --summaryCoverageThreshold 90 --summaryCoverageThreshold 100 --summaryCoverageThreshold 150 --minBaseQuality 15 --minMappingQuality 30 --start 1 --stop 1000 --nBins 999 -dt NONE -o $*.coverage

###############################################
# how to make the vcf
###############################################
%.vcf: %.realignedbam
	samtools mpileup -ugf $(BWA_reference) $*_realigned.bam | bcftools view -N -cvg - > $*.vcf 


###############################################
# how to make the "full" vcf - output everything not just variant positions 
# example : make -f ResequencingWF1.3.mk -d --no-builtin-rules /dataset/AG_1000_bulls/scratch/JJLNZLM000000000402/JJLNZLM000000000402.fullvcf
# (consider adding fullvcf as a dependency of "all", or else just change the standard 
# vcf build so it includes all positions ?)
###############################################
%.fullvcf: %.realignedbam
	samtools mpileup -ugf $(BWA_reference) $*_realigned.bam | bcftools view -N -cg - > $*.fullvcf



###############################################
# how to make a patched realignedbam file. This was needed because the flexbar processing used min length of 35
# rather than 50. This assumes that the build folder has been moved to the archive filesystem , and that 
# a new empty build folder has been created, containing (a shortcut to ) the original realignedbam that 
# is to be patched. 
###############################################
%.patchedrealignedbam: %.realignedbam
	#patch it (the calling script has set up a soft link from the build folder back to the realigned bam)
	$(RUN_SAMBAMBA) view -t 8 -h -f bam -o $*.patchedrealignedbam -F "sequence_length >=50" $*.realignedbam 

	#unlink the original realigned bam 
	unlink $*.realignedbam 
    
	#link to the new one and index it - iso after the build of this target is complete , we can make the orginal top level target 
        # - i.e. rebuild the vcf and coverage
	ln -fs $*.patchedrealignedbam $*_realigned.bam
	ln -fs $*.patchedrealignedbam $*.realignedbam
	bamtools index -in $*_realigned.bam



###############################################
# how to make the "realigned" bam 
# (if an intervals argument is specified then we do not 
# depend on building the intervals target. For example you could specify
# targetIntervals Chr1:100-200 for testing (a run which builds the  
# intervals is very slow even for a small input file - so for testing need to hard-code intervals)
###############################################
ifndef $(targetIntervalsvarname)
%.realignedbam: %.samplemergedbam %.intervals
	# if we have samplemergedbam then we should also have _samplemerged.bam plus index
	#$(RUN_JAVA) -Xmx10g -Djava.io.tmpdir=$(mytmp) -jar $(GATK) -T IndelRealigner -R $(BWA_reference) -targetIntervals $*.intervals  -I $*_samplemerged.bam  -o $*.realignedbam -S LENIENT
	$(RUN_JAVA) -Xmx10g -Djava.io.tmpdir=$(mytmp) -jar $(GATK) -T IndelRealigner --read_filter NotPrimaryAlignment -R $(BWA_reference) -targetIntervals $*.intervals  -I $*_samplemerged.bam  -o $*.realignedbam -S LENIENT
else
%.realignedbam: %.samplemergedbam 
	#$(RUN_JAVA) -Xmx10g -Djava.io.tmpdir=$(mytmp) -jar $(GATK) -T IndelRealigner -R $(BWA_reference) -targetIntervals $(targetIntervals) -I $*_samplemerged.bam -o $*.realignedbam -S LENIENT
	$(RUN_JAVA) -Xmx10g -Djava.io.tmpdir=$(mytmp) -jar $(GATK) -T IndelRealigner --read_filter NotPrimaryAlignment -R $(BWA_reference) -targetIntervals $(targetIntervals) -I $*_samplemerged.bam -o $*.realignedbam -S LENIENT
endif
	# make a shortcut with suffix "bam" - many of the tools do not like custom suffixes like realignedbam etc
	ln -fs $*.realignedbam $*_realigned.bam 
	bamtools index -in $*_realigned.bam 

###############################################
# how to make the "intervals" files
###############################################
%.intervals: %.samplemergedbam
	# if we have samplemergedbam then we should also have _samplemerged.bam plus index
	$(RUN_JAVA) -Xmx4g -Djava.io.tmpdir=$(mytmp) -jar $(GATK) -T RealignerTargetCreator -nt 16 -R $(BWA_reference) -I $*_samplemerged.bam -o $*.intervals -S LENIENT -rbs 5000000


#############################################################################
# how to make the sample-merged BAM (optionally including removal of duplicates at the sample as opposed to lane level)
#############################################################################
%.samplemergedbam: $(lanemergedBAMIncludeList)
	$(RUN_SAMBAMBA) merge -t 8  $*.samplemergedbam_with_duplicates $(lanemergedBAMIncludeList)
	# if removeSubjectDuplicates = n , make the target (and the associated links with .bam suffix) simply by linking to the above.
	# (GATK only likes .bam or .sam suffices)
ifeq ($(subst N,n,$(strip $(removeSubjectDuplicates))),n)
	ln -fs $*.samplemergedbam_with_duplicates $*.samplemergedbam
	ln -fs $*.samplemergedbam_with_duplicates $*_samplemerged.bam
else
	# or, if required to removeSubjectDuplicates - process the with-duplicates merged file to remove duplicates, and link to that
	$(RUN_SAMBAMBA) markdup -t 8 -r --tmpdir=$(mytmp) --overflow-list-size=800000 --io-buffer-size=256 $*.samplemergedbam_with_duplicates $*.samplemergedbam_no_duplicates
	ln -fs $*.samplemergedbam_no_duplicates $*.samplemergedbam
	ln -fs $*.samplemergedbam_no_duplicates $*_samplemerged.bam
endif
	bamtools index -in $*_samplemerged.bam 


#############################################################################
# how to make the lane-merged BAM (optionally including removal of duplicates at this level)
#############################################################################
$(prestr)%$(poststr).lanemergedbam: $(prestr)%$(poststr).sortedbam
	$(RUN_SAMBAMBA) merge -t 8  $(*D)/$(prestr)$(*F)$(poststr).lanemergedbam_with_duplicates $(*D)/$(prestr)$(*F)$(poststr).sortedbam $(*D)/$(prestr)$(*F)$(midstr)$(p1)$(poststr)_single.sortedbam $(*D)/$(prestr)$(*F)$(midstr)$(p2)$(poststr)_single.sortedbam
	# if removeLaneDuplicates = n , make the target (and the associated links with .bam suffix) simply by linking to the above.
	# (GATK only likes .bam or .sam suffices)
ifeq ($(subst N,n,$(strip $(removeLaneDuplicates))),n)
	ln -fs $(*D)/$(prestr)$(*F)$(poststr).lanemergedbam_with_duplicates $(*D)/$(prestr)$(*F)$(poststr).lanemergedbam
	ln -fs $(*D)/$(prestr)$(*F)$(poststr).lanemergedbam_with_duplicates $(*D)/$(prestr)$(*F)$(poststr)_merged.bam
else
	# or, if required to removeLaneDuplicates - process the with-duplicates merged file to remove duplicates, and link to that
	$(RUN_SAMBAMBA) markdup -t 8 -r --tmpdir=$(mytmp) --overflow-list-size=800000 --io-buffer-size=256 $(*D)/$(prestr)$(*F)$(poststr).lanemergedbam_with_duplicates $(*D)/$(prestr)$(*F)$(poststr).lanemergedbam_no_duplicates
	ln -fs $(*D)/$(prestr)$(*F)$(poststr).lanemergedbam_no_duplicates $(*D)/$(prestr)$(*F)$(poststr).lanemergedbam
	ln -fs $(*D)/$(prestr)$(*F)$(poststr).lanemergedbam_no_duplicates $(*D)/$(prestr)$(*F)$(poststr)_merged.bam
endif
	bamtools index -in $(*D)/$(prestr)$(*F)$(poststr)_merged.bam



#############################################################################
# how to make the merged BAM (optionally including removal of duplicates)
#############################################################################
$(prestr)%$(poststr).lanemergedbam: $(prestr)%$(poststr).sortedbam
	$(RUN_SAMBAMBA) merge -t 8  $(*D)/$(prestr)$(*F)$(poststr).lanemergedbam_with_duplicates $(*D)/$(prestr)$(*F)$(poststr).sortedbam $(*D)/$(prestr)$(*F)$(midstr)$(p1)$(poststr)_single.sortedbam $(*D)/$(prestr)$(*F)$(midstr)$(p2)$(poststr)_single.sortedbam
	# if removeLaneDuplicates = n , make the target (and the associated links with .bam suffix) simply by linking to the above. 
	# (GATK only likes .bam or .sam suffices)  
ifeq ($(subst N,n,$(strip $(removeLaneDuplicates))),n)
	ln -fs $(*D)/$(prestr)$(*F)$(poststr).lanemergedbam_with_duplicates $(*D)/$(prestr)$(*F)$(poststr).lanemergedbam 
	ln -fs $(*D)/$(prestr)$(*F)$(poststr).lanemergedbam_with_duplicates $(*D)/$(prestr)$(*F)$(poststr)_merged.bam
else
	# or, if required to removeLaneDuplicates - process the with-duplicates merged file to remove duplicates, and link to that
	$(RUN_SAMBAMBA)	markdup -t 8 -r --tmpdir=$(mytmp) --overflow-list-size=800000 --io-buffer-size=256 $(*D)/$(prestr)$(*F)$(poststr).lanemergedbam_with_duplicates $(*D)/$(prestr)$(*F)$(poststr).lanemergedbam_no_duplicates
	ln -fs $(*D)/$(prestr)$(*F)$(poststr).lanemergedbam_no_duplicates $(*D)/$(prestr)$(*F)$(poststr).lanemergedbam 
	ln -fs $(*D)/$(prestr)$(*F)$(poststr).lanemergedbam_no_duplicates $(*D)/$(prestr)$(*F)$(poststr)_merged.bam
endif
	bamtools index -in $(*D)/$(prestr)$(*F)$(poststr)_merged.bam



###############################################
# how to make the sorted BAMs
###############################################
$(prestr)%$(poststr).sortedbam: $(prestr)%$(poststr).filteredbam
	$(RUN_SAMBAMBA) sort --tmpdir=$(mytmp) -t 8 -m 10G $(*D)/$(prestr)$(*F)$(poststr).filteredbam  -o $(*D)/$(prestr)$(*F)$(poststr).sortedbam
	$(RUN_SAMBAMBA) sort --tmpdir=$(mytmp) -t 8 -m 10G $(*D)/$(prestr)$(*F)$(midstr)$(p1)$(poststr)_single.filteredbam -o $(*D)/$(prestr)$(*F)$(midstr)$(p1)$(poststr)_single.sortedbam 
	$(RUN_SAMBAMBA) sort --tmpdir=$(mytmp) -t 8 -m 10G $(*D)/$(prestr)$(*F)$(midstr)$(p2)$(poststr)_single.filteredbam -o $(*D)/$(prestr)$(*F)$(midstr)$(p2)$(poststr)_single.sortedbam 

###############################################
# how to make the filtered  BAMs. We insert a dependency on bamstats 
# at this point to get those done as well 
###############################################
$(prestr)%$(poststr).filteredbam: $(prestr)%$(poststr).bam  $(prestr)%$(poststr).bamstats
	$(RUN_SAMBAMBA) view -t8 -f bam -F "mapping_quality >=20" $(*D)/$(prestr)$(*F)$(poststr).bam  -o $(*D)/$(prestr)$(*F)$(poststr).filteredbam 
	$(RUN_SAMBAMBA) view -t8 -f bam -F "mapping_quality >=20" $(*D)/$(prestr)$(*F)$(midstr)$(p1)$(poststr)_single.bam -o  $(*D)/$(prestr)$(*F)$(midstr)$(p1)$(poststr)_single.filteredbam 
	$(RUN_SAMBAMBA) view -t8 -f bam -F "mapping_quality >=20" $(*D)/$(prestr)$(*F)$(midstr)$(p2)$(poststr)_single.bam -o  $(*D)/$(prestr)$(*F)$(midstr)$(p2)$(poststr)_single.filteredbam 


###############################################
# how to get stats for the bam files
###############################################
$(prestr)%$(poststr).bamstats: $(prestr)%$(poststr).bam
	echo "===================================" > $(*D)/$(prestr)$(*F)$(poststr).bamstats 
	echo $(*D)/$(prestr)$(*F)$(poststr).bam  >> $(*D)/$(prestr)$(*F)$(poststr).bamstats
	echo "===================================" >> $(*D)/$(prestr)$(*F)$(poststr).bamstats
	bamtools stats -in $(*D)/$(prestr)$(*F)$(poststr).bam  -insert >>  $(*D)/$(prestr)$(*F)$(poststr).bamstats
	echo "===================================" >> $(*D)/$(prestr)$(*F)$(poststr).bamstats
	echo  $(*D)/$(prestr)$(*F)$(midstr)$(p1)$(poststr)_single.bam >> $(*D)/$(prestr)$(*F)$(poststr).bamstats
	echo "===================================" >> $(*D)/$(prestr)$(*F)$(poststr).bamstats
	bamtools stats -in  $(*D)/$(prestr)$(*F)$(midstr)$(p1)$(poststr)_single.bam  >> $(*D)/$(prestr)$(*F)$(poststr).bamstats
	echo "===================================" >> $(*D)/$(prestr)$(*F)$(poststr).bamstats
	echo $(*D)/$(prestr)$(*F)$(midstr)$(p2)$(poststr)_single.bam  >> $(*D)/$(prestr)$(*F)$(poststr).bamstats
	echo "===================================" >> $(*D)/$(prestr)$(*F)$(poststr).bamstats
	bamtools stats -in $(*D)/$(prestr)$(*F)$(midstr)$(p2)$(poststr)_single.bam  >> $(*D)/$(prestr)$(*F)$(poststr).bamstats


###############################################
# how to make  bam files. Made from 
# flexbar output from both pairs, and indirectly also 
# FastQC output 
# Note that we only mention 1 of the pairs in the rule - but process both
# We also optionally inject readgroup info at this point
###############################################
$(prestr)%$(poststr).bam: $(prestr)%$(midstr)$(p1)$(poststr)_fastqc.zip $(prestr)%$(midstr)$(p1)$(poststr).fastq.flex.gz 
ifndef $(rgvarname) 
	bwa mem -t 8 -M $(BWA_reference) $(*D)/$(prestr)$(*F)$(midstr)$(p1)$(poststr).fastq.flex.gz $(*D)/$(prestr)$(*F)$(midstr)$(p2)$(poststr).fastq.flex.gz | samtools view -bS - > $(*D)/$(prestr)$(*F)$(poststr).bam  
	bwa mem -t 8 -M $(BWA_reference) $(*D)/$(prestr)$(*F)$(midstr)$(p1)$(poststr)_single.fastq.flex.gz | samtools view -bS - > $(*D)/$(prestr)$(*F)$(midstr)$(p1)$(poststr)_single.bam 
	bwa mem -t 8 -M $(BWA_reference) $(*D)/$(prestr)$(*F)$(midstr)$(p2)$(poststr)_single.fastq.flex.gz | samtools view -bS - > $(*D)/$(prestr)$(*F)$(midstr)$(p2)$(poststr)_single.bam 
else
	bwa mem -t 8 -M -R '$(rgprefix)\tID:$(prestr)$(*F)$(poststr)\tPU:$(prestr)$(*F)$(poststr)' $(BWA_reference) $(*D)/$(prestr)$(*F)$(midstr)$(p1)$(poststr).fastq.flex.gz $(*D)/$(prestr)$(*F)$(midstr)$(p2)$(poststr).fastq.flex.gz | samtools view -bS - > $(*D)/$(prestr)$(*F)$(poststr).bam  
	bwa mem -t 8 -M -R '$(rgprefix)\tID:$(prestr)$(*F)$(poststr)\tPU:$(prestr)$(*F)$(poststr)' $(BWA_reference) $(*D)/$(prestr)$(*F)$(midstr)$(p1)$(poststr)_single.fastq.flex.gz | samtools view -bS - > $(*D)/$(prestr)$(*F)$(midstr)$(p1)$(poststr)_single.bam 
	bwa mem -t 8 -M -R '$(rgprefix)\tID:$(prestr)$(*F)$(poststr)\tPU:$(prestr)$(*F)$(poststr)' $(BWA_reference) $(*D)/$(prestr)$(*F)$(midstr)$(p2)$(poststr)_single.fastq.flex.gz | samtools view -bS - > $(*D)/$(prestr)$(*F)$(midstr)$(p2)$(poststr)_single.bam 
endif

###############################################
# how to make flexbar intermediate files  
# Note that we only mention 1 of the pairs in the rule - but process both
###############################################
$(prestr)%$(midstr)$(p1)$(poststr).fastq.flex.gz:
	$(RUN_TARDIS) -w -c $(TARDIS_chunksize) -d $(TARDIS_workdir) flexbar --threads 4 --reads _condition_paired_fastq_input_$(dd)/$(prestr)$(*F)$(midstr)$(p1)$(poststr).fastq.gz --reads2 _condition_paired_fastq_input_$(dd)/$(prestr)$(*F)$(midstr)$(p2)$(poststr).fastq.gz  --format sanger --max-uncalled 2 --min-read-length 50 --pre-trim-phred 20 --adapters $(adaptersFile) --adapter-trim-end ANY --adapter-min-overlap 4 --adapter-threshold 3 --single-reads --target _condition_output_flexbartag _condition_fastq_product_{}_1.fastq,$(*D)/$(prestr)$(*F)$(midstr)$(p1)$(poststr).fastq.flex _condition_fastq_product_{}_2.fastq,$(*D)/$(prestr)$(*F)$(midstr)$(p2)$(poststr).fastq.flex _condition_fastq_product_{}_1_single.fastq,$(*D)/$(prestr)$(*F)$(midstr)$(p1)$(poststr)_single.fastq.flex _condition_fastq_product_{}_2_single.fastq,$(*D)/$(prestr)$(*F)$(midstr)$(p2)$(poststr)_single.fastq.flex \> _condition_text_output_$(*D)/$(prestr)$(*F)$(poststr).stdout


###############################################
# how to make fastqc files
###############################################
$(prestr)%$(midstr)$(p1)$(poststr)_fastqc.zip:
	$(RUN_FASTQC) $(dd)/$(prestr)$(*F)$(midstr)$(p1)$(poststr).fastq.gz -o $(*D) ; $(RUN_FASTQC) $(dd)/$(prestr)$(*F)$(midstr)$(p2)$(poststr).fastq.gz -o $(*D)

##############################################
# specify the intermediate files to keep 
##############################################
.PRECIOUS: %.coverage.sample_summary %.intervals %.samplemergedbam %.vcf %.realignedbam %.intervals $(prestr)%$(poststr).lanemergedbam $(prestr)%$(poststr).bam $(prestr)%$(poststr).sortedbam $(prestr)%$(poststr).filteredbam $(prestr)%$(poststr).bamstats $(prestr)%$(midstr)$(p1)$(poststr)_fastqc.zip $(prestr)%$(midstr)$(p1)$(poststr).fastq.flex.gz  

##############################################
# cleaning - not yet doing this using make  
##############################################
clean:
	echo "no clean for now" 

