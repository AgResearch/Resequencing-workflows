#
# Resequencing workflow version 5. 
#***************************************************************************************
# changes 
#***************************************************************************************
# 31/5/2015 changed to use quadtrim instead of flexbar
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
builddir=/not set
tardis_chunksize=1000000

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
# another example
# 932594_GTGAAA_L001_R2_001_single.fastq.gz
# 932594_GTGAAA_L001_R2_002_single.fastq.gz
#
# 932594_GTGAAA_  + L001 + _R + 
# 
# p1=
#
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
# Another example: - file names like C7N0GANXX-1804-04-7-1_L004_R2.fastq.gz, here 
# leave the poststr variable unset, so you would invoke the makefile wrapper script as  e.g.  
# ./ResequencingWF1.5.sh -n -S ${SAMPLE} -C condor -D /bifo/archive/fungal_data/nzgl/Raw -G /genomes/saro1.fasta  -e fungal_set -B /dataset/public_invermay_scratch/scratch/fungal_data/bwa -T /dataset/public_invermay_scratch/scratch/fungal_data/bwa -p C7N0GANXX -q 1 -r 2 -s _R  
# which would invoke this makefile as e.g. 
# make -f ResequencingWF1.5.mk -d --no-builtin-rules -j 8 -n BWA_reference=/genomes/saro1.fasta quadtrim_option_set=fungal_set dd=/bifo/archive/fungal_data/nzgl/Raw p1=1 p2=2 prestr=C7N0GANXX midstr=_R poststr= 'rgprefix=@RG\\tSM:sample1\\tPL:illumina\\tLB:sample1\\tCN:AgResearch' mytmp=/dataset/public_invermay_scratch/scratch/fungal_data/bwa/sample1tmp builddir=/dataset/public_invermay_scratch/scratch/fungal_data/bwa/sample1 removeSubjectDuplicates=y removeLaneDuplicates=n 'lanemergedBAMIncludeList= /dataset/public_invermay_scratch/scratch/fungal_data/bwa/sample1/C7N0GANXX-1804-01-4-1_L004.lanemergedbam [........ etc for other files in the data folder.......]' /dataset/public_invermay_scratch/scratch/fungal_data/bwa/sample1/sample1.all
#
# 
# Sometimes samples are spread across completely different runs or lanes and a "link farm" approach
# is needed to created a set of input files that can be matched. For example 
# consider that we want to pull files from ....
#
#/dataset/BLGsheep/archive/NZGL01418_1/Raw/C6KJ6ANXX-1418-1-05-01_ATCACG_L008_R1_001.fastq.gz
#/dataset/BLGsheep/archive/NZGL01418_1/Raw/C6KJ6ANXX-1418-1-05-01_ATCACG_L008_R2_001.fastq.gz
#/dataset/BLGsheep/archive/NZGL01418_5/Raw/H2L55BCXX-1418-1_ATCACG_L002_R1_001.fastq.gz
#/dataset/BLGsheep/archive/NZGL01418_5/Raw/H2L55BCXX-1418-1_ATCACG_L002_R2_001.fastq.gz
#
#we can achieve this by creating a folder containing shortcuts to the actual files, with a
#consistent prefix in the link name...
#
#mkdir /dataset/AG_1000_sheep/scratch/general_processing_062015/linkfarms
#mkdir /dataset/AG_1000_sheep/scratch/general_processing_062015/linkfarms/NZUKNM100017918713
#for file in /dataset/BLGsheep/archive/NZGL01418_1/Raw/*ATCACG*.fastq.gz; do
#   base=`basename $file`
#   ln -s $file /dataset/AG_1000_sheep/scratch/general_processing_062015/linkfarms/NZUKNM100017918713/prefix$base
#done
#for file in /dataset/BLGsheep/archive/NZGL01418_5/Raw/*ATCACG*.fastq.gz; do
#   base=`basename $file`
#   ln -s $file /dataset/AG_1000_sheep/scratch/general_processing_062015/linkfarms/NZUKNM100017918713/prefix$base
#done
# 
#
#./ResequencingWF1.5.sh -S NZUKNM100017918713 -D /dataset/AG_1000_sheep/scratch/general_processing_062015/linkfarms/NZUKNM100017918713  -G /dataset/datacache/scratch/ensembl/oar3.1/indexes/Ovis_aries.Oar_v3.1.dna_sm.toplevel.fa -B /dataset/AG_1000_sheep/scratch/general_processing_062015 -T /dataset/AG_1000_sheep/scratch/general_processing_062015 -p prefix -q 1 -r 2 -s _R -t _0
# 
# 
# ******************************************************************************************
# other variables (not project specific)
# ******************************************************************************************
RUN_TARDIS=/home/mccullocha/galaxy/hpc/dev/tardis.py
RUN_FASTQC=fastqc
RUN_SAMBAMBA=sambamba
RUN_SAMTOOLS=samtools
RUN_JAVA=java
RUN_QUADTRIM=quadtrim
GET_FASTQ=./get_fastq.sh
# quadtrim option sets  - we are passed (via the variable quadtrim_option_set) set names like "sheep_set" etc, rather than -d sheep ,
sheep_set=-d sheep
cattle_set=-d bulls

# variables for tardis and other apps
tardis_workdir=$(builddir)
BWA_reference=not set


#*****************************************************************************************
# from here are the targets and rules
#*****************************************************************************************


###############################################
# top level phony test targets  (used for various testing / debugging)
###############################################

###############################################
# top level target "logprecis" . This extracts from the log file the 
# relevant commands that were run, in a readable order 
# - run this after the actual build has been completed
###############################################
%.logprecis: %.log
	echo "get_fastq" > $*.logprecis
	echo "---------" >> $*.logprecis
	egrep "^$(GET_FASTQ)" $*.log >> $*.logprecis

	echo "fastqc" >> $*.logprecis
	echo "------" >> $*.logprecis
	egrep "^$(RUN_FASTQC)" $*.log >> $*.logprecis

	echo "quadtrim" >> $*.logprecis
	echo "-------" >> $*.logprecis
	grep "$(RUN_QUADTRIM)" $*.log  | grep -v "tool args" | egrep -v "^#" >> $*.logprecis

	echo "bwa" >> $*.logprecis
	echo "---" >> $*.logprecis
	grep "bwa mem" $*.log  | grep -v "tool args" | egrep -v "^#"  | grep -v "CMD" >> $*.logprecis

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
	egrep "^$(RUN_SAMTOOLS) rmdup" $*.log >> $*.logprecis

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

	echo "file linking" >> $*.logprecis
	echo "------------" >> $*.logprecis
	egrep "^ln -fs" $*.log >> $*.logprecis

	echo "waiting" >> $*.logprecis
	echo "-------" >> $*.logprecis
	egrep "^while" $*.log >> $*.logprecis



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
%.versions:
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
	echo $(RUN_JAVA) -Xmx4G -jar -Djava.io.tmpdir=$(mytmp) -jar $(GATK_LITE_JAR) >> versions.log
	$(RUN_JAVA) -Xmx4G -jar -Djava.io.tmpdir=$(mytmp) -jar $(GATK_LITE_JAR) -h | head -5 >> versions.log  2>&1
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
	echo "quadtrim"  >> versions.log
	echo "-------"  >> versions.log
	/dataset/hiseq/active/bin/quadtrim -h | grep Version >> versions.log


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
	$(RUN_JAVA) -Xmx4G -jar -Djava.io.tmpdir=$(mytmp) -jar $(GATK_LITE_JAR) -T DepthOfCoverage -R $(BWA_reference) -I $*_realigned.bam --omitDepthOutputAtEachBase --logging_level ERROR --summaryCoverageThreshold 10 --summaryCoverageThreshold 20 --summaryCoverageThreshold 30 --summaryCoverageThreshold 40 --summaryCoverageThreshold 50 --summaryCoverageThreshold 80 --summaryCoverageThreshold 90 --summaryCoverageThreshold 100 --summaryCoverageThreshold 150 --minBaseQuality 15 --minMappingQuality 30 --start 1 --stop 1000 --nBins 999 -dt NONE -o $*.coverage

###############################################
# how to make the vcf
###############################################
%.vcf: %.realignedbam
	#samtools mpileup -ugf $(BWA_reference) $*_realigned.bam | bcftools view -N -cvg - > $*.vcf 
	samtools mpileup -q 20 -ugf $(BWA_reference) $*_realigned.bam | bcftools view -N -cvg - > $*.vcf 


###############################################
# how to make the "full" vcf - output everything not just variant positions 
# example : make -f ResequencingWF1.3.mk -d --no-builtin-rules /dataset/AG_1000_bulls/scratch/JJLNZLM000000000402/JJLNZLM000000000402.fullvcf
# (consider adding fullvcf as a dependency of "all", or else just change the standard 
# vcf build so it includes all positions ?)
###############################################
%.fullvcf: %.realignedbam
	#samtools mpileup -ugf $(BWA_reference) $*_realigned.bam | bcftools view -N -cg - > $*.fullvcf
	samtools mpileup -q 20 -ugf $(BWA_reference) $*_realigned.bam | bcftools view -N -cg - > $*.fullvcf



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
	$(RUN_JAVA) -Xmx10g -Djava.io.tmpdir=$(mytmp) -jar $(GATK_LITE_JAR) -T IndelRealigner --read_filter NotPrimaryAlignment -R $(BWA_reference) -targetIntervals $*.intervals  -I $*_samplemerged.bam  -o $*.realignedbam -S LENIENT
else
%.realignedbam: %.samplemergedbam 
	#$(RUN_JAVA) -Xmx10g -Djava.io.tmpdir=$(mytmp) -jar $(GATK) -T IndelRealigner -R $(BWA_reference) -targetIntervals $(targetIntervals) -I $*_samplemerged.bam -o $*.realignedbam -S LENIENT
	$(RUN_JAVA) -Xmx10g -Djava.io.tmpdir=$(mytmp) -jar $(GATK_LITE_JAR) -T IndelRealigner --read_filter NotPrimaryAlignment -R $(BWA_reference) -targetIntervals $(targetIntervals) -I $*_samplemerged.bam -o $*.realignedbam -S LENIENT
endif
	# make a shortcut with suffix "bam" - many of the tools do not like custom suffixes like realignedbam etc
	ln -fs $*.realignedbam $*_realigned.bam 
	bamtools index -in $*_realigned.bam 

###############################################
# how to make the "intervals" files
###############################################
%.intervals: %.samplemergedbam
	# if we have samplemergedbam then we should also have _samplemerged.bam plus index
	$(RUN_JAVA) -Xmx4g -Djava.io.tmpdir=$(mytmp) -jar $(GATK_LITE_JAR) -T RealignerTargetCreator -nt 16 -R $(BWA_reference) -I $*_samplemerged.bam -o $*.intervals -S LENIENT -rbs 5000000


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
	# problem - sambamba fails for very fragmented references ("More than 16383 reference sequences are unsupported")
	#$(RUN_SAMBAMBA) markdup -t 8 -r --tmpdir=$(mytmp) --overflow-list-size=800000 --io-buffer-size=256 $*.samplemergedbam_with_duplicates $*.samplemergedbam_no_duplicates
	$(RUN_SAMTOOLS) rmdup  $*.samplemergedbam_with_duplicates $*.samplemergedbam_no_duplicates
	ln -fs $*.samplemergedbam_no_duplicates $*.samplemergedbam
	ln -fs $*.samplemergedbam_no_duplicates $*_samplemerged.bam
endif
	bamtools index -in $*_samplemerged.bam 


#############################################################################
# how to make the lane-merged BAM (optionally including removal of duplicates at this level)
#############################################################################
ifeq ($(strip $(poststr)),)
.SECONDEXPANSION:
%.lanemergedbam: %.sortedbam
	#%.lanemergedbam: %.sortedbam %_RX.singlessortedbam  
	#$(RUN_SAMBAMBA) merge -t 8  $(basename $@).lanemergedbam_with_duplicates $? 
	$(RUN_SAMBAMBA) merge -t 8  $(basename $@).lanemergedbam_with_duplicates $< $*.singlessortedbam 
else
.SECONDEXPANSION:
%.lanemergedbam: %.sortedbam    
	#%.lanemergedbam: %.sortedbam $(builddir)/$$(basename $$(subst $(poststr),$(midstr)X$(poststr),$$*)).singlessortedbam  
	#$(RUN_SAMBAMBA) merge -t 8  $(basename $@).lanemergedbam_with_duplicates $< $(builddir)/$(basename $(subst $(poststr),$(midstr)X$(poststr),$(notdir $*))).singlessortedbam
	$(RUN_SAMBAMBA) merge -t 8  $(basename $@).lanemergedbam_with_duplicates $< $*.singlessortedbam
endif
	# if removeLaneDuplicates = n , make the target (and the associated links with .bam suffix) simply by linking to the above.
	# (GATK only likes .bam or .sam suffices)
ifeq ($(subst N,n,$(strip $(removeLaneDuplicates))),n)
	ln -fs $(basename $@).lanemergedbam_with_duplicates $@ 
	ln -fs $(basename $@).lanemergedbam_with_duplicates $(basename $@)_merged.bam
else
	# or, if required to removeLaneDuplicates - process the with-duplicates merged file to remove duplicates, and link to that
	# problem - sambamba fails for very fragmented references ("More than 16383 reference sequences are unsupported")
	$(RUN_SAMTOOLS) rmdup $(basename $@).lanemergedbam_with_duplicates $(basename $@).lanemergedbam_no_duplicates 
	ln -fs $(basename $@).lanemergedbam_no_duplicates $@ 
	ln -fs $(basename $@).lanemergedbam_no_duplicates $(basename $@)_merged.bam
endif
	bamtools index -in $(basename $@)_merged.bam 


###############################################
# how to make the sorted BAMs
###############################################
%.sortedbam: %.filteredbam
	$(RUN_SAMBAMBA) sort --tmpdir=$(mytmp) -t 8 -m 10G $<  -o $@
	$(RUN_SAMBAMBA) sort --tmpdir=$(mytmp) -t 8 -m 10G $*.singlesfilteredbam -o $*.singlessortedbam


# this rule deprecated now handled as side-effect above
#%.singlessortedbam: %.singlesfilteredbam
#	$(RUN_SAMBAMBA) sort --tmpdir=$(mytmp) -t 8 -m 10G $<  -o $@

###############################################
# how to make the filtered  BAMs. 
###############################################
%.filteredbam: %.bam 
	# do e.g. this if we want to filter - but no longer filtering here so just pass the unfiltered through 
	#$(RUN_SAMBAMBA) view -t8 -f bam -F "mapping_quality >=20" $< -o $@ 
	ln -s $< $@ 
	ln -s $*.singlesbam $*.singlesfilteredbam
	

# this rule deprecated now handled as side effect above
#%.singlesfilteredbam: %.singlesbam
#	#$(RUN_SAMBAMBA) view -t8 -f bam -F "mapping_quality >=20" $< -o $@ 
#	ln -s $< $@ 


###############################################
# how to get stats for the bam files
###############################################
%.bamstats: %.bam 
	bamtools stats -in $<  > $@
	bamtools stats -in $*.singlesbam  > $*.singlesbamstats
	
# this rule deprecated now handled as side-effect above
#%.singlesbamstats: %.singlesbam 
#	#bamtools stats -in $(builddir)/$<  > $(builddir)/$@
#	bamtools stats -in $<  > $@


# Deprecated , singles now handled as side effect  
#singlesbam depend on bam in the sense that once the quadtrim prerequisite of 
# bam has been run, we can then do singlesbam , because the quadtrim step also 
# delivers the pre-requisite of singlesbam
# Note we use tardis to ensure taht we wait for the output  (** need to check this)
#%.singlesbam:  %.singlefastq.quadtrim
#ifndef $(rgvarname) 
#	$(RUN_TARDIS) -w -d $(tardis_workdir) bwa mem -t 4 -M $(BWA_reference) $(basename $@).singlefastq.quadtrim \| samtools view -bS - \> _condition_wait_output_$@ 
#else
#	$(RUN_TARDIS) -w -d $(tardis_workdir) bwa mem -t 4 -M -R \'$($(rgvarname))\\tID:$(*F)\\tPU:$(*F)\' $(BWA_reference) $(basename $@).singlefastq.quadtrim \| samtools view -bS - \> _condition_wait_output_$@ 
#endif


###############################################
# how to make bam files. Made from 
# quadtrim output from both pairs, and indirectly also 
# FastQC output 
# Note that we only mention 1 of the pairs in the rule - but process both
# We also optionally inject readgroup info at this point
###############################################
#  The point of the dependency construction in the next rule is that from 
#  a sample moniker like 
# 932594_GTGAAA_L001_001
# we need to construct the paired end filenames - e.g. 
#932594_GTGAAA_L001_R1_001
#932594_GTGAAA_L001_R2_001
# - hence the approach of passing in filename "landmarks" p1, p2, prestr, midstr, poststr
# e.g. in the above example we pass in 
# p1=1 p2=2 prestr=9 midstr=_R poststr=_001 
# KNOWN BUG: 
# for some reason the dependency expansion in the next line strips off the directory part of the 
# target - e.g. if we are constructing the paired-end file dependencies from /builddir/mine.bam
# then the dependencies are constucted as munged_mine_fastqc.zip etc - rather than /builddir/munged_mine_fastqc.zip
# Currrent work-around is to add back in the build dir in the rule for these dependencies - i.e. 
# $(builddir)
.SECONDEXPANSION:
%.bam:  $(builddir)/$$(basename $$(subst $(poststr),$(midstr)$(p1)$(poststr),$$*)).fastq.quadtrim
ifndef $(rgvarname) 
	#$(RUN_TARDIS) -w -d $(tardis_workdir) bwa mem -t 4 -M $(BWA_reference) $(dir $*)/$(basename $(subst $(poststr),$(midstr)$(p1)$(poststr),$(notdir $*))).fastq.quadtrim $(dir $*)/$(basename $(subst $(poststr),$(midstr)$(p2)$(poststr),$(notdir $*))).fastq.quadtrim \| samtools view -bS - \> _condition_wait_output_$@  
	$(RUN_TARDIS) -w -d $(tardis_workdir) -c $(tardis_chunksize) bwa mem -t 4 -M $(BWA_reference) _condition_paired_fastq_input_$(dir $*)/$(basename $(subst $(poststr),$(midstr)$(p1)$(poststr),$(notdir $*))).fastq.quadtrim $(dir $*)/$(basename _condition_paired_fastq_input_$(subst $(poststr),$(midstr)$(p2)$(poststr),$(notdir $*))).fastq.quadtrim \> _condition_sam_output_$@  
else
	#$(RUN_TARDIS) -w -d $(tardis_workdir) bwa mem -t 4 -M -R \'$($(rgvarname))\\tID:$(*F)\\tPU:$(*F)\'  $(BWA_reference) $(dir $*)/$(basename $(subst $(poststr),$(midstr)$(p1)$(poststr),$(notdir $*))).fastq.quadtrim $(dir $*)/$(basename $(subst $(poststr),$(midstr)$(p2)$(poststr),$(notdir $*))).fastq.quadtrim \| samtools view -bS - \> _condition_wait_output_$@  
	$(RUN_TARDIS) -w -d $(tardis_workdir) -c $(tardis_chunksize) bwa mem -t 4 -M -R \'$($(rgvarname))\\tID:$(*F)\\tPU:$(*F)\'  $(BWA_reference) _condition_paired_fastq_input_$(dir $*)/$(basename $(subst $(poststr),$(midstr)$(p1)$(poststr),$(notdir $*))).fastq.quadtrim _condition_paired_fastq_input_$(dir $*)/$(basename $(subst $(poststr),$(midstr)$(p2)$(poststr),$(notdir $*))).fastq.quadtrim \> _condition_sam_output_$@  
endif
# also handle singletons in this rule 
ifndef $(rgvarname) 
	#$(RUN_TARDIS) -w -d $(tardis_workdir) bwa mem -t 4 -M $(BWA_reference) $*$(midstr)X.singlefastq.quadtrim \| samtools view -bS - \> _condition_wait_output_$*.singlesbam 
	$(RUN_TARDIS) -w -d $(tardis_workdir) -c $(tardis_chunksize) bwa mem -t 4 -M $(BWA_reference) _condition_fastq_input_$*$(midstr)X.singlefastq.quadtrim \> _condition_sam_output_$*.singlesbam 
else
	#$(RUN_TARDIS) -w -d $(tardis_workdir) bwa mem -t 4 -M -R \'$($(rgvarname))\\tID:$(*F)\\tPU:$(*F)\' $(BWA_reference) $*$(midstr)X.singlefastq.quadtrim \| samtools view -bS - \> _condition_wait_output_$*.singlesbam 
	$(RUN_TARDIS) -w -d $(tardis_workdir) -c $(tardis_chunksize) bwa mem -t 4 -M -R \'$($(rgvarname))\\tID:$(*F)\\tPU:$(*F)\' $(BWA_reference) _condition_fastq_input_$*$(midstr)X.singlefastq.quadtrim \> _condition_sam_output_$*.singlesbam 
endif





###############################################
# Deprecated how to make quadtrim singles intermediate files
# (now handled as side-effect for simplicity)
###############################################
#%.singlefastq.quadtrim: %.fastq.quadtrim
#	echo $@ done as side effect of  $<


###############################################
# how to make quadtrim intermediate files  
# Note that we only mention 1 of the pairs in the rule - but process both
# Other output that is made as a side-effect of this rule : 
# %_fastqc.zip 			(for simplicity ) 
# %.singlefastq.quadtrim  	(side effect of quadtrim)
###############################################
%.fastq.quadtrim: %.fastq.gz
ifeq ($(strip $(poststr)),)
	#$(RUN_TARDIS) -w -c $(tardis_chunksize) -d $(tardis_workdir) -batonfile $*.baton $(RUN_QUADTRIM) $($(quadtrim_option_set)) _condition_paired_fastq_input_$(dd)/$(subst fastq.quadtrim,fastq.gz,$(notdir $@)) _condition_paired_fastq_input_$(dd)/$(subst fastq.quadtrim,fastq.gz,$(subst $(midstr)$(p1),$(midstr)$(p2),$(notdir $@))) '_condition_uncompressedfastq_product_\S*?_R1\.\d{5}-pass.fq,$@'  '_condition_uncompressedfastq_product_\S*?_R2\.\d{5}-pass.fq,$(*D)/$(subst $(midstr)$(p1),$(midstr)$(p2),$(notdir $@))' '_condition_uncompressedfastq_product_\S*?_RX\.\d{5}-singleton.fq,$(*D)/$(subst $(midstr)$(p1),$(midstr)X,$(subst fastq,singlefastq,$(notdir $@)))' \> _condition_text_output_$(subst .fastq.quadtrim,.stdout,$(subst $(midstr)$(p1),,$@))
	#$(RUN_TARDIS) -w -c $(tardis_chunksize) -d $(tardis_workdir) -batonfile $*.baton $(RUN_QUADTRIM) $($(quadtrim_option_set)) _condition_paired_fastq_input_$(subst fastq.quadtrim,fastq.gz,$@) _condition_paired_fastq_input_$(subst fastq.quadtrim,fastq.gz,$(subst $(midstr)$(p1),$(midstr)$(p2),$@)) '_condition_uncompressedfastq_product_\S*?_R1\.\d{5}-pass.fq,$@'  '_condition_uncompressedfastq_product_\S*?_R2\.\d{5}-pass.fq,$(*D)/$(subst $(midstr)$(p1),$(midstr)$(p2),$(notdir $@))' '_condition_uncompressedfastq_product_\S*?_RX\.\d{5}-singleton.fq,$(*D)/$(subst $(midstr)$(p1),$(midstr)X,$(subst fastq,singlefastq,$(notdir $@)))' \> _condition_text_output_$(subst .fastq.quadtrim,.stdout,$(subst $(midstr)$(p1),,$@))
	#$(RUN_TARDIS) -w -c $(tardis_chunksize) -d $(tardis_workdir) -batonfile $*.baton $(RUN_QUADTRIM) $($(quadtrim_option_set)) _condition_paired_fastq_input_$(subst fastq.quadtrim,fastq.gz,$@) _condition_paired_fastq_input_$(subst fastq.quadtrim,fastq.gz,$(subst $(midstr)$(p1),$(midstr)$(p2),$@)) '_condition_uncompressedfastq_product_\S*?$(midstr)1\.\d{5}-pass.fq,$@'  '_condition_uncompressedfastq_product_\S*?$(midstr)2\.\d{5}-pass.fq,$(*D)/$(subst $(midstr)$(p1),$(midstr)$(p2),$(notdir $@))' '_condition_uncompressedfastq_product_\S*?$(midstr)X\.\d{5}-singleton.fq,$(*D)/$(subst $(midstr)$(p1),$(midstr)X,$(subst fastq,singlefastq,$(notdir $@)))' \> _condition_text_output_$(subst .fastq.quadtrim,.stdout,$(subst $(midstr)$(p1),,$@))
	$(RUN_TARDIS) -w -c $(tardis_chunksize) -d $(tardis_workdir) -batonfile $*.baton $(RUN_QUADTRIM) $($(quadtrim_option_set)) _condition_paired_fastq_input_$(subst fastq.quadtrim,fastq.gz,$@) _condition_paired_fastq_input_$(*D)/$(subst fastq.quadtrim,fastq.gz,$(subst $(midstr)$(p1),$(midstr)$(p2),$(notdir $@))) '_condition_uncompressedfastq_product_\S*?$(midstr)1\.\d{5}-pass.fq,$@'  '_condition_uncompressedfastq_product_\S*?$(midstr)2\.\d{5}-pass.fq,$(*D)/$(subst $(midstr)$(p1),$(midstr)$(p2),$(notdir $@))' '_condition_uncompressedfastq_product_\S*?$(midstr)X\.\d{5}-singleton.fq,$(*D)/$(subst $(midstr)$(p1),$(midstr)X,$(subst fastq,singlefastq,$(notdir $@)))' \> _condition_text_output_$(*D)/$(subst .fastq.quadtrim,.stdout,$(subst $(midstr)$(p1),,$(notdir $@)))
else
	#$(RUN_TARDIS) -w -c $(tardis_chunksize) -d $(tardis_workdir) -batonfile $*.baton $(RUN_QUADTRIM) $($(quadtrim_option_set)) _condition_paired_fastq_input_$(dd)/$(subst fastq.quadtrim,fastq.gz,$(notdir $@)) _condition_paired_fastq_input_$(dd)/$(subst fastq.quadtrim,fastq.gz,$(subst $(midstr)$(p1),$(midstr)$(p2),$(notdir $@))) '_condition_uncompressedfastq_product_\S*?_R1_\S*?\d{5}-pass.fq,$@'  '_condition_uncompressedfastq_product_\S*?_R2_\S*?\d{5}-pass.fq,$(*D)/$(subst $(midstr)$(p1),$(midstr)$(p2),$(notdir $@))' '_condition_uncompressedfastq_product_\S*?_RX_\S*?\d{5}-singleton.fq,$(*D)/$(subst $(midstr)$(p1),$(midstr)X,$(subst fastq,singlefastq,$(notdir $@)))' \> _condition_text_output_$(subst .fastq.quadtrim,.stdout,$(subst $(midstr)$(p1),,$@))
	#$(RUN_TARDIS) -w -c $(tardis_chunksize) -d $(tardis_workdir) -batonfile $*.baton $(RUN_QUADTRIM) $($(quadtrim_option_set)) _condition_paired_fastq_input_$(subst fastq.quadtrim,fastq.gz,$@) _condition_paired_fastq_input_$(subst fastq.quadtrim,fastq.gz,$(subst $(midstr)$(p1),$(midstr)$(p2),$@)) '_condition_uncompressedfastq_product_\S*?_R1_\S*?\d{5}-pass.fq,$@'  '_condition_uncompressedfastq_product_\S*?_R2_\S*?\d{5}-pass.fq,$(*D)/$(subst $(midstr)$(p1),$(midstr)$(p2),$(notdir $@))' '_condition_uncompressedfastq_product_\S*?_RX_\S*?\d{5}-singleton.fq,$(*D)/$(subst $(midstr)$(p1),$(midstr)X,$(subst fastq,singlefastq,$(notdir $@)))' \> _condition_text_output_$(subst .fastq.quadtrim,.stdout,$(subst $(midstr)$(p1),,$@))
	$(RUN_TARDIS) -w -c $(tardis_chunksize) -d $(tardis_workdir) -batonfile $*.baton $(RUN_QUADTRIM) $($(quadtrim_option_set)) _condition_paired_fastq_input_$(subst fastq.quadtrim,fastq.gz,$@) _condition_paired_fastq_input_$(*D)/$(subst fastq.quadtrim,fastq.gz,$(subst $(midstr)$(p1),$(midstr)$(p2),$(notdir $@))) '_condition_uncompressedfastq_product_\S*?$(midstr)1_\S*?\d{5}-pass.fq,$@'  '_condition_uncompressedfastq_product_\S*?$(midstr)2_\S*?\d{5}-pass.fq,$(*D)/$(subst $(midstr)$(p1),$(midstr)$(p2),$(notdir $@))' '_condition_uncompressedfastq_product_\S*?$(midstr)X_\S*?\d{5}-singleton.fq,$(*D)/$(subst $(midstr)$(p1),$(midstr)X,$(subst fastq,singlefastq,$(notdir $@)))' \> _condition_text_output_$(*D)/$(subst .fastq.quadtrim,.stdout,$(subst $(midstr)$(p1),,$(notdir $@)))
endif
	$(RUN_FASTQC) $< -o $(builddir) ; $(RUN_FASTQC) $(*D)/$(subst $(midstr)$(p1),$(midstr)$(p2),$(notdir $<)) -o $(builddir) 		


###############################################
# Deprectaed seperate rule - how to make fastqc files
###############################################
#%_fastqc.zip: %.fastq.gz
#	#$(RUN_FASTQC) $(dd)/$(subst _fastqc.zip,.fastq.gz,$(notdir $@)) -o $(builddir) ; $(RUN_FASTQC) $(dd)/$(subst $(midstr)$(p1),$(midstr)$(p2),$(subst _fastqc.zip,.fastq.gz,$(notdir $@))) -o $(builddir) 
#	$(RUN_FASTQC) $(subst _fastqc.zip,.fastq.gz,$@) -o $(builddir) ; $(RUN_FASTQC) $(subst $(midstr)$(p1),$(midstr)$(p2),$(subst _fastqc.zip,.fastq.gz,$@)) -o $(builddir) 

###############################################
# how to make fastq files if we don't have them
###############################################
%.fastq.gz:
ifeq ($(input_method),sra)
	# (this will extract fastq from sra archive)
	$(GET_FASTQ) $(input_method) $(dd)/$(notdir $@) $@
else 
	# (this will just make links to the existing files)
	$(GET_FASTQ) $(input_method) $(dd)/$(notdir $@) $@ 
	$(GET_FASTQ) $(input_method) $(dd)/$(subst $(midstr)$(p1),$(midstr)$(p2),$(notdir $@)) $(subst $(midstr)$(p1),$(midstr)$(p2),$@) 
endif


##############################################
# specify the intermediate files to keep 
##############################################
.PRECIOUS: %.coverage.sample_summary %.intervals %.samplemergedbam %.vcf %.realignedbam %.intervals %.lanemergedbam %.bam %.singlesbam %.sortedbam %.singlessortedbam %.filteredbam %.singlesfilteredbam %.bamstats %.singlesbamstats %_fastqc.zip %.fastq.quadtrim %.singlefastq.quadtrim %.fastq.gz

##############################################
# cleaning - not yet doing this using make  
##############################################
clean:
	echo "no clean for now" 

