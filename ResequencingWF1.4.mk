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
builddir=/not set

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

# ******************************************************************************************
# other variables (not project specific)
# ******************************************************************************************
RUN_TARDIS=tardis.py
RUN_FASTQC=fastqc
RUN_SAMBAMBA=/dataset/AG_1000_bulls/active/bin/sambamba
RUN_SAMTOOLS=samtools
RUN_JAVA=/usr/bin/java
GATK=/dataset/AFC_dairy_cows/archive/GenomeAnalysisTKLite-2.3-9-gdcdccbb/GenomeAnalysisTKLite.jar

# variables for tardis and other apps
TARDIS_chunksize=20000
# for testing small files
#TARDIS_chunksize=1000
TARDIS_workdir=$(builddir)
BWA_reference=not set


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
	egrep "^ln -s" $*.log >> $*.logprecis

	echo "waiting" >> $*.logprecis
	echo "-------" >> $*.logprecis
	egrep "^while >> $*.logprecis



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
	# problem - sambamba fails for very fragmented references ("More than 16383 reference sequences are unsupported")
	#$(RUN_SAMBAMBA) markdup -t 8 -r --tmpdir=$(mytmp) --overflow-list-size=800000 --io-buffer-size=256 $*.samplemergedbam_with_duplicates $*.samplemergedbam_no_duplicates
	$(RUN_SAMTOOLS) rmdup  $*.samplemergedbam_with_duplicates $*.samplemergedbam_no_duplicates
	ln -fs $*.samplemergedbam_no_duplicates $*.samplemergedbam
	ln -fs $*.samplemergedbam_no_duplicates $*_samplemerged.bam
endif
	bamtools index -in $*_samplemerged.bam 


#############################################################################
# how to make the lane-merged BAM (optionally including removal of duplicates at this level)
# - see below for notes on the way we need to expand the dependency list at this point to generate the 
# paired-file related names - i.e. the singles bams - from the sample moniker that is matched in the 
# target. Note that although the lanemerged bam depends on both the R1 and R2  singles bam, we
# only mentioned the R1 bam in the dependency list because the same process generates both the R1 and R2 
# intermediate files.
#############################################################################
#$(prestr)%$(poststr).lanemergedbam: $(prestr)%$(poststr).sortedbam  
.SECONDEXPANSION:
#%.lanemergedbam: %.sortedbam $(builddir)/$$(basename $$(subst $(poststr),$(midstr)$(p1)$(poststr),$$*)).singlessortedbam  $(builddir)/$$(basename $$(subst $(poststr),$(midstr)$(p2)$(poststr),$$*)).singlessortedbam 
%.lanemergedbam: %.sortedbam $(builddir)/$$(basename $$(subst $(poststr),$(midstr)$(p1)$(poststr),$$*)).singlessortedbam  
	#$(RUN_SAMBAMBA) merge -t 8  $(*D)/$(prestr)$(*F)$(poststr).lanemergedbam_with_duplicates $(*D)/$(prestr)$(*F)$(poststr).sortedbam $(*D)/$(prestr)$(*F)$(midstr)$(p1)$(poststr)_single.sortedbam $(*D)/$(prestr)$(*F)$(midstr)$(p2)$(poststr)_single.sortedbam
        # for some reason the directories are missing from the second and theird dependencies so can't just use $+
	$(RUN_SAMBAMBA) merge -t 8  $(basename $@).lanemergedbam_with_duplicates $< $(basename $(subst $(poststr),$(midstr)$(p1)$(poststr),$(*))).singlessortedbam $(basename $(subst $(poststr),$(midstr)$(p2)$(poststr),$(*))).singlessortedbam 
	# if removeLaneDuplicates = n , make the target (and the associated links with .bam suffix) simply by linking to the above.
	# (GATK only likes .bam or .sam suffices)
ifeq ($(subst N,n,$(strip $(removeLaneDuplicates))),n)
	ln -fs $(basename $@).lanemergedbam_with_duplicates $@ 
	ln -fs $(basename $@).lanemergedbam_with_duplicates $(basename $@)_merged.bam
else
	# or, if required to removeLaneDuplicates - process the with-duplicates merged file to remove duplicates, and link to that
	# problem - sambamba fails for very fragmented references ("More than 16383 reference sequences are unsupported")
	#$(RUN_SAMBAMBA) markdup -t 8 -r --tmpdir=$(mytmp) --overflow-list-size=800000 --io-buffer-size=256 $(*D)/$(prestr)$(*F)$(poststr).lanemergedbam_with_duplicates $(*D)/$(prestr)$(*F)$(poststr).lanemergedbam_no_duplicates
	$(RUN_SAMTOOLS) rmdup $(basename $@).lanemergedbam_with_duplicates $(basename $@).lanemergedbam_no_duplicates 
	ln -fs $(basename $@).lanemergedbam_no_duplicates $@ 
	ln -fs $(basename $@).lanemergedbam_no_duplicates $(basename $@)_merged.bam
endif
	bamtools index -in $(basename $@)_merged.bam 


###############################################
# how to make the sorted BAMs
###############################################
#$(prestr)%$(poststr).sortedbam: $(prestr)%$(poststr).filteredbam
%.sortedbam: %.filteredbam
	$(RUN_SAMBAMBA) sort --tmpdir=$(mytmp) -t 8 -m 10G $<  -o $@
%.singlessortedbam: %.singlesfilteredbam
	#$(RUN_SAMBAMBA) sort --tmpdir=$(mytmp) -t 8 -m 10G $(builddir)/$<  -o $(builddir)/$@
	$(RUN_SAMBAMBA) sort --tmpdir=$(mytmp) -t 8 -m 10G $<  -o $@

###############################################
# how to make the filtered  BAMs. We insert a dependency on bamstats 
# at this point to get those done as well 
###############################################
%.filteredbam: %.bam 
	$(RUN_SAMBAMBA) view -t8 -f bam -F "mapping_quality >=20" $< -o $@ 
%.singlesfilteredbam: %.singlesbam
	#$(RUN_SAMBAMBA) view -t8 -f bam -F "mapping_quality >=20" $(builddir)/$< -o $(builddir)/$@ 
	$(RUN_SAMBAMBA) view -t8 -f bam -F "mapping_quality >=20" $< -o $@ 


###############################################
# how to get stats for the bam files
###############################################
%.bamstats: %.bam 
	bamtools stats -in $<  > $@
%.singlesbamstats: %.singlesbam 
	#bamtools stats -in $(builddir)/$<  > $(builddir)/$@
	bamtools stats -in $<  > $@


# singlesbam depend on bam in the sense that once the flexbar prerequisite of 
# bam has been run, we can then do singlesbam , because the flexbar step also 
# delivers the pre-requisite of singlesbam
%.singlesbam:  %.singlefastq.flex.gz
ifndef $(rgvarname) 
	#bwa mem -t 8 -M $(BWA_reference) $(builddir)/$(basename $@).singlefastq.flex.gz | samtools view -bS - > $(builddir)/$@ 
	#bwa mem -t 8 -M $(BWA_reference) $(basename $@).singlefastq.flex.gz | samtools view -bS - > $@ 
	#bwa mem -t 8 -M -R '$($(rgvarname))\tID:$(*F)\tPU:$(*F)' $(BWA_reference) $(builddir)/$(basename $@).singlefastq.flex.gz | samtools view -bS - > $(builddir)/$@ 
	#bwa mem -t 8 -M -R '$($(rgvarname))\tID:$(*F)\tPU:$(*F)' $(BWA_reference) $(basename $@).singlefastq.flex.gz | samtools view -bS - > $@ 
	# we really need to run this on the farm , with the correct memory requirements set, to avoid e.g. out of RAM errors 
        # if run natively 
	$(RUN_TARDIS) -w -d $(TARDIS_workdir) bwa mem -t 4 -M $(BWA_reference) $(basename $@).singlefastq.flex.gz \| samtools view -bS - \> _condition_wait_output_$@ 
else
	#bwa mem -t 8 -M -R '$($(rgvarname))\tID:$(*F)\tPU:$(*F)' $(BWA_reference) $(basename $@).singlefastq.flex.gz | samtools view -bS - > $@ 
	$(RUN_TARDIS) -w -d $(TARDIS_workdir) bwa mem -t 4 -M -R \'$($(rgvarname))\\tID:$(*F)\\tPU:$(*F)\' $(BWA_reference) $(basename $@).singlefastq.flex.gz \| samtools view -bS - \> _condition_wait_output_$@ 
endif


###############################################
# how to make bam files. Made from 
# flexbar output from both pairs, and indirectly also 
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
%.bam: $(builddir)/$$(basename $$(subst $(poststr),$(midstr)$(p1)$(poststr),$$*))_fastqc.zip $(builddir)/$$(basename $$(subst $(poststr),$(midstr)$(p1)$(poststr),$$*)).fastq.flex.gz
ifndef $(rgvarname) 
	#bwa mem -t 8 -M $(BWA_reference) $(basename $(subst $(poststr),$(midstr)$(p1)$(poststr),$*)).fastq.flex.gz $(basename $(subst $(poststr),$(midstr)$(p2)$(poststr),$*)).fastq.flex.gz | samtools view -bS - > $@  
	#$(RUN_TARDIS) -w -d $(TARDIS_workdir) bwa mem -t 4 -M $(BWA_reference) $(basename $(subst $(poststr),$(midstr)$(p1)$(poststr),$*)).fastq.flex.gz $(basename $(subst $(poststr),$(midstr)$(p2)$(poststr),$*)).fastq.flex.gz \| samtools view -bS - \> _condition_wait_output_$@  
	$(RUN_TARDIS) -w -d $(TARDIS_workdir) bwa mem -t 4 -M $(BWA_reference) $(dir $*)/$(basename $(subst $(poststr),$(midstr)$(p1)$(poststr),$(notdir $*))).fastq.flex.gz $(dir $*)/$(basename $(subst $(poststr),$(midstr)$(p2)$(poststr),$(notdir $*))).fastq.flex.gz \| samtools view -bS - \> _condition_wait_output_$@  
else
	#bwa mem -t 8 -M -R '$($(rgvarname))\tID:$(*F)\tPU:$(*F)'  $(BWA_reference) $(basename $(subst $(poststr),$(midstr)$(p1)$(poststr),$*)).fastq.flex.gz $(basename $(subst $(poststr),$(midstr)$(p2)$(poststr),$*)).fastq.flex.gz | samtools view -bS - > $@  
	#$(RUN_TARDIS) -w -d $(TARDIS_workdir) bwa mem -t 4 -M -R \'$($(rgvarname))\\tID:$(*F)\\tPU:$(*F)\'  $(BWA_reference) $(basename $(subst $(poststr),$(midstr)$(p1)$(poststr),$*)).fastq.flex.gz $(basename $(subst $(poststr),$(midstr)$(p2)$(poststr),$*)).fastq.flex.gz \| samtools view -bS - \> _condition_wait_output_$@  
	$(RUN_TARDIS) -w -d $(TARDIS_workdir) bwa mem -t 4 -M -R \'$($(rgvarname))\\tID:$(*F)\\tPU:$(*F)\'  $(BWA_reference) $(dir $*)/$(basename $(subst $(poststr),$(midstr)$(p1)$(poststr),$(notdir $*))).fastq.flex.gz $(dir $*)/$(basename $(subst $(poststr),$(midstr)$(p2)$(poststr),$(notdir $*))).fastq.flex.gz \| samtools view -bS - \> _condition_wait_output_$@  
endif


###############################################
# how to make flexbar intermediate files  
# Note that we only mention 1 of the pairs in the rule - but process both
###############################################
%.fastq.flex.gz:
	#$(RUN_TARDIS) -w -c $(TARDIS_chunksize) -d $(TARDIS_workdir) flexbar --threads 4 --reads _condition_paired_fastq_input_$(dd)/$(subst fastq.flex.gz,fastq.gz,$(notdir $@)) --reads2 _condition_paired_fastq_input_$(dd)/$(subst fastq.flex.gz,fastq.gz,$(subst $(midstr)$(p1),$(midstr)$(p2),$(notdir $@))) --format sanger --max-uncalled 2 --min-read-length 50 --pre-trim-phred 20 --adapters $(adaptersFile) --adapter-trim-end ANY --adapter-min-overlap 4 --adapter-threshold 3 --single-reads --target _condition_output_flexbartag _condition_fastq_product_{}_1.fastq,$(builddir)/$(basename $@) _condition_fastq_product_{}_2.fastq,$(builddir)/$(subst $(midstr)$(p1),$(midstr)$(p2),$(basename $@)) _condition_fastq_product_{}_1_single.fastq,$(builddir)/$(subst .fastq.flex,.singlefastq.flex,$(basename $@)) _condition_fastq_product_{}_2_single.fastq,$(builddir)/$(subst $(midstr)$(p1),$(midstr)$(p2),$(subst .fastq.flex,.singlefastq.flex,$(basename $@))) \> _condition_text_output_$(builddir)/$(subst .fastq.flex.gz,.stdout,$(subst $(midstr)$(p1),,$@))
	$(RUN_TARDIS) -w -c $(TARDIS_chunksize) -d $(TARDIS_workdir) -batonfile $*.baton flexbar --threads 2 --reads _condition_paired_fastq_input_$(dd)/$(subst fastq.flex.gz,fastq.gz,$(notdir $@)) --reads2 _condition_paired_fastq_input_$(dd)/$(subst fastq.flex.gz,fastq.gz,$(subst $(midstr)$(p1),$(midstr)$(p2),$(notdir $@))) --format sanger --max-uncalled 2 --min-read-length 50 --pre-trim-phred 20 --adapters $(adaptersFile) --adapter-trim-end ANY --adapter-min-overlap 4 --adapter-threshold 3 --single-reads --target _condition_output_flexbartag _condition_fastq_product_{}_1.fastq,$(basename $@) _condition_fastq_product_{}_2.fastq,$(subst $(midstr)$(p1),$(midstr)$(p2),$(basename $@)) _condition_fastq_product_{}_1_single.fastq,$(subst .fastq.flex,.singlefastq.flex,$(basename $@)) _condition_fastq_product_{}_2_single.fastq,$(subst $(midstr)$(p1),$(midstr)$(p2),$(subst .fastq.flex,.singlefastq.flex,$(basename $@))) \> _condition_text_output_$(subst .fastq.flex.gz,.stdout,$(subst $(midstr)$(p1),,$@))

###############################################
# how to make flexbar singles intermediate files
# - we wait for these to be made by the flexbar process
# The flexbar process will write out a "batonfile" when its done
# (i.e. as in pass the baton)
###############################################
%.singlefastq.flex.gz:
	while [ ! -f $*.baton ]; do  sleep 1; done 

###############################################
# how to make fastqc files
###############################################
%_fastqc.zip:
	#$(RUN_FASTQC) $(dd)/$(prestr)$(*F)$(midstr)$(p1)$(poststr).fastq.gz -o $(*D) ; $(RUN_FASTQC) $(dd)/$(prestr)$(*F)$(midstr)$(p2)$(poststr).fastq.gz -o $(*D)
	$(RUN_FASTQC) $(dd)/$(subst _fastqc.zip,.fastq.gz,$(notdir $@)) -o $(builddir) ; $(RUN_FASTQC) $(dd)/$(subst $(midstr)$(p1),$(midstr)$(p2),$(subst _fastqc.zip,.fastq.gz,$(notdir $@))) -o $(builddir) 

##############################################
# specify the intermediate files to keep 
##############################################
.PRECIOUS: %.coverage.sample_summary %.intervals %.samplemergedbam %.vcf %.realignedbam %.intervals %.lanemergedbam %.bam %.singlesbam %.sortedbam %.singlessortedbam %.filteredbam %.singlesfilteredbam %.bamstats %.singlesbamstats %_fastqc.zip %.fastq.flex.gz %.singlefastq.flex.gz 

##############################################
# cleaning - not yet doing this using make  
##############################################
clean:
	echo "no clean for now" 

