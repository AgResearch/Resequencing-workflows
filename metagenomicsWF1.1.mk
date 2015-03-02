#
# Metagenomics workflow version 1 - taxonomy assignment using blast via qiime scripts. 
# This workflow is experimental !  
#***************************************************************************************
# changes 
#***************************************************************************************
#
#***************************************************************************************
# known bugs 
#***************************************************************************************
# > loading the Qiime module up front causes the tardis based step to fail  
#   because that environment does not include biopython (needed by tardis)
#   (workaround - call make twice - the first time without the qiime 
#   module, then again for the summary phase after loading the qiime module)
# > the qiime summary script crashes on large datasets. The workaround is to summarise
#   each sample seperately, and then paste the individual sample summaries 
#   together "manually" using a custom script. Example of the manual merge : 
#   ./tax_summary_merge.py `find /dataset/Kittelmann_Buccal_Ill/scratch/nzgl01005/build1/processed* -name "*.combined_L4.txt" -print` 
#   - this final step could be added to this workflow but has not been done yet.
# 
#***************************************************************************************
# references:
#***************************************************************************************
# make: 
#     http://www.gnu.org/software/make/manual/make.html
# qiime:
#     http://qiime.org/scripts/assign_taxonomy.html
#     http://www.wernerlab.org/teaching/qiime/overview/c
#
# Note - you should avoid making your targets in the same folder as your source data (in case you overwrite or clean up
# something important). Create an empty build folder, and then your target to make is empty_build_folder/animal_name.all    
#
# Note - the folder where you run the make should contain a tardis config file with a section like this: 
#
#[tardis_engine]
#shell_template_name=qiime_shell
#job_template_name=condor_send_env_job
#seqlength_min=200   # if you want to restrict the run to include sequences >= 200
#
# - the shell template ensures that the Qiime scripts will work.
# The job template ensures that the module environment variables
# (needed to load the qiime module) are available in the shell that 
# runs the script
#
#*****************************************************************************************
# variables 
#*****************************************************************************************
# The following variables can be passed in on the command-line
#
#*****************************************************************************************
# examples  
#*****************************************************************************************

# ******************************************************************************************
# initialise lists of targets , derived from the input files
# ******************************************************************************************
mytmp=/tmp

R1join_names :=  $(addsuffix .join, $(basename $(notdir $(wildcard $(run1dir)/*_R1_*.trimmed))))
# example : processed_Undetermined_S0_L001_R1_001.fastq.trimmed.join
joined_pair_monikers := $(addprefix $(builddir)/, $(foreach name, $(R1join_names), $(subst _R1_,_,$(name))))
# example : /dataset/Kittelmann_Buccal_Ill/scratch/nzgl01005/build1/processed_S398Rumen-75MG_S87_L001_001.fastq.join
joined_pair_tax_assignment_monikers := $(addsuffix .tax_assignments, $(joined_pair_monikers))
# example : /dataset/Kittelmann_Buccal_Ill/scratch/nzgl01005/build1/processed_S398Rumen-75MG_S87_L001_001.fastq.join.tax_assignments


R1tax_assignments :=  $(addprefix $(builddir)/, $(addsuffix .combined.tax_assignments, $(notdir $(wildcard $(run1dir)/*_R1_*.trimmed))))
# example : processed_Undetermined_S0_L001_R1_001.fastq.trimmed.combined.tax_assignments
R2tax_assignments :=  $(addprefix $(builddir)/, $(addsuffix .combined.tax_assignments, $(notdir $(wildcard $(run1dir)/*_R2_*.trimmed))))
# example : processed_Undetermined_S0_L001_R2_001.fastq.trimmed.combined.tax_assignments


454tax_assignments :=  $(addprefix $(builddir)/, $(addsuffix .454tax_assignments, $(notdir $(wildcard $(run454dir)/*repset.txt))))
# example : /dataset/Kittelmann_Buccal_Ill/scratch/nzgl01005/build454/All_repset.txt.tax_assignments 

454otu_map :=  $(addprefix $(builddir)/, $(addsuffix .otu_map, $(notdir $(wildcard $(run454dir)/*repset.txt))))
# example : /dataset/Kittelmann_Buccal_Ill/scratch/nzgl01005/build454/All_repset.txt.otu_map 

R1tax_assignment_summaries :=  $(addprefix $(builddir)/, $(addsuffix .combined.tax_assignments.summary, $(notdir $(wildcard $(run1dir)/*_R1_*.trimmed))))
# example : processed_Undetermined_S0_L001_R1_001.fastq.trimmed.combined.tax_assignments.summary
R2tax_assignment_summaries :=  $(addprefix $(builddir)/, $(addsuffix .combined.tax_assignments.summary, $(notdir $(wildcard $(run1dir)/*_R2_*.trimmed))))
# example : processed_Undetermined_S0_L001_R2_001.fastq.trimmed.combined.tax_assignments.summary
joined_pair_tax_assignment_summaries :=  $(addsuffix .tax_assignments.summary, $(joined_pair_monikers))
# example : processed_S398Rumen-75MG_S87_L001_001.fastq.join.tax_assignments.summary

R1tax_assignments_fake_otu_map :=  $(addprefix $(builddir)/, $(addsuffix .combined.tax_assignments.fake_otu_map, $(notdir $(wildcard $(run1dir)/*_R1_*.trimmed))))
# example : processed_Undetermined_S0_L001_R1_001.fastq.trimmed.combined.tax_assignments.fake_otu_map
R2tax_assignments_fake_otu_map :=  $(addprefix $(builddir)/, $(addsuffix .combined.tax_assignments.fake_otu_map, $(notdir $(wildcard $(run1dir)/*_R2_*.trimmed))))
# example : processed_Undetermined_S0_L001_R2_001.fastq.trimmed.combined.tax_assignments.fake_otu_map
joined_pair_tax_assignments_fake_otu_map :=  $(addsuffix .tax_assignments.fake_otu_map, $(joined_pair_monikers))
# example : processed_S398Rumen-75MG_S87_L001_001.fastq.join.tax_assignments.fake_otu_map
joined_pair_tax_assignments_patched_with_otuid :=  $(addsuffix .tax_assignments.patched_with_otuid, $(joined_pair_monikers))
# example : processed_S398Rumen-75MG_S87_L001_001.fastq.join.tax_assignments.patched_with_otuid

R1tax_assignments_patched_with_otuid := $(addprefix $(builddir)/, $(addsuffix .combined.tax_assignments.patched_with_otuid, $(notdir $(wildcard $(run1dir)/*_R1_*.trimmed))))
# example : processed_Undetermined_S0_L001_R1_001.fastq.trimmed.combined.tax_assignments.patched_with_otuid
R2tax_assignments_patched_with_otuid := $(addprefix $(builddir)/, $(addsuffix .combined.tax_assignments.patched_with_otuid, $(notdir $(wildcard $(run1dir)/*_R2_*.trimmed))))
# example : processed_Undetermined_S0_L001_R2_001.fastq.trimmed.combined.tax_assignments.patched_with_otuid

R1tax_assignments_sampled :=  $(addprefix $(builddir)/, $(addsuffix .combined.tax_assignments_sampled, $(notdir $(wildcard $(run1dir)/*_R1_*.trimmed))))
# example : processed_Undetermined_S0_L001_R1_001.fastq.trimmed.combined.tax_assignments_sampled

R1tax_assignments_sampled_fake_otu_map :=  $(addprefix $(builddir)/, $(addsuffix .combined.tax_assignments_sampled.fake_otu_map, $(notdir $(wildcard $(run1dir)/*_R1_*.trimmed))))
# example : processed_Undetermined_S0_L001_R1_001.fastq.trimmed.combined.tax_assignments_sampled.fake_otu_map

R1tax_assignments_sampled_patched_with_otuid := $(addprefix $(builddir)/, $(addsuffix .combined.tax_assignments_sampled.patched_with_otuid, $(notdir $(wildcard $(run1dir)/*_R1_*.trimmed))))
# example : processed_Undetermined_S0_L001_R1_001.fastq.trimmed.combined.tax_assignments_sampled.patched_with_otuid

# ******************************************************************************************
# other variables (not project specific)
# ******************************************************************************************
#RUN_TARDIS=tardis.py
RUN_TARDIS=/home/mccullocha/galaxy/hpc/dev/tardis.py 
RUN_FASTQC=fastqc

# variables for tardis and other apps
TARDIS_chunksize=300000
TARDIS_workdir=$(mytmp)


#*****************************************************************************************
# from here are the targets and rules
#*****************************************************************************************


###############################################
# top level target "logprecis" . This extracts from the log file the 
# relevant commands that were run, in a readable order 
# - run this after the actual build has been completed
# - note , you need to tell make to ignore errors
# for this target (for when grep fails to find a pattern)
# (make -i )
###############################################
%.logprecis: %.log
	echo "fastqc" > $*.logprecis
	echo "------" >> $*.logprecis
	egrep "^$(RUN_FASTQC)" $*.log >> $*.logprecis
	echo "fastq-join" >> $*.logprecis
	echo "----------" >> $*.logprecis
	egrep "^fastq-join" $*.log >> $*.logprecis
	echo "tardis" >> $*.logprecis
	echo "-------" >> $*.logprecis
	egrep "tardis" $*.log  | grep -v "logging this session" | egrep -v "^\#" >> $*.logprecis
	echo "setup" >> $*.logprecis
	echo "-----" >> $*.logprecis
	egrep "^ln\s+" $*.log | grep combined >> $*.logprecis
	echo "combine" >> $*.logprecis
	echo "-------" >> $*.logprecis
	egrep "^cat" $*.log | grep combined >> $*.logprecis
	echo "make_otu_table" >> $*.logprecis
	echo "--------------" >> $*.logprecis
	egrep "make_otu_table" $*.log  >> $*.logprecis
	echo "summarize_taxa" >> $*.logprecis
	echo "--------------" >> $*.logprecis
	egrep "summarize_taxa" $*.log  >> $*.logprecis

###############################################
# top level phony target "versions"  - output versions of all tools 
# - note , you need to tell make to ignore errors 
# for this target - for some tools , the only way to get a version
# string is to run them with options that provoke an error
# (make -i )
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
	echo "fastqc"  >> versions.log
	echo "------"  >> versions.log
	echo $(RUN_FASTQC) -version  >> versions.log
	$(RUN_FASTQC) -version  >> versions.log  2>&1
	echo rpm -q fastqc >> versions.log
	rpm -q fastqc >> versions.log 2>&1
	echo "fastq-join"  >> versions.log
	echo "------"  >> versions.log
	echo fastq-join  >> versions.log
	fastq-join  >> versions.log  2>&1
	echo rpm -q ea-utils >> versions.log
	rpm -q ea-utils >> versions.log 2>&1


###############################################
# how to make a taxonomy assignment table ( for all samples) based on R1 files
# NOTE this does NOT work due to a memory fault  - we need to make a table for each sample
###############################################
%R1tax_table: $(R1tax_assignments_fake_otu_map) $(R1tax_assignments_patched_with_otuid)
	# save the lists of files that will be concatenated so we can verify they are in the correct order
	# (they should be !)
	echo $(R1tax_assignments_fake_otu_map) > $(*D)/R1tax_table.fake_otu_map.filelist
	echo $(R1tax_assignments_patched_with_otuid) > $(*D)/R1tax_table.patched_with_otuid.filelist
	# make the single concatenated files , re-writing the otu id as we go
	cat $(R1tax_assignments_fake_otu_map) | awk -F "\t" '{printf("%s\t%s\n", NR-1,$$2);}' - > $(*D)/R1tax_table.fake_otu_map.tmp
	cat $(R1tax_assignments_patched_with_otuid) | awk -F "\t" '{printf("%s\t%s\t%s\t%s\n",NR-1,$$2,$$3,$$4);}' > $(*D)/R1tax_table.patched_with_otuid.tmp
	# run the summaries
	make_otu_table.py -i $(*D)/R1tax_table.fake_otu_map.tmp -t $(*D)/R1tax_table.patched_with_otuid.tmp -o $(*D)/R1tax_table.otu_table
	summarize_taxa.py  -i $(*D)/R1tax_table.otu_table -o $@/ -a -L 4,6,7,10

###############################################
# how to make a taxonomy assignment table (for all samples) based on R2 files
# NOTE this does NOT work due to a memory fault  - we need to make a table for each sample
###############################################
%R2tax_table: $(R2tax_assignments_fake_otu_map) $(R2tax_assignments_patched_with_otuid)
	# save the lists of files that will be concatenated so we can verify they are in the correct order
	# (they should be !)
	echo $(R2tax_assignments_fake_otu_map) > $(*D)/R2tax_table.fake_otu_map.filelist
	echo $(R2tax_assignments_patched_with_otuid) > $(*D)/R2tax_table.patched_with_otuid.filelist
	# make the single concatenated files , re-writing the otu id as we go
	cat $(R2tax_assignments_fake_otu_map) | awk -F "\t" '{printf("%s\t%s\n", NR-1,$$2);}' - > $(*D)/R2tax_table.fake_otu_map.tmp
	cat $(R2tax_assignments_patched_with_otuid) | awk -F "\t" '{printf("%s\t%s\t%s\t%s\n",NR-1,$$2,$$3,$$4);}' > $(*D)/R2tax_table.patched_with_otuid.tmp
	# run the summaries
	make_otu_table.py -i $(*D)/R2tax_table.fake_otu_map.tmp -t $(*D)/R2tax_table.patched_with_otuid.tmp -o $(*D)/R2tax_table.otu_table
	summarize_taxa.py  -i $(*D)/R2tax_table.otu_table -o $@/ -a -L 4,6,7,10

###############################################
# how to make a taxonomy assignment table (for all samples) based on sampled R1 files
# this DOES work
###############################################
%R1tax_table_sampled: $(R1tax_assignments_sampled_fake_otu_map) $(R1tax_assignments_sampled_patched_with_otuid)
	# save the lists of files that will be concatenated so we can verify they are in the correct order
	# (they should be !)
	echo $(R1tax_assignments_sampled_fake_otu_map) > $(*D)/R1tax_table_sampled.fake_otu_map.filelist
	echo $(R1tax_assignments_sampled_patched_with_otuid) > $(*D)/R1tax_table_sampled.patched_with_otuid.filelist
	# make the single concatenated files , re-writing the otu id as we go
	cat $(R1tax_assignments_sampled_fake_otu_map) | awk -F "\t" '{printf("%s\t%s\n", NR-1,$$2);}' - > $(*D)/R1tax_table_sampled.fake_otu_map.tmp
	cat $(R1tax_assignments_sampled_patched_with_otuid) | awk -F "\t" '{printf("%s\t%s\t%s\t%s\n",NR-1,$$2,$$3,$$4);}' > $(*D)/R1tax_table_sampled.patched_with_otuid.tmp
	# run the summaries        
	make_otu_table.py -i $(*D)/R1tax_table_sampled.fake_otu_map.tmp -t $(*D)/R1tax_table_sampled.patched_with_otuid.tmp -o $(*D)/R1tax_table_sampled.otu_table
	summarize_taxa.py  -i $(*D)/R1tax_table_sampled.otu_table -o $@/ -a -L 4,6,7,10 


###############################################
# how to make a taxonomy assignment summary table (for all samples) based on 454 data (where we are given the OTU map file, we don't make a fake one)
###############################################
%.454tax_summary_table: $(454tax_assignments)
	ln -f -s $(run454dir)/$(*F).otu_map $*.otu_map 
	make_otu_table.py -i $*.otu_map  -t $< -o $*.otu_table
	summarize_taxa.py  -i $*.otu_table -o $@/ -a -L 4,6,7,10


###############################################
# how to fake OTU's - this target is used by R1, R2 and joins but not 454, where the OTU's are given 
###############################################
%.fake_otu_map: %
	cat $*  | awk -F "\t" '{printf("%s\t%s^^^^%s\n",NR-1,"$(*F)",$$1);}' -  | sed 's/_/-/g' - | sed 's/\^\^\^\^/_/g' > $@

%.patched_with_otuid: %
	cat $*  | awk -F "\t" '{printf("%s\t%s\t%s\t%s\n",NR-1,$$2,$$3,$$4);}' -  > $@


###############################################
# how to make the R1 tax assignment summaries. This generates a summarey for each sample
# (which is what you need for large files)
###############################################
.PHONY : R1tax_assignment_summaries
R1tax_assignment_summaries: $(R1tax_assignment_summaries)

###############################################
# how to make the R2 tax assignment summaries. This generates a summarey for each sample
# (which is what you need for large files)
###############################################
.PHONY : R2tax_assignment_summaries
R2tax_assignment_summaries: $(R2tax_assignment_summaries)

###############################################
# how to make the joined pair tax assignment summaries. This generates a summarey for each sample
# (which is what you need for large files)
###############################################
.PHONY : joined_pair_tax_assignment_summaries
joined_pair_tax_assignment_summaries: $(joined_pair_tax_assignment_summaries)


###############################################
# how to make taxonomy assignments  based on the 454 files
###############################################
.PHONY : 454tax_assignments
454tax_assignments: $(454tax_assignments)

###############################################
# how to make taxonomy assignments  based on the R1 files
###############################################
.PHONY : R1tax_assignments
R1tax_assignments: $(R1tax_assignments)

###############################################
# how to make taxonomy assignments  based on the R2 files
###############################################
.PHONY : R2tax_assignments
R2tax_assignments: $(R2tax_assignments)

###############################################
# how to make taxonomy assignments based on the join files
###############################################
.PHONY : joined_pair_tax_assignments
joined_pair_tax_assignments: $(joined_pair_tax_assignment_monikers)

###############################################
# how to make taxonomy assignments  based on sampled R1 files
###############################################
.PHONY : R1tax_assignments_sampled
R1tax_assignments_sampled: $(R1tax_assignments_sampled)


###############################################
# how to make pair-joins 
###############################################
.PHONY : joined_pairs
joined_pairs: $(joined_pair_monikers)

##############################################
# how to make  a summary of an R1 / R2  tax assignment file 
# for smaller datasets, we make a single summary covering all the 
# samples, and this target is not neeeded then. However for 
# larger datasets we need to run the summary script for each sample
# (and then "manually" combine them)
##############################################
%.combined.tax_assignments.summary: %.combined.tax_assignments.fake_otu_map %.combined.tax_assignments.patched_with_otuid
	make_otu_table.py -i $(*D)/$(*F).combined.tax_assignments.fake_otu_map -t $(*D)/$(*F).combined.tax_assignments.patched_with_otuid -o $(*D)/$(*F).combined.otu_table
	summarize_taxa.py  -i $(*D)/$(*F).combined.otu_table  -o $@/ -a -L 4,6,7,10

##############################################
# how to make  a summary of a join-pair tax assignment file.
# its identical to the above except for the suffices 
##############################################
%.join.tax_assignments.summary: %.join.tax_assignments.fake_otu_map %.join.tax_assignments.patched_with_otuid
	make_otu_table.py -i $(*D)/$(*F).join.tax_assignments.fake_otu_map -t $(*D)/$(*F).join.tax_assignments.patched_with_otuid -o $(*D)/$(*F).join.otu_table
	summarize_taxa.py  -i $(*D)/$(*F).join.otu_table  -o $@/ -a -L 4,6,7,10


##############################################
# how to make  a tax assignment on a join file 
##############################################
%.join.tax_assignments: %.join
	#~/zeromq/tardis.py -w -c 1000 -d /dataset/Kittelmann_Buccal_Ill/scratch/tmp  /stash.local/qiime/1.8.0/bin/assign_taxonomy.py -i _condition_fastq2fasta_input_$*.join -m blast -t /dataset/Kittelmann_Buccal_Ill/active/bin/combined.taxonomy -b /dataset/Kittelmann_Buccal_Ill/active/bin/combined.fasta -o _condition_output_assign_tax_temp _condition_uncompressedtext_product_{}/$(*F).join_tax_assignments.txt,$*.join.tax_assignments  _condition_uncompressedtext_product_{}/$(*F).join_tax_assignments.log,$*.join.tax_assignments.log
	# these next version keeps the intermediate fasta files 
	~/zeromq/tardis.py -k -w -c 500 -d /dataset/Kittelmann_Buccal_Ill/scratch/tmp  -t /dataset/Kittelmann_Buccal_Ill/active/bin/job_template.txt /stash.local/qiime/1.8.0/bin/assign_taxonomy.py -i _condition_fastq2fasta_input_$*.join -m blast -t /dataset/Kittelmann_Buccal_Ill/active/bin/combined.taxonomy -b /dataset/Kittelmann_Buccal_Ill/active/bin/combined.fasta -o _condition_output_assign_tax_temp _condition_uncompressedtext_product_{}/$(*F).join_tax_assignments.txt,$*.join.tax_assignments  _condition_uncompressedtext_product_{}/$(*F).join_tax_assignments.log,$*.join.tax_assignments.log

##############################################
# how to make  a tax assignment on a 454  file
##############################################
%.454tax_assignments: 
	ln -f -s $(run454dir)/$(*F) $*
	# tried hpctype local due to compute farm problems  but for some unknown reason getting no blast hits. 
	#~/zeromq/tardis.py -w -c 200 -hpctype local -d /dataset/Kittelmann_Buccal_Ill/scratch/tmp  /stash.local/qiime/1.8.0/bin/assign_taxonomy.py -i _condition_fasta_input_$* -m blast -t /dataset/Kittelmann_Buccal_Ill/active/bin/combined.taxonomy -b /dataset/Kittelmann_Buccal_Ill/active/bin/combined.fasta -o _condition_output_assign_tax_temp _condition_uncompressedtext_product_{}/$(*F)_tax_assignments.txt,$*.454tax_assignments  _condition_uncompressedtext_product_{}/$(*F)_tax_assignments.log,$*.454tax_assignments.log
	# OK going back to compute farm , using special template 
	~/zeromq/tardis.py -w -c 200 -d /dataset/Kittelmann_Buccal_Ill/scratch/tmp  -t /dataset/Kittelmann_Buccal_Ill/active/bin/job_template.txt  /stash.local/qiime/1.8.0/bin/assign_taxonomy.py -i _condition_fasta_input_$* -m blast -t /dataset/Kittelmann_Buccal_Ill/active/bin/combined.taxonomy -b /dataset/Kittelmann_Buccal_Ill/active/bin/combined.fasta -o _condition_output_assign_tax_temp _condition_uncompressedtext_product_{}/$(*F)_tax_assignments.txt,$*.454tax_assignments  _condition_uncompressedtext_product_{}/$(*F)_tax_assignments.log,$*.454tax_assignments.log


##############################################
# how to make  a tax assignment on an R1 or R2 file
##############################################
%.combined.tax_assignments: %.combined 
	#~/zeromq/tardis.py -w -c 1000 -d /dataset/Kittelmann_Buccal_Ill/scratch/tmp  /stash.local/qiime/1.8.0/bin/assign_taxonomy.py -i _condition_fastq2fasta_input_$*.combined -m blast -t /dataset/Kittelmann_Buccal_Ill/active/bin/combined.taxonomy -b /dataset/Kittelmann_Buccal_Ill/active/bin/combined.fasta -o _condition_output_assign_tax_temp _condition_uncompressedtext_product_{}/$(*F).combined_tax_assignments.txt,$*.combined.tax_assignments  _condition_uncompressedtext_product_{}/$(*F).combined_tax_assignments.log,$*.combined.tax_assignments.log
	# this next version keeps the intermediate fasta files 
	#~/zeromq/tardis.py -k -w -c 500 -d /dataset/Kittelmann_Buccal_Ill/scratch/tmp  /stash.local/qiime/1.8.0/bin/assign_taxonomy.py -i _condition_fastq2fasta_input_$*.combined -m blast -t /dataset/Kittelmann_Buccal_Ill/active/bin/combined.taxonomy -b /dataset/Kittelmann_Buccal_Ill/active/bin/combined.fasta -o _condition_output_assign_tax_temp _condition_uncompressedtext_product_{}/$(*F).combined_tax_assignments.txt,$*.combined.tax_assignments  _condition_uncompressedtext_product_{}/$(*F).combined_tax_assignments.log,$*.combined.tax_assignments.log
	# this next version uses a custom template to exclude some problematic nodes 
	~/zeromq/tardis.py -k -w -c 500 -d /dataset/Kittelmann_Buccal_Ill/scratch/tmp -t /dataset/Kittelmann_Buccal_Ill/active/bin/job_template.txt  /stash.local/qiime/1.8.0/bin/assign_taxonomy.py -i _condition_fastq2fasta_input_$*.combined -m blast -t /dataset/Kittelmann_Buccal_Ill/active/bin/combined.taxonomy -b /dataset/Kittelmann_Buccal_Ill/active/bin/combined.fasta -o _condition_output_assign_tax_temp _condition_uncompressedtext_product_{}/$(*F).combined_tax_assignments.txt,$*.combined.tax_assignments  _condition_uncompressedtext_product_{}/$(*F).combined_tax_assignments.log,$*.combined.tax_assignments.log


##############################################
# how to make  a tax assignment on a sampled R1 or R2 file 
##############################################
%.combined.tax_assignments_sampled: %.combined 
	#~/zeromq/tardis.py -w -c 1000 -s .1 -d /dataset/Kittelmann_Buccal_Ill/scratch/tmp  /stash.local/qiime/1.8.0/bin/assign_taxonomy.py -i _condition_fastq2fasta_input_$*.combined -m blast -t /dataset/Kittelmann_Buccal_Ill/active/bin/combined.taxonomy -b /dataset/Kittelmann_Buccal_Ill/active/bin/combined.fasta -o _condition_output_assign_tax_temp _condition_uncompressedtext_product_{}/$(*F).combined_tax_assignments.txt,$*.combined.tax_assignments_sampled  _condition_uncompressedtext_product_{}/$(*F).combined_tax_assignments.log,$*.combined.tax_assignments_sampled.log
	~/zeromq/tardis.py -k -w -c 1000 -s .1 -d /dataset/Kittelmann_Buccal_Ill/scratch/tmp  /stash.local/qiime/1.8.0/bin/assign_taxonomy.py -i _condition_fastq2fasta_input_$*.combined -m blast -t /dataset/Kittelmann_Buccal_Ill/active/bin/combined.taxonomy -b /dataset/Kittelmann_Buccal_Ill/active/bin/combined.fasta -o _condition_output_assign_tax_temp _condition_uncompressedtext_product_{}/$(*F).combined_tax_assignments.txt,$*.combined.tax_assignments_sampled  _condition_uncompressedtext_product_{}/$(*F).combined_tax_assignments.log,$*.combined.tax_assignments_sampled.log


##############################################
# how to make  a joined pair 
##############################################
%L001_001.fastq.join: %L001_R1_001.fastq.trimmed.combined  %L001_R2_001.fastq.trimmed.combined
	fastq-join -v / $*L001_R1_001.fastq.trimmed.combined $*L001_R2_001.fastq.trimmed.combined -o  $*L001_001.fastq.un1 -o $*L001_001.fastq.un2 -o $*L001_001.fastq.join 


##############################################
# how to make  a combined file 
##############################################
%.combined: 
	cat $(run1dir)/$(*F) $(run2dir)/$(*F) > $@
	$(RUN_FASTQC) $@ -o $(builddir)


##############################################
# specify the intermediate files to keep 
##############################################
.PRECIOUS: %.combined %.join  

##############################################
# cleaning - not yet doing this using make  
##############################################
clean:
	echo "no clean for now" 

