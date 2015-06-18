#
# Makefile to archive build from Resequencing workflow version 5. 
#***************************************************************************************
# changes 
#***************************************************************************************
# 
#***************************************************************************************
# references:
#***************************************************************************************
# make: 
#     http://www.gnu.org/software/make/manual/make.html
#
#
#*****************************************************************************************
# examples  
#*****************************************************************************************
#
#
# ******************************************************************************************
# initialise project specific variables - these are usually overridden by command-line settings as above
# ******************************************************************************************
builddir=/not set
archivedir=/not set
tempdir=/not set


# ******************************************************************************************
# initialise lists of file targets 
# ******************************************************************************************
vcf_files := $(addprefix $(archivedir)/, $(notdir $(wildcard $(builddir)/*.vcf)))
samplemerged_bam_files := $(addprefix $(archivedir)/, $(notdir $(wildcard $(builddir)/*_samplemerged.bam*)))
realigned_bam_files := $(addprefix $(archivedir)/, $(notdir $(wildcard $(builddir)/*_realigned.bam*)))
coverage_files :=  $(addprefix $(archivedir)/, $(notdir $(wildcard $(builddir)/*.coverage.*)))
intervals_files :=  $(addprefix $(archivedir)/, $(notdir $(wildcard $(builddir)/*.intervals)))
fastqc_files := $(addprefix $(archivedir)/, $(notdir $(wildcard $(builddir)/*fastqc*)))
original_files_listing := $(archivedir)/original_files_listing.txt



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
%.archivelogprecis: %.archivelog
	echo "file copying" > $*.archivelogprecis
	echo "------------" >> $*.archivelogprecis
	egrep "^cp -p" $*.archivelog >> $*.archivelogprecis

	echo "file listing" >> $*.archivelogprecis
	echo "------------" >> $*.archivelogprecis
	egrep "^ls -slt" $*.archivelog >> $*.archivelogprecis

	echo "waiting" >> $*.archivelogprecis
	echo "-------" >> $*.archivelogprecis
	egrep "^while" $*.archivelog >> $*.archivelogprecis

###############################################
# top level phony target "all"  - this is the one thats usually used 
###############################################
.PHONY : %.all 
%.all: $(original_files_listing) $(vcf_files) $(samplemerged_bam_files) $(realigned_bam_files) $(coverage_files) $(intervals_files) $(fastqc_files) %.logging
	echo "making all targets"

###############################################
# how to make the listing of original files in the build dir
###############################################
$(archivedir)/original_files_listing.txt:
	ls -slt $(builddir) > $@

###############################################
# how to archive data files 
###############################################
$(archivedir)/%:  $(builddir)/%
	cp -pRL $< $@ 

###############################################
# how to archive log files
###############################################
.PHONY : %.logging
%.logging:
	cp -p $(tempdir)/*.log $(archivedir)
	cp -p $(tempdir)/*.mk $(archivedir)

##############################################
# specify the "intermediate" files to keep - i.e. actually these are the archived files
##############################################
.PRECIOUS: $(archivedir)/%  $(archivedir)/original_files_listing.txt



##############################################
# cleaning - not yet doing this using make  
##############################################
clean:
	echo "no clean for now" 

