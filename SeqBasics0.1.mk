##############################################################################################
# This makefie can be used to process either a single file, or a whole folder of files, through 
# one or more of Blast, FastQC. It will handle either compressed or uncompressed files 
# (compressed files are left compressed), and allows you to take a random sample of 
# the files (e.g. for a contamination check blast run). You can run one of the steps and 
# then later come back and run more steps, or restart a workflow part way through 
# (make will work out what is already done and not redo it)   
#
# This workflow is experimental !
#
##############################################################################################
# known bugs / limitations
##############################################################################################
#
# * haven't added the FastQC target yet
# * some settings such as tardis chunksize , sample rate and blast database should probably either be passed in or calculated, rather than
#   being set in here. 
# * should sanity check -j option to prevent accidentlaly launching a huge number of processes 
# * when restarting a workflow after a crash of some sort - need to watch out for a target that was only 
#   partly completed  - e.g. an incomplete output file (make will only require the target file to 
#   exist and have the correct date-stamp - it can't tell its not complete). 
#
##############################################################################################
# examples:
##############################################################################################
#
# * note that it is a good idea to do a dry run using the -n option as shown , before starting your actual run. You can check the 
# verbose output and make sure that the files you expected to be processed will be processed , in the way you expected )
# * its also a very good idea to capture stdout and stderr from the run as shown
#
# blast a single file
# ===================
#
# --first create your output folder. Its usually best to use a new empty folder for 
# --each run of this workflow (because its easy to accidentllay overwrite an existing file otherwise) 
# mkdir /dataset/mydata/scratch/blastruns
#
# --The basic idea is that you "make" a target of myfile.blastout, in your new empty folder , and 
# --also tell make where to find "myfile"
#
# --dry run : 
# make -n -f SeqBasics0.1.mk -d --no-builtin-rules data_source_dir=/dataset/mydata/archive/140627_M02810_0023_000000000-A856J/processed_trimmed  /dataset/mydata/scratch/blastruns/processed_S395Buccal-PG-25MG_S18_L001_R1_001.fastq.trimmed.blastout  > SeqBasics.log 2>&1
# --actual run : 
# make -f SeqBasics0.1.mk -d --no-builtin-rules data_source_dir=/dataset/mydata/archive/140627_M02810_0023_000000000-A856J/processed_trimmed  /dataset/mydata/scratch/blastruns/processed_S395Buccal-PG-25MG_S18_L001_R1_001.fastq.trimmed.blastout  > SeqBasics.log 2>&1
#
# blast all files in a folder
# ==========================
#
# --first create your output folder. Its usually best to use a new empty folder for
# --each run of this workflow
# mkdir /dataset/mydata/scratch/blastruns
#
# --The basic idea is that you "make" a ("phony") target of "all_in_folder", tell make where the output should go, and  
# --also tell make where to find all the files to process 
#
# --(note we use -j 8 to tell make to start processing up to 8 files concurrently. Each of these 8 will in turn be split by tardis to run on the 
# --compute farm)
# --For this target we need to specify the outputdir (for the single file target it is implicit)
# --dry run :
# make -n -f SeqBasics0.1.mk -d --no-builtin-rules -j 8 data_source_dir=/dataset/mydata/archive/140627_M02810_0023_000000000-A856J/processed_trimmed  outputdir=/dataset/mydata/scratch/blastruns all_in_folder  > SeqBasics.log 2>&1
# --actual run :
# make -f SeqBasics0.1.mk -d --no-builtin-rules -j 8 data_source_dir=/dataset/mydata/archive/140627_M02810_0023_000000000-A856J/processed_trimmed  outputdir=/dataset/mydata/scratch/blastruns all_in_folder  > SeqBasics.log 2>&1

#############################################################################################
# parameters  
#############################################################################################
tardis_chunk_size=500
blast_db=/dataset/mydata/active/bin/combined.fasta 
other_tardis_options=
# example
# to process a 1/100 random sample of input 
# tardis_options=-s .001


#############################################################################################
# target lists (to support processing all files in a folder -see the "all_in_folder" phony target below)
# for the given folder , we will process 
# *.fastq
# *.fastq.gz
# *.fastq.trimmed
# *.fastq.trimmed.gz
#
# - to add additional patterns , just add additional variables as below, and include
# these in the rule dependency as below
#
# if you want to process files from multiple folders, the easiest way
# is probably to create a "shortcut-farm-folder" - i.e. make a new
# empty folder somewhere, create links to all the files you want to process, 
# and then use this as your source folder. (This approach could also be used 
# to process a subset of files in a folder)  
#############################################################################################
uncompressed_fastq_blast_targets=$(addprefix $(outputdir)/, $(addsuffix .blastout, $(notdir $(wildcard $(data_source_dir)/*.fastq))))
# e.g. from a source file "processed_Undetermined_S0_L001_R2_001.fastq", this would generate a target (say)  
# /dataset/mydata/scratch/blastruns/processed_Undetermined_S0_L001_R2_001.fastq.blastout
uncompressed_trimmed_fastq_blast_targets=$(addprefix $(outputdir)/, $(addsuffix .blastout, $(notdir $(wildcard $(data_source_dir)/*.fastq.trimmed))))
compressed_fastq_blast_targets=$(addprefix $(outputdir)/, $(addsuffix .blastout, $(notdir $(wildcard $(data_source_dir)/*.fastq.gz))))
compressed_trimmed_fastq_blast_targets=$(addprefix $(outputdir)/, $(addsuffix .blastout, $(notdir $(wildcard $(data_source_dir)/*.fastq.trimmed.gz))))

.PHONY : all_in_folder
all_in_folder: $(uncompressed_fastq_blast_targets) $(compressed_fastq_blast_targets) $(uncompressed_trimmed_fastq_blast_targets) $(compressed_trimmed_fastq_blast_targets)


%.blastout: 
	~/zeromq/tardis.py -w -c $(tardis_chunk_size) -d $(*D) $(other_tardis_options) blastn -query _condition_fastq2fasta_input_$(data_source_dir)/$(*F) -task blastn -num_threads 2 -db $(blast_db) -evalue 1.0e-10 -outfmt \"7 qseqid sseqid sacc stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore\" -max_target_seqs 1 -out _condition_text_output_$@ \> _condition_text_output_$@.stdout  2>&1

##############################################
# specify the intermediate files to keep 
##############################################
.PRECIOUS: %.blastout   

##############################################
# cleaning - not yet doing this using make  
##############################################
clean:
	echo "no clean for now" 
