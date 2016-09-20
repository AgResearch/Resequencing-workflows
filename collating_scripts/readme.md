This folder contains  scripting to process "raw"  BAM files from N animals 
to yield a single VCF file per chromosome that assays variants at each 
position based on evidence from all the BAMS.

1. gatk_run.sh sets up a link-farm with links to all of the BAM files to 
   process. It is given an "intervals" file contains a split of the 
   genome into ~ 5,000 overlapping intervals, and runs the tardis 
   utility to queue ~ 5,000 jobs onto the cluster, each processing
   1 genomic interval and all BAM files. Each of these jobs runs the 
   "gatk_wrapper.sh" script.

2. finalise_vcfs.sh contains methods to collate the ~ 5,000 VCF files
   covering overlapping genomic sub-regions from (1) , into a single 
   final filtered VCF file.

Each of these steps is run manually, in some cases uncommenting 
and commenting bash function calls in sequence - i.e. this is not 
yet an automated process. The scripts contain many hard-coded paths
as they are a snapshot of scripting used to run a specific VCF build.


