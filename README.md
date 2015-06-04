This repository contains a variety of make-based bioinformatics 
workflows.  

Usually foo.mk is accompanied by a script to drive it 

Usually the makefile includes targets which can be used to make 
a precis of the verbose log that is generated, and to list tool 
versions. Examples of a verbose log, its precis and a versions.log
are included.

The scripts and makefiles include hard-coded paths, and are quite
specific to particular tool choices - e.g. the makefile that uses
quadtrim needs to deal with the way quadtrim likes to organise
its output, which is quite different to the way flexbar organises 
its output.

The makefiles use a command-markup utility tardis. Markup is used 
to indicate the input and output files for a given command, which 
enables tardis to condition the command to run on a cluster.
(https://bitbucket.org/agr-bifo/tardis)

