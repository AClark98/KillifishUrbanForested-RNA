#!/bin/bash

###Created Feb 9th, 2024
###This was done to look at the quality of the RNA data obtained from Novagene

###First we load de program fastqc
module load fastqc

##execute the program with all the fq.gz files (*.fq.gz)
##have all the output files redirected to a folder calles Results2

fastqc /home/aubaaf001/TSCAR/*.fq.gz --outdir=./Result2

exit
