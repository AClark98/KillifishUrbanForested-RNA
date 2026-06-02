#!/bin/bash

###Created Sept 9th, 2024
###This is used to look at the quality of the RNA data

###First we load  program fastqc
module load fastqc

##execute the program with all the fq.gz files (*.fq.gz)
##can choose to have them redirected to another folder by placing the directory or just leave it to keep in same folder

fastqc *fq.gz

exit
