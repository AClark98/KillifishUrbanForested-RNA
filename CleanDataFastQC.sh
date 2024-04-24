#! /bin/bash

###This script will run Fastqc on the sequences that have been trimmed

###First we load de program fastqc
source /apps/profiles/modules_asax.sh.dyn
module load fastqc

##Make variables to find the data

CD=/home/aubaaf001/TSCAR/Trimmomatic
PCQ=PostCleanQuality

mkdir -p ${CD}/${PCQ}
cd ${CD}

fastqc TSCARfiltered_1paired.fq.gz --outdir=${CD}/${PCQ}
fastqc TSCARfiltered_2paired.fq.gz --outdir=${CD}/${PCQ}

echo "Fastqc analysis with clead data done!"

##move to the directory where the fastqc results are in and make a tarball
cd ${CD}/${PCQ}

##Make everything into a tarball

tar cvzf ${PCQ}.tar.gz *
echo "tarball done"
exit
