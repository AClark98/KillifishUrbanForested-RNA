#!/bin/sh
#
# Load the module
source /apps/profiles/modules_asax.sh.dyn
module load samtools/1.18
# Delete the unsorted bam file so it's not taking up unnecessary space
	# rm ${sample}.bam
# !!!!!deleting the unsorted bam file was not done because we need the unsorted bam file to identify SNPs within the data!!!!!
#
# Index our sorted bam file, which is necessary for the next program we want to use
samtools index L1C1_sorted.bam
samtools index L1C2_sorted.bam
samtools index L1C3_sorted.bam
samtools index L1C4_sorted.bam
samtools index L1C5_sorted.bam
samtools index L2U1_sorted.bam
samtools index L2U2_sorted.bam
samtools index L2U3_sorted.bam
samtools index L2U4_sorted.bam
samtools index L2U5_sorted.bam
samtools index L3C1_sorted.bam
samtools index L3C2_sorted.bam
samtools index L3C3_sorted.bam
samtools index L3C4_sorted.bam
samtools index L3C5_sorted.bam
samtools index L4U1_sorted.bam
samtools index L4U2_sorted.bam
samtools index L4U3_sorted.bam
samtools index L4U4_sorted.bam
samtools index L4U5_sorted.bam
samtools index W1C1_sorted.bam
samtools index W1C2_sorted.bam
samtools index W1C3_sorted.bam
samtools index W1C4_sorted.bam
samtools index W2U1_sorted.bam
samtools index W2U2_sorted.bam
samtools index W2U3_sorted.bam
samtools index W2U4_sorted.bam
samtools index W2U5_sorted.bam
samtools index W3C1_sorted.bam
samtools index W3C2_sorted.bam
samtools index W3C3_sorted.bam
samtools index W3C4_sorted.bam
samtools index W4U1_sorted.bam
samtools index W4U2_sorted.bam
samtools index W4U3_sorted.bam
samtools index W4U4_sorted.bam

exit
