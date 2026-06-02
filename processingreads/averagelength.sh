
#!/bin/sh
#
# Load the module
source /apps/profiles/modules_asax.sh.dyn
module load samtools/1.18
#this code is used to identify the average length of the samples, which is needed for the prepDE.py script to format the read count data for DESeq2
samtools stats L1C1.bam | grep "average length"
samtools stats L1C2.bam | grep "average length"
samtools stats L1C3.bam | grep "average length"
samtools stats L1C4.bam | grep "average length"
samtools stats L1C5.bam | grep "average length"
samtools stats L2U1.bam | grep "average length"
samtools stats L2U2.bam | grep "average length"
samtools stats L2U3.bam | grep "average length"
samtools stats L2U4.bam | grep "average length"
samtools stats L2U5.bam | grep "average length"
samtools stats L3C1.bam | grep "average length"
samtools stats L3C2.bam | grep "average length"
samtools stats L3C3.bam | grep "average length"
samtools stats L3C4.bam | grep "average length"
samtools stats L3C5.bam | grep "average length"
samtools stats L4U1.bam | grep "average length"
samtools stats L4U2.bam | grep "average length"
samtools stats L4U3.bam | grep "average length"
samtools stats L4U4.bam | grep "average length"
samtools stats L4U5.bam | grep "average length"
samtools stats W1C1.bam | grep "average length"
samtools stats W1C2.bam | grep "average length"
samtools stats W1C3.bam | grep "average length"
samtools stats W1C4.bam | grep "average length"
samtools stats W2U1.bam | grep "average length"
samtools stats W2U2.bam | grep "average length"
samtools stats W2U3.bam | grep "average length"
samtools stats W2U4.bam | grep "average length"
samtools stats W2U5.bam | grep "average length"
samtools stats W3C1.bam | grep "average length"
samtools stats W3C2.bam | grep "average length"
samtools stats W3C3.bam | grep "average length"
samtools stats W3C4.bam | grep "average length"
samtools stats W4U1.bam | grep "average length"
samtools stats W4U2.bam | grep "average length"
samtools stats W4U3.bam | grep "average length"
samtools stats W4U4.bam | grep "average length"

exit
