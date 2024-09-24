
#!/bin/sh
#
# Load the module
source /apps/profiles/modules_asax.sh.dyn
module load samtools/1.18
samtools sort -@ 8 -o L1C1_sorted.bam L1C1.bam
samtools sort -@ 8 -o L1C2_sorted.bam L1C2.bam
samtools sort -@ 8 -o L1C3_sorted.bam L1C3.bam
samtools sort -@ 8 -o L1C4_sorted.bam L1C4.bam
samtools sort -@ 8 -o L1C5_sorted.bam L1C5.bam
samtools sort -@ 8 -o L2U1_sorted.bam L2U1.bam
samtools sort -@ 8 -o L2U2_sorted.bam L2U2.bam
samtools sort -@ 8 -o L2U3_sorted.bam L2U3.bam
samtools sort -@ 8 -o L2U4_sorted.bam L2U4.bam
samtools sort -@ 8 -o L2U5_sorted.bam L2U5.bam
samtools sort -@ 8 -o L3C1_sorted.bam L3C1.bam
samtools sort -@ 8 -o L3C2_sorted.bam L3C2.bam
samtools sort -@ 8 -o L3C3_sorted.bam L3C3.bam
samtools sort -@ 8 -o L3C4_sorted.bam L3C4.bam
samtools sort -@ 8 -o L3C5_sorted.bam L3C5.bam
samtools sort -@ 8 -o L4U1_sorted.bam L4U1.bam
samtools sort -@ 8 -o L4U2_sorted.bam L4U2.bam
samtools sort -@ 8 -o L4U3_sorted.bam L4U3.bam
samtools sort -@ 8 -o L4U4_sorted.bam L4U4.bam
samtools sort -@ 8 -o L4U5_sorted.bam L4U5.bam
samtools sort -@ 8 -o W1C1_sorted.bam W1C1.bam
samtools sort -@ 8 -o W1C2_sorted.bam W1C2.bam
samtools sort -@ 8 -o W1C3_sorted.bam W1C3.bam
samtools sort -@ 8 -o W1C4_sorted.bam W1C4.bam
samtools sort -@ 8 -o W2U1_sorted.bam W2U1.bam
samtools sort -@ 8 -o W2U2_sorted.bam W2U2.bam
samtools sort -@ 8 -o W2U3_sorted.bam W2U3.bam
samtools sort -@ 8 -o W2U4_sorted.bam W2U4.bam
samtools sort -@ 8 -o W2U5_sorted.bam W2U5.bam
samtools sort -@ 8 -o W3C1_sorted.bam W3C1.bam
samtools sort -@ 8 -o W3C2_sorted.bam W3C2.bam
samtools sort -@ 8 -o W3C3_sorted.bam W3C3.bam
samtools sort -@ 8 -o W3C4_sorted.bam W3C4.bam
samtools sort -@ 8 -o W4U1_sorted.bam W4U1.bam
samtools sort -@ 8 -o W4U2_sorted.bam W4U2.bam
samtools sort -@ 8 -o W4U3_sorted.bam W4U3.bam
samtools sort -@ 8 -o W4U4_sorted.bam W4U4.bam


exit
