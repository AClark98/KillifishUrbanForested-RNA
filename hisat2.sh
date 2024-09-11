#!/bin/sh


# Load the module
module load gcc/9.3.0
module load hisat2/2.2.1
module load samtools/Jun2022

# This first step involves locating splice sites and exons within your genome. It requires that your genome has an annotation, which must be in GTF format. 
# Most genomes that have come out recently will have their annotation information in GFF format, which is slightly different. You can convert between GFF and GTF format using the program gffread. There's great tutorials for this on the internet!

