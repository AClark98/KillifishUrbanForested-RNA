#!/bin/sh
#
# Load the module
#load the GNU compiling software environment
source /apps/profiles/modules_asax.sh.dyn
module load gcc/9.3.0
what is the command line?

#!/bin/sh
#
#load the HISAT2 alignment program for mapping ngs environment
module load hisat2/2.2.0
what is the command line?

#!/bin/sh
#
#load the SAMtools environment
source /apps/profiles/modules_asax.sh.dyn
module load samtools/1.18
what is the command line?

# This first step involves locating splice sites and exons within your genome. It requires that your genome has an annotation, which must be in GTF format. 
# Most genomes that have come out recently will have their annotation information in GFF format, which is slightly different. You can convert between GFF and GTF format using the program gffread. There's great tutorials for this on the internet!

# Identify your splice sites by running:
# hisat2_extract_splice_sites.py F_heteroclitus_genannotation.gtf > splice_sites.hisat

# Identify your exons by running:
# hisat2_extract_exons.py F_heteroclitus_genannotation.gtf > exons.hisat

# Now that we have our splice sites and exons, we can use that information to help build a splice site-informed index. This means that the program will take into account splice sites when it is mapping transcripts to the (unspliced) whole genome. 
# Run the following code to build your index:
# The -p flag specifies we want to run this code on 8 threads
# The -ss flag is used to specify the file we just created that contains the known splice sites in our genome
# The --exon flag is used to specify the file we just created that contains the known exons in our genome
# Then, we just list the file containing our genome (in fasta format), and the desired "base name" for our index. The program is going to create several files, but they'll all start with this name. 
# hisat2-build -p 8 -ss splice_sites.hisat --exon exons.hisat Fhetero_genome.fa Fhetero_genome
