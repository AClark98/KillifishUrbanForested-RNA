#!/bin/sh
#
# Load the module
source /apps/profiles/modules_asax.sh.dyn
# module load hisat2/2.2.0
# module load samtools/1.18
#^save as hisat2.sh

# This first step involves locating splice sites and exons within your genome. It requires that your genome has an annotation, which must be in GTF format. 
# Most genomes that have come out recently will have their annotation information in GFF format, which is slightly different. You can convert between GFF and GTF format using the program gffread. There's great tutorials for this on the internet!

#First, place "module load hisat2/2.2.0" into the terminal 
# Identify your splice sites by running:
# hisat2_extract_splice_sites.py Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gtf > splice_sites.hisat

# Identify your exons by running:
# hisat2_extract_exons.py Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gtf > exons.hisat

# Now that we have our splice sites and exons, we can use that information to help build a splice site-informed index. This means that the program will take into account splice sites when it is mapping transcripts to the (unspliced) whole genome. 
# Run the following code to build your index:
# The -p flag specifies we want to run this code on 8 threads
# The -ss flag is used to specify the file we just created that contains the known splice sites in our genome
# The --exon flag is used to specify the file we just created that contains the known exons in our genome
# Then, we just list the file containing our genome (in fasta format), and the desired "base name" for our index. The program is going to create several files, but they'll all start with this name. 
# hisat2-build -p 8 --ss splice_sites.hisat --exon exons.hisat Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gtf Fhetero_genome

# For each sample in our list:
      # We'll run hisat2 to align our cleaned reads to the F.heteroclitus genome
	    # -p tells the program to run on 8 threads
	    # -k tells the program to only search for a maximum of 3 alignments per read. We're only going to continue our analysis with reads that mapped uniquely (i.e., 1 time), so this makes sure the program doesn't spend a bunch of extra time aligning multi-mapping reads
	    # -x is the flag that we use to specify the "base name" of our index files
	    # -1 is the flag that we use to specify the fastq file containing our paired1 reads
	    # -2 is the flag that we use to specify the fastq file containing our paired2 reads
	    # After we've created our alignments, we're going to pipe the output of that command to the program samtools, which will convert our output from a .SAM format to a binary .BAM format. The file will encode the same information, but it takes up less space. 
     
