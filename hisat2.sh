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

# Now that we have our index built, let's run HISAT2 to map the reads from each sample to our genome. 
# We'll loop through all of our samples like we did with the trimmomatic script. 

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
hisat2 -p 8 -k 3 -x Fhetero_genome -1 L1C1B1filtered_1paired.fq -2 L1C1B2filtered_2paired.fq | samtools view -@ 8 -Sbh > L1C1.bam
hisat2 -p 8 -k 3 -x Fhetero_genome -1 L1C2B1filtered_1paired.fq -2 L1C2B2filtered_2paired.fq | samtools view -@ 8 -Sbh > L1C2.bam
hisat2 -p 8 -k 3 -x Fhetero_genome -1 L1C3B1filtered_1paired.fq -2 L1C3B2filtered_2paired.fq | samtools view -@ 8 -Sbh > L1C3.bam
hisat2 -p 8 -k 3 -x Fhetero_genome -1 L1C4B1filtered_1paired.fq -2 L1C4B2filtered_2paired.fq | samtools view -@ 8 -Sbh > L1C4.bam
hisat2 -p 8 -k 3 -x Fhetero_genome -1 L1C5B1filtered_1paired.fq -2 L1C5B2filtered_2paired.fq | samtools view -@ 8 -Sbh > L1C5.bam
hisat2 -p 8 -k 3 -x Fhetero_genome -1 L2U1B21filtered_1paired.fq -2 L2U1B22filtered_2paired.fq | samtools view -@ 8 -Sbh > L2U1.bam
hisat2 -p 8 -k 3 -x Fhetero_genome -1 L2U2B1filtered_1paired.fq -2 L2U2B2filtered_2paired.fq | samtools view -@ 8 -Sbh > L2U2.bam
hisat2 -p 8 -k 3 -x Fhetero_genome -1 L2U3B1filtered_1paired.fq -2 L2U3B2filtered_2paired.fq | samtools view -@ 8 -Sbh > L2U3.bam
hisat2 -p 8 -k 3 -x Fhetero_genome -1 L2U4B1filtered_1paired.fq -2 L2U4B2filtered_2paired.fq | samtools view -@ 8 -Sbh > L2U4.bam
hisat2 -p 8 -k 3 -x Fhetero_genome -1 L2U5B1filtered_1paired.fq -2 L2U5B2filtered_2paired.fq | samtools view -@ 8 -Sbh > L2U5.bam
hisat2 -p 8 -k 3 -x Fhetero_genome -1 L3C1B21filtered_1paired.fq -2 L3C1B22filtered_2paired.fq | samtools view -@ 8 -Sbh > L3C1.bam
hisat2 -p 8 -k 3 -x Fhetero_genome -1 L3C2B1filtered_1paired.fq -2 L3C2B2filtered_2paired.fq | samtools view -@ 8 -Sbh > L3C2.bam
hisat2 -p 8 -k 3 -x Fhetero_genome -1 L3C3B1filtered_1paired.fq -2 L3C3B2filtered_2paired.fq | samtools view -@ 8 -Sbh > L3C3.bam
hisat2 -p 8 -k 3 -x Fhetero_genome -1 L3C4B1filtered_1paired.fq -2 L3C4B2filtered_2paired.fq | samtools view -@ 8 -Sbh > L3C4.bam
hisat2 -p 8 -k 3 -x Fhetero_genome -1 L3C5B1filtered_1paired.fq -2 L3C5B2filtered_2paired.fq | samtools view -@ 8 -Sbh > L3C5.bam
hisat2 -p 8 -k 3 -x Fhetero_genome -1 L4U1B1filtered_1paired.fq -2 L4U1B2filtered_2paired.fq | samtools view -@ 8 -Sbh > L4U1.bam
hisat2 -p 8 -k 3 -x Fhetero_genome -1 L4U2B1filtered_1paired.fq -2 L4U2B2filtered_2paired.fq | samtools view -@ 8 -Sbh > L4U2.bam
hisat2 -p 8 -k 3 -x Fhetero_genome -1 L4U3B1filtered_1paired.fq -2 L4U3B2filtered_2paired.fq | samtools view -@ 8 -Sbh > L4U3.bam
hisat2 -p 8 -k 3 -x Fhetero_genome -1 L4U4B1filtered_1paired.fq -2 L4U4B2filtered_2paired.fq | samtools view -@ 8 -Sbh > L4U4.bam
hisat2 -p 8 -k 3 -x Fhetero_genome -1 L4U5B1filtered_1paired.fq -2 L4U5B2filtered_2paired.fq | samtools view -@ 8 -Sbh > L4U5.bam
hisat2 -p 8 -k 3 -x Fhetero_genome -1 W1C1B1filtered_1paired.fq -2 W1C1B2filtered_2paired.fq | samtools view -@ 8 -Sbh > W1C1.bam
hisat2 -p 8 -k 3 -x Fhetero_genome -1 W1C2B1filtered_1paired.fq -2 W1C2B2filtered_2paired.fq | samtools view -@ 8 -Sbh > W1C2.bam
hisat2 -p 8 -k 3 -x Fhetero_genome -1 W1C3B21filtered_1paired.fq -2 W1C3B22filtered_2paired.fq | samtools view -@ 8 -Sbh > W1C3.bam
hisat2 -p 8 -k 3 -x Fhetero_genome -1 W1C4B1filtered_1paired.fq -2 W1C4B2filtered_2paired.fq | samtools view -@ 8 -Sbh > W1C4.bam
hisat2 -p 8 -k 3 -x Fhetero_genome -1 W2U1B1filtered_1paired.fq -2 W2U1B2filtered_2paired.fq | samtools view -@ 8 -Sbh > W2U1.bam
hisat2 -p 8 -k 3 -x Fhetero_genome -1 W2U2B1filtered_1paired.fq -2 W2U2B2filtered_2paired.fq | samtools view -@ 8 -Sbh > W2U2.bam
hisat2 -p 8 -k 3 -x Fhetero_genome -1 W2U3B1filtered_1paired.fq -2 W2U3B2filtered_2paired.fq | samtools view -@ 8 -Sbh > W2U3.bam
hisat2 -p 8 -k 3 -x Fhetero_genome -1 W2U4B1filtered_1paired.fq -2 W2U4B2filtered_2paired.fq | samtools view -@ 8 -Sbh > W2U4.bam
hisat2 -p 8 -k 3 -x Fhetero_genome -1 W2U5B1filtered_1paired.fq -2 W2U5B2filtered_2paired.fq | samtools view -@ 8 -Sbh > W2U5.bam
hisat2 -p 8 -k 3 -x Fhetero_genome -1 W3C1B1filtered_1paired.fq -2 W3C1B2filtered_2paired.fq | samtools view -@ 8 -Sbh > W3C1.bam
hisat2 -p 8 -k 3 -x Fhetero_genome -1 W3C2B1filtered_1paired.fq -2 W3C2B2filtered_2paired.fq | samtools view -@ 8 -Sbh > W3C2.bam
hisat2 -p 8 -k 3 -x Fhetero_genome -1 W3C3B1filtered_1paired.fq -2 W3C3B2filtered_2paired.fq | samtools view -@ 8 -Sbh > W3C3.bam
hisat2 -p 8 -k 3 -x Fhetero_genome -1 W3C4B1filtered_1paired.fq -2 W3C4B2filtered_2paired.fq | samtools view -@ 8 -Sbh > W3C4.bam
hisat2 -p 8 -k 3 -x Fhetero_genome -1 W4U1B1filtered_1paired.fq -2 W4U1B2filtered_2paired.fq | samtools view -@ 8 -Sbh > W4U1.bam
hisat2 -p 8 -k 3 -x Fhetero_genome -1 W4U2B1filtered_1paired.fq -2 W4U2B2filtered_2paired.fq | samtools view -@ 8 -Sbh > W4U2.bam
hisat2 -p 8 -k 3 -x Fhetero_genome -1 W4U3B1filtered_1paired.fq -2 W4U3B2filtered_2paired.fq | samtools view -@ 8 -Sbh > W4U3.bam
hisat2 -p 8 -k 3 -x Fhetero_genome -1 W4U4B1filtered_1paired.fq -2 W4U4B2filtered_2paired.fq | samtools view -@ 8 -Sbh > W4U4.bam

exit
