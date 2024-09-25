#!/bin/sh
#
# Load the module
source /apps/profiles/modules_asax.sh.dyn
module load stringtie/2.2.1
#
#
# make a subdirectory where we are going to output the results file
# mkdir stringtie_results
# Now we will run the program stringtie, which will use our alignments (the sorted .bam files) to estimate expression levels for each gene in our genome

# The -p flag tells the program to run on 8 threads
# -e tells stringtie to operate in expression estimation mode, limiting the processing of read alignments to estimating the coverage of transcripts given in the -G option
# -G allows you to specify a gff file containing annotation information for the genome you mapped your reads to
# -B specifies the output format of tables containing coverage data for the reference transcripts in the GFF file
# for sample in ${samplelist}
# do
	# mkdir stringtie_results/${sample}
# stringtie -p 8 -e -B -G genomic.gff -o ~/stringtie_results/${sample}/${sample}.gtf ${sample}_sorted.bam
stringtie -p 8 -e -B -G Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gff3 -o ~/stringtie_results/L1C1/L1C1.gtf L1C1_sorted.bam
stringtie -p 8 -e -B -G Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gff3 -o /home/aubadc002/killifish/fastqc/paired/stringtie_results/L1C2.gtf L1C2_sorted.bam
stringtie -p 8 -e -B -G Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gff3 -o /home/aubadc002/killifish/fastqc/paired/stringtie_results/L1C3.gtf L1C3_sorted.bam
stringtie -p 8 -e -B -G Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gff3 -o /home/aubadc002/killifish/fastqc/paired/stringtie_results/L1C4.gtf L1C4_sorted.bam
stringtie -p 8 -e -B -G Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gff3 -o /home/aubadc002/killifish/fastqc/paired/stringtie_results/L1C5.gtf L1C5_sorted.bam
stringtie -p 8 -e -B -G Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gff3 -o /home/aubadc002/killifish/fastqc/paired/stringtie_results/L2U1.gtf L2U1_sorted.bam
stringtie -p 8 -e -B -G Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gff3 -o /home/aubadc002/killifish/fastqc/paired/stringtie_results/L2U2.gtf L2U2_sorted.bam
stringtie -p 8 -e -B -G Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gff3 -o /home/aubadc002/killifish/fastqc/paired/stringtie_results/L2U3.gtf L2U3_sorted.bam
stringtie -p 8 -e -B -G Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gff3 -o /home/aubadc002/killifish/fastqc/paired/stringtie_results/L2U4.gtf L2U4_sorted.bam
stringtie -p 8 -e -B -G Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gff3 -o /home/aubadc002/killifish/fastqc/paired/stringtie_results/L2U5.gtf L2U5_sorted.bam
stringtie -p 8 -e -B -G Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gff3 -o /home/aubadc002/killifish/fastqc/paired/stringtie_results/L3C1.gtf L3C1_sorted.bam
stringtie -p 8 -e -B -G Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gff3 -o /home/aubadc002/killifish/fastqc/paired/stringtie_results/L3C2.gtf L3C2_sorted.bam
stringtie -p 8 -e -B -G Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gff3 -o /home/aubadc002/killifish/fastqc/paired/stringtie_results/L3C3.gtf L3C3_sorted.bam
stringtie -p 8 -e -B -G Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gff3 -o /home/aubadc002/killifish/fastqc/paired/stringtie_results/L3C4.gtf L3C4_sorted.bam
stringtie -p 8 -e -B -G Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gff3 -o /home/aubadc002/killifish/fastqc/paired/stringtie_results/L3C5.gtf L3C5_sorted.bam
stringtie -p 8 -e -B -G Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gff3 -o /home/aubadc002/killifish/fastqc/paired/stringtie_results/L4U1.gtf L4U1_sorted.bam
stringtie -p 8 -e -B -G Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gff3 -o /home/aubadc002/killifish/fastqc/paired/stringtie_results/L4U2.gtf L4U2_sorted.bam
stringtie -p 8 -e -B -G Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gff3 -o /home/aubadc002/killifish/fastqc/paired/stringtie_results/L4U3.gtf L4U3_sorted.bam
stringtie -p 8 -e -B -G Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gff3 -o /home/aubadc002/killifish/fastqc/paired/stringtie_results/L4U4.gtf L4U4_sorted.bam
stringtie -p 8 -e -B -G Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gff3 -o /home/aubadc002/killifish/fastqc/paired/stringtie_results/L4U5.gtf L4U5_sorted.bam
stringtie -p 8 -e -B -G Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gff3 -o /home/aubadc002/killifish/fastqc/paired/stringtie_results/W1C1.gtf W1C1_sorted.bam
stringtie -p 8 -e -B -G Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gff3 -o /home/aubadc002/killifish/fastqc/paired/stringtie_results/W1C2.gtf W1C2_sorted.bam
stringtie -p 8 -e -B -G Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gff3 -o /home/aubadc002/killifish/fastqc/paired/stringtie_results/W1C3.gtf W1C3_sorted.bam
stringtie -p 8 -e -B -G Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gff3 -o /home/aubadc002/killifish/fastqc/paired/stringtie_results/W1C4.gtf W1C4_sorted.bam
stringtie -p 8 -e -B -G Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gff3 -o /home/aubadc002/killifish/fastqc/paired/stringtie_results/W2U1.gtf W2U1_sorted.bam
stringtie -p 8 -e -B -G Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gff3 -o /home/aubadc002/killifish/fastqc/paired/stringtie_results/W2U2.gtf W2U2_sorted.bam
stringtie -p 8 -e -B -G Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gff3 -o /home/aubadc002/killifish/fastqc/paired/stringtie_results/W2U3.gtf W2U3_sorted.bam
stringtie -p 8 -e -B -G Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gff3 -o /home/aubadc002/killifish/fastqc/paired/stringtie_results/W2U4.gtf W2U4_sorted.bam
stringtie -p 8 -e -B -G Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gff3 -o /home/aubadc002/killifish/fastqc/paired/stringtie_results/W2U5.gtf W2U5_sorted.bam
stringtie -p 8 -e -B -G Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gff3 -o /home/aubadc002/killifish/fastqc/paired/stringtie_results/W3C1.gtf W3C1_sorted.bam
stringtie -p 8 -e -B -G Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gff3 -o /home/aubadc002/killifish/fastqc/paired/stringtie_results/W3C2.gtf W3C2_sorted.bam
stringtie -p 8 -e -B -G Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gff3 -o /home/aubadc002/killifish/fastqc/paired/stringtie_results/W3C3.gtf W3C3_sorted.bam
stringtie -p 8 -e -B -G Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gff3 -o /home/aubadc002/killifish/fastqc/paired/stringtie_results/W3C4.gtf W3C4_sorted.bam
stringtie -p 8 -e -B -G Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gff3 -o /home/aubadc002/killifish/fastqc/paired/stringtie_results/W4U1.gtf W4U1_sorted.bam
stringtie -p 8 -e -B -G Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gff3 -o /home/aubadc002/killifish/fastqc/paired/stringtie_results/W4U2.gtf W4U2_sorted.bam
stringtie -p 8 -e -B -G Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gff3 -o /home/aubadc002/killifish/fastqc/paired/stringtie_results/W4U3.gtf W4U3_sorted.bam
stringtie -p 8 -e -B -G Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gff3 -o /home/aubadc002/killifish/fastqc/paired/stringtie_results/W4U4.gtf W4U4_sorted.bam

exit
