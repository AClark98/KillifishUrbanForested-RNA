################# GO_MWU for Gulf Killifish #######################


# Antrelle D. Clark, Auburn University Department of Biological Sciences, Bernal Lab, Fall 2024; adc0064@auburn.edu
# R-script version 2024.09.0+375

# Mock script from Wright, R. M., Aglyamova, G. V., Meyer, E. and Matz, M. V. Gene expression associated with white syndromes in a reef-building coral, Acropora hyacinthus. BMC Genomics 2015, 16: 371. ( http://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-015-1540-2 )


# GO_MWU uses continuous measure of significance (such as fold-change or -log(p-value)) to identify GO categories that are
# significantly enriches with either up- or down- regulated genes. The advantage - no need to impose arbitrary signifcance cutoff.

# If the measure is binary (0 or 1) the script will perform a typical "GO enrichment" analysis based FIsher's exact test: it will show GO categories over-represented among the genes that have 1 as their measure.
# On the plot, different fonts are used to indicate significance and color indicates enrichment with either up (red) or down (blue) regulated genes. No colors are shown for binary measure analysis.
# The tree on the plot is hierarchical clustering of GO categories based on shared genes. Categories with no branch length between them are subsets of each other.
# The fraction next to GO category name indicates the fraction of "good" genes in it; "good" genes being the ones exceeding the arbitrary absValye cutoff (option in gomwwuPlot). FOr Fisher's based test, specify absValue=0.5. THis value does not affect statistics and is used for plotting only.
# Stretch the plot manually to match tree to text

##################################################################
# First, insure that "perl" is on your local computer and installed properly (see "Perl_troubleshoot.txt" in Box folder for help) -AND- set working directory
system("C:/Windows/strawberry-perl-5.40.0.1-64bit-PDL/perl/bin/perl.exe -v") # should output a paragraph starting with "This is perl 5, version..."
setwd("/Users/antre/Box/Bernal_lab/Antrelle/GulfKillifish")

#Load the library
library(dplyr)
library(GO.db)
library(ape)


#### LongTide vs LongTideRunoff #### - "LTvsLTR_ELFC_AG.csv"
#### LongTide vs WeakleyTide #### - "LTvsWT_ELFC_AG.csv"
#### LongTide vs WeakleyTideRunoff #### - "LTvsWTR_ELFC_AG.csv"
#### Weakley Tide vs WeaklyTideRunoff #### - "WTvsWTR_ELFC_AG.csv"
#### LongTideRunoff vs WeakleyTideRunoff #### - "LTRvsWTR_ELFC_AG.csv"
#### LongTideRunoff vs WeakleyTide #### - "WTvsLTR_ELFC_AG.csv"


# the "MWU" tables that are automatically outputted to your designated folder are the ones you want to use to identify the first and last 20 categories that are up and down regulated between PWCs

#################################################################################




#### LongTide vs WeakleyTide #### - "LTvsWT_ELFC_AG.csv"
# Edit these to match your data file names: 
# MF - LTvsWT
setwd("C:/Users/antre/Box/Bernal_lab/Antrelle/GulfKillifish")
input="LTvsWT_MB.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either significant or not).
goAnnotations="MummichogGOMB.txt" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="MF" # either MF, or BP, or CC
source("gomwu.functions.R")

# ------------- Calculating stats
# It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previous runs first.

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="C:/Windows/strawberry-perl-5.40.0.1-64bit-PDL/perl/bin/perl.exe", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1, #this was originally set to 0.1 but kept getting error so had to relax the threshold # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,  #this was originally set to 5 but kept getting error so had to relax the threshold# a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           #  Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #  Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)

# do not continue if the printout shows that no GO terms pass 10% FDR.
#output: 24 GO terms at 10% FDR

# ----------- Plotting results

x11()  # use windows() or x11() for WINDOWS -and- quartz() for MACS
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #	absValue=1, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue. Originally was 0.05
                  level2=0.05, # FDR cutoff to print in regular (not italic) font. Originally was 0.01
                  level3=0.01, # FDR cutoff to print in large bold font. Originally was 0.001, but changed to 0.01
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001. Originally set at 0.1, 0.05, and 0.01 

# text representation of results, with actual adjusted p-values
results[[1]]


# BP - LTvsWT

setwd("C:/Users/antre/Box/Bernal_lab/Antrelle/GulfKillifish")
input="LTvsWT_MB.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either significant or not).
goAnnotations="MummichogGOMB.txt" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="BP" # either MF, or BP, or CC
source("gomwu.functions.R")

# ------------- Calculating stats
# It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previous runs first.

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="C:/Windows/strawberry-perl-5.40.0.1-64bit-PDL/perl/bin/perl.exe", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1, #this was originally set to 0.1 but kept getting error so had to relax the threshold # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,  #this was originally set to 5 but kept getting error so had to relax the threshold# a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           #  Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #  Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)

# do not continue if the printout shows that no GO terms pass 10% FDR.
#output: 43 GO terms at 10% FDR

# ----------- Plotting results

x11()  # use windows() or x11() for WINDOWS -and- quartz() for MACS
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #	absValue=1, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue. Originally was 0.05
                  level2=0.05, # FDR cutoff to print in regular (not italic) font. Originally was 0.01
                  level3=0.01, # FDR cutoff to print in large bold font. Originally was 0.001, but changed to 0.01
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001. Originally set at 0.1, 0.05, and 0.01 

# text representation of results, with actual adjusted p-values
results[[1]]


# CC - LTvsWT

setwd("C:/Users/antre/Box/Bernal_lab/Antrelle/GulfKillifish")
input="LTvsWT_MB.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either significant or not).
goAnnotations="MummichogGOMB.txt" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="CC" # either MF, or BP, or CC
source("gomwu.functions.R")

# ------------- Calculating stats
# It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previous runs first.

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="C:/Windows/strawberry-perl-5.40.0.1-64bit-PDL/perl/bin/perl.exe", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1, #this was originally set to 0.1 but kept getting error so had to relax the threshold # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,  #this was originally set to 5 but kept getting error so had to relax the threshold# a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           #  Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #  Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)

# do not continue if the printout shows that no GO terms pass 10% FDR.
#output: 34 GO terms at 10% FDR

# ----------- Plotting results

x11()  # use windows() or x11() for WINDOWS -and- quartz() for MACS
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #	absValue=1, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue. Originally was 0.05
                  level2=0.05, # FDR cutoff to print in regular (not italic) font. Originally was 0.01
                  level3=0.01, # FDR cutoff to print in large bold font. Originally was 0.001, but changed to 0.01
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001. Originally set at 0.1, 0.05, and 0.01 

# text representation of results, with actual adjusted p-values
results[[1]]


#################################################################################


#### LongTide vs LongTideRunoff #### - "LTvsLTR_ELFC_AG.csv"

# Edit these to match your data file names: 
# MF - LTvsLTR
setwd("C:/Users/antre/Box/Bernal_lab/Antrelle/GulfKillifish")
input="LTvsLTR_ELFC_AG.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either significant or not).
goAnnotations="MummichogGOMB.txt" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="MF" # either MF, or BP, or CC
source("gomwu.functions.R")

# ------------- Calculating stats
# It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previous runs first.

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="C:/Windows/strawberry-perl-5.40.0.1-64bit-PDL/perl/bin/perl.exe", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1, #this was originally set to 0.1 but kept getting error so had to relax the threshold # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,  #this was originally set to 5 but kept getting error so had to relax the threshold# a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           #  Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #  Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)

# do not continue if the printout shows that no GO terms pass 10% FDR.
#output: 85 GO terms at 10% FDR

# ----------- Plotting results

x11()  # use windows() or x11() for WINDOWS -and- quartz() for MACS
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #	absValue=1, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue. Originally was 0.05
                  level2=0.05, # FDR cutoff to print in regular (not italic) font. Originally was 0.01
                  level3=0.01, # FDR cutoff to print in large bold font. Originally was 0.001, but changed to 0.01
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001. Originally set at 0.1, 0.05, and 0.01 

# text representation of results, with actual adjusted p-values
results[[1]]



# BP - LTvsLTR
setwd("C:/Users/antre/Box/Bernal_lab/Antrelle/GulfKillifish")
input="LTvsLTR_ELFC_AG.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either significant or not).
goAnnotations="MummichogGOMB.txt" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="BP" # either MF, or BP, or CC
source("gomwu.functions.R")

# ------------- Calculating stats
# It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previous runs first.

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="C:/Windows/strawberry-perl-5.40.0.1-64bit-PDL/perl/bin/perl.exe", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1, #this was originally set to 0.1 but kept getting error so had to relax the threshold # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,  #this was originally set to 5 but kept getting error so had to relax the threshold# a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           #  Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #  Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)

# do not continue if the printout shows that no GO terms pass 10% FDR.
#output: 112 GO terms at 10% FDR

# ----------- Plotting results

x11()  # use windows() or x11() for WINDOWS -and- quartz() for MACS
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #	absValue=1, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue. Originally was 0.05
                  level2=0.05, # FDR cutoff to print in regular (not italic) font. Originally was 0.01
                  level3=0.01, # FDR cutoff to print in large bold font. Originally was 0.001, but changed to 0.01
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001. Originally set at 0.1, 0.05, and 0.01 

# text representation of results, with actual adjusted p-values
results[[1]]



# CC - LTvsLTR

setwd("C:/Users/antre/Box/Bernal_lab/Antrelle/GulfKillifish")
input="LTvsLTR_ELFC_AG.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either significant or not).
goAnnotations="MummichogGOMB.txt" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="CC" # either MF, or BP, or CC
source("gomwu.functions.R")

# ------------- Calculating stats
# It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previous runs first.

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="C:/Windows/strawberry-perl-5.40.0.1-64bit-PDL/perl/bin/perl.exe", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1, #this was originally set to 0.1 but kept getting error so had to relax the threshold # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,  #this was originally set to 5 but kept getting error so had to relax the threshold# a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           #  Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #  Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)

# do not continue if the printout shows that no GO terms pass 10% FDR.
#output: 101 GO terms at 10% FDR

# ----------- Plotting results

x11()  # use windows() or x11() for WINDOWS -and- quartz() for MACS
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #	absValue=1, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue. Originally was 0.05
                  level2=0.05, # FDR cutoff to print in regular (not italic) font. Originally was 0.01
                  level3=0.01, # FDR cutoff to print in large bold font. Originally was 0.001, but changed to 0.01
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001. Originally set at 0.1, 0.05, and 0.01 

# text representation of results, with actual adjusted p-values
results[[1]]



#################################################################################


#### LongTide vs WeakleyTideRunoff #### - "LTvsWTR_ELFC_AG.csv"

# MF - LTvsWTR
setwd("C:/Users/antre/Box/Bernal_lab/Antrelle/GulfKillifish")
input="LTvsWTR_ELFC_AG.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either significant or not).
goAnnotations="MummichogGOMB.txt" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="MF" # either MF, or BP, or CC
source("gomwu.functions.R")

# ------------- Calculating stats
# It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previous runs first.

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="C:/Windows/strawberry-perl-5.40.0.1-64bit-PDL/perl/bin/perl.exe", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1, #this was originally set to 0.1 but kept getting error so had to relax the threshold # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,  #this was originally set to 5 but kept getting error so had to relax the threshold# a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           #  Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #  Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)

# do not continue if the printout shows that no GO terms pass 10% FDR.
#output: 47 GO terms at 10% FDR

# ----------- Plotting results

x11()  # use windows() or x11() for WINDOWS -and- quartz() for MACS
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #	absValue=1, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue. Originally was 0.05
                  level2=0.05, # FDR cutoff to print in regular (not italic) font. Originally was 0.01
                  level3=0.01, # FDR cutoff to print in large bold font. Originally was 0.001, but changed to 0.01
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001. Originally set at 0.1, 0.05, and 0.01 

# text representation of results, with actual adjusted p-values
results[[1]]



# BP - LTvsWTR
setwd("C:/Users/antre/Box/Bernal_lab/Antrelle/GulfKillifish")
input="LTvsWTR_ELFC_AG.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either significant or not).
goAnnotations="MummichogGOMB.txt" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="BP" # either MF, or BP, or CC
source("gomwu.functions.R")

# ------------- Calculating stats
# It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previous runs first.

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="C:/Windows/strawberry-perl-5.40.0.1-64bit-PDL/perl/bin/perl.exe", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1, #this was originally set to 0.1 but kept getting error so had to relax the threshold # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,  #this was originally set to 5 but kept getting error so had to relax the threshold# a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           #  Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #  Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)

# do not continue if the printout shows that no GO terms pass 10% FDR.
#output: 93 GO terms at 10% FDR

# ----------- Plotting results

x11()  # use windows() or x11() for WINDOWS -and- quartz() for MACS
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #	absValue=1, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue. Originally was 0.05
                  level2=0.05, # FDR cutoff to print in regular (not italic) font. Originally was 0.01
                  level3=0.01, # FDR cutoff to print in large bold font. Originally was 0.001, but changed to 0.01
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001. Originally set at 0.1, 0.05, and 0.01 

# text representation of results, with actual adjusted p-values
results[[1]]



# CC - LTvsWTR
setwd("C:/Users/antre/Box/Bernal_lab/Antrelle/GulfKillifish")
input="LTvsWTR_ELFC_AG.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either significant or not).
goAnnotations="MummichogGOMB.txt" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="CC" # either MF, or BP, or CC
source("gomwu.functions.R")

# ------------- Calculating stats
# It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previous runs first.

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="C:/Windows/strawberry-perl-5.40.0.1-64bit-PDL/perl/bin/perl.exe", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1, #this was originally set to 0.1 but kept getting error so had to relax the threshold # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,  #this was originally set to 5 but kept getting error so had to relax the threshold# a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           #  Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #  Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)

# do not continue if the printout shows that no GO terms pass 10% FDR.
#output: 41 GO terms at 10% FDR

# ----------- Plotting results

x11()  # use windows() or x11() for WINDOWS -and- quartz() for MACS
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #	absValue=1, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue. Originally was 0.05
                  level2=0.05, # FDR cutoff to print in regular (not italic) font. Originally was 0.01
                  level3=0.01, # FDR cutoff to print in large bold font. Originally was 0.001, but changed to 0.01
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001. Originally set at 0.1, 0.05, and 0.01 

# text representation of results, with actual adjusted p-values
results[[1]]



#################################################################################

#### Weakley Tide vs WeaklyTideRunoff #### - "WTvsWTR_ELFC_AG.csv"

# MF - WTvsWTR
setwd("C:/Users/antre/Box/Bernal_lab/Antrelle/GulfKillifish")
input="WTvsWTR_ELFC_AG.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either significant or not).
goAnnotations="MummichogGOMB.txt" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="MF" # either MF, or BP, or CC
source("gomwu.functions.R")

# ------------- Calculating stats
# It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previous runs first.

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="C:/Windows/strawberry-perl-5.40.0.1-64bit-PDL/perl/bin/perl.exe", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1, #this was originally set to 0.1 but kept getting error so had to relax the threshold # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,  #this was originally set to 5 but kept getting error so had to relax the threshold# a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           #  Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #  Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)

# do not continue if the printout shows that no GO terms pass 10% FDR.
#output: 102 GO terms at 10% FDR

# ----------- Plotting results

x11()  # use windows() or x11() for WINDOWS -and- quartz() for MACS
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #	absValue=1, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue. Originally was 0.05
                  level2=0.05, # FDR cutoff to print in regular (not italic) font. Originally was 0.01
                  level3=0.01, # FDR cutoff to print in large bold font. Originally was 0.001, but changed to 0.01
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001. Originally set at 0.1, 0.05, and 0.01 

# text representation of results, with actual adjusted p-values
results[[1]]



# BP - WTvsWTR
setwd("C:/Users/antre/Box/Bernal_lab/Antrelle/GulfKillifish")
input="WTvsWTR_ELFC_AG.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either significant or not).
goAnnotations="MummichogGOMB.txt" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="BP" # either MF, or BP, or CC
source("gomwu.functions.R")

# ------------- Calculating stats
# It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previous runs first.

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="C:/Windows/strawberry-perl-5.40.0.1-64bit-PDL/perl/bin/perl.exe", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1, #this was originally set to 0.1 but kept getting error so had to relax the threshold # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,  #this was originally set to 5 but kept getting error so had to relax the threshold# a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           #  Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #  Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)

# do not continue if the printout shows that no GO terms pass 10% FDR.
#output: 172 GO terms at 10% FDR

# ----------- Plotting results

x11()  # use windows() or x11() for WINDOWS -and- quartz() for MACS
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #	absValue=1, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue. Originally was 0.05
                  level2=0.05, # FDR cutoff to print in regular (not italic) font. Originally was 0.01
                  level3=0.01, # FDR cutoff to print in large bold font. Originally was 0.001, but changed to 0.01
                  txtsize=0.6,    # 1.2 decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=1.0, # 0.5 height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001. Originally set at 0.1, 0.05, and 0.01 

# text representation of results, with actual adjusted p-values
results[[1]]



# CC - WTvsWTR
setwd("C:/Users/antre/Box/Bernal_lab/Antrelle/GulfKillifish")
input="WTvsWTR_ELFC_AG.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either significant or not).
goAnnotations="MummichogGOMB.txt" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="CC" # either MF, or BP, or CC
source("gomwu.functions.R")

# ------------- Calculating stats
# It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previous runs first.

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="C:/Windows/strawberry-perl-5.40.0.1-64bit-PDL/perl/bin/perl.exe", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1, #this was originally set to 0.1 but kept getting error so had to relax the threshold # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,  #this was originally set to 5 but kept getting error so had to relax the threshold# a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           #  Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #  Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)

# do not continue if the printout shows that no GO terms pass 10% FDR.
#output: 75 GO terms at 10% FDR

# ----------- Plotting results

x11()  # use windows() or x11() for WINDOWS -and- quartz() for MACS
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #	absValue=1, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue. Originally was 0.05
                  level2=0.05, # FDR cutoff to print in regular (not italic) font. Originally was 0.01
                  level3=0.01, # FDR cutoff to print in large bold font. Originally was 0.001, but changed to 0.01
                  txtsize=1.0,    # 1.2 decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # 0.5 height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001. Originally set at 0.1, 0.05, and 0.01 

# text representation of results, with actual adjusted p-values
results[[1]]



#################################################################################


#### LongTideRunoff vs WeakleyTideRunoff #### - "LTRvsWTR_ELFC_AG.csv"

# MF - LTRvsWTR
setwd("C:/Users/antre/Box/Bernal_lab/Antrelle/GulfKillifish")
input="LTRvsWTR_ELFC_AG.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either significant or not).
goAnnotations="MummichogGOMB.txt" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="MF" # either MF, or BP, or CC
source("gomwu.functions.R")

# ------------- Calculating stats
# It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previous runs first.

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="C:/Windows/strawberry-perl-5.40.0.1-64bit-PDL/perl/bin/perl.exe", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1, #this was originally set to 0.1 but kept getting error so had to relax the threshold # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,  #this was originally set to 5 but kept getting error so had to relax the threshold# a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           #  Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #  Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)

# do not continue if the printout shows that no GO terms pass 10% FDR.
#output: 143 GO terms at 10% FDR

# ----------- Plotting results

x11()  # use windows() or x11() for WINDOWS -and- quartz() for MACS
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #	absValue=1, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue. Originally was 0.05
                  level2=0.05, # FDR cutoff to print in regular (not italic) font. Originally was 0.01
                  level3=0.01, # FDR cutoff to print in large bold font. Originally was 0.001, but changed to 0.01
                  txtsize=0.9,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001. Originally set at 0.1, 0.05, and 0.01 

# text representation of results, with actual adjusted p-values
results[[1]]



# BP - LTRvsWTR
setwd("C:/Users/antre/Box/Bernal_lab/Antrelle/GulfKillifish")
input="LTRvsWTR_ELFC_AG.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either significant or not).
goAnnotations="MummichogGOMB.txt" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="BP" # either MF, or BP, or CC
source("gomwu.functions.R")

# ------------- Calculating stats
# It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previous runs first.

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="C:/Windows/strawberry-perl-5.40.0.1-64bit-PDL/perl/bin/perl.exe", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1, #this was originally set to 0.1 but kept getting error so had to relax the threshold # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,  #this was originally set to 5 but kept getting error so had to relax the threshold# a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           #  Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #  Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)

# do not continue if the printout shows that no GO terms pass 10% FDR.
#output: 195 GO terms at 10% FDR

# ----------- Plotting results

x11()  # use windows() or x11() for WINDOWS -and- quartz() for MACS
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #	absValue=1, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue. Originally was 0.05
                  level2=0.05, # FDR cutoff to print in regular (not italic) font. Originally was 0.01
                  level3=0.01, # FDR cutoff to print in large bold font. Originally was 0.001, but changed to 0.01
                  txtsize=0.4,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001. Originally set at 0.1, 0.05, and 0.01 

# text representation of results, with actual adjusted p-values
results[[1]]



# CC - LTRvsWTR
setwd("C:/Users/antre/Box/Bernal_lab/Antrelle/GulfKillifish")
input="LTRvsWTR_ELFC_AG.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either significant or not).
goAnnotations="MummichogGOMB.txt" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="CC" # either MF, or BP, or CC
source("gomwu.functions.R")

# ------------- Calculating stats
# It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previous runs first.

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="C:/Windows/strawberry-perl-5.40.0.1-64bit-PDL/perl/bin/perl.exe", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1, #this was originally set to 0.1 but kept getting error so had to relax the threshold # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,  #this was originally set to 5 but kept getting error so had to relax the threshold# a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           #  Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #  Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)

# do not continue if the printout shows that no GO terms pass 10% FDR.
#output: 108 GO terms at 10% FDR

# ----------- Plotting results

x11()  # use windows() or x11() for WINDOWS -and- quartz() for MACS
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #	absValue=1, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue. Originally was 0.05
                  level2=0.05, # FDR cutoff to print in regular (not italic) font. Originally was 0.01
                  level3=0.01, # FDR cutoff to print in large bold font. Originally was 0.001, but changed to 0.01
                  txtsize=0.9,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001. Originally set at 0.1, 0.05, and 0.01 

# text representation of results, with actual adjusted p-values
results[[1]]



#################################################################################


#### LongTideRunoff vs WeakleyTide #### - "LTRvsWT_ELFC_AG.csv"


# MF - LTRvsWT
setwd("C:/Users/antre/Box/Bernal_lab/Antrelle/GulfKillifish")
input="LTRvsWT_ELFC_AG.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either significant or not).
goAnnotations="MummichogGOMB.txt" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="MF" # either MF, or BP, or CC
source("gomwu.functions.R")

# ------------- Calculating stats
# It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previous runs first.

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="C:/Windows/strawberry-perl-5.40.0.1-64bit-PDL/perl/bin/perl.exe", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1, #this was originally set to 0.1 but kept getting error so had to relax the threshold # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,  #this was originally set to 5 but kept getting error so had to relax the threshold# a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           #  Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #  Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)

# do not continue if the printout shows that no GO terms pass 10% FDR.
#output: 51 GO terms at 10% FDR

# ----------- Plotting results

x11()  # use windows() or x11() for WINDOWS -and- quartz() for MACS
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #	absValue=1, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue. Originally was 0.05
                  level2=0.05, # FDR cutoff to print in regular (not italic) font. Originally was 0.01
                  level3=0.01, # FDR cutoff to print in large bold font. Originally was 0.001, but changed to 0.01
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001. Originally set at 0.1, 0.05, and 0.01 

# text representation of results, with actual adjusted p-values
results[[1]]



# BP - LTRvsWT
setwd("C:/Users/antre/Box/Bernal_lab/Antrelle/GulfKillifish")
input="LTRvsWT_ELFC_AG.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either significant or not).
goAnnotations="MummichogGOMB.txt" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="BP" # either MF, or BP, or CC
source("gomwu.functions.R")

# ------------- Calculating stats
# It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previous runs first.

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="C:/Windows/strawberry-perl-5.40.0.1-64bit-PDL/perl/bin/perl.exe", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1, #this was originally set to 0.1 but kept getting error so had to relax the threshold # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,  #this was originally set to 5 but kept getting error so had to relax the threshold# a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           #  Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #  Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)

# do not continue if the printout shows that no GO terms pass 10% FDR.
#output: 56 GO terms at 10% FDR

# ----------- Plotting results

x11()  # use windows() or x11() for WINDOWS -and- quartz() for MACS
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #	absValue=1, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue. Originally was 0.05
                  level2=0.05, # FDR cutoff to print in regular (not italic) font. Originally was 0.01
                  level3=0.01, # FDR cutoff to print in large bold font. Originally was 0.001, but changed to 0.01
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001. Originally set at 0.1, 0.05, and 0.01 

# text representation of results, with actual adjusted p-values
results[[1]]



# CC - LTRvsWT
setwd("C:/Users/antre/Box/Bernal_lab/Antrelle/GulfKillifish")
input="WTvsLTR_ELFC_AG.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either significant or not).
goAnnotations="MummichogGOMB.txt" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="CC" # either MF, or BP, or CC
source("gomwu.functions.R")

# ------------- Calculating stats
# It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previous runs first.

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="C:/Windows/strawberry-perl-5.40.0.1-64bit-PDL/perl/bin/perl.exe", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1, #this was originally set to 0.1 but kept getting error so had to relax the threshold # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,  #this was originally set to 5 but kept getting error so had to relax the threshold# a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           #  Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #  Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)

# do not continue if the printout shows that no GO terms pass 10% FDR.
#output: 57 GO terms at 10% FDR

# ----------- Plotting results

x11()  # use windows() or x11() for WINDOWS -and- quartz() for MACS
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #	absValue=1, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue. Originally was 0.05
                  level2=0.05, # FDR cutoff to print in regular (not italic) font. Originally was 0.01
                  level3=0.01, # FDR cutoff to print in large bold font. Originally was 0.001, but changed to 0.01
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001. Originally set at 0.1, 0.05, and 0.01 

# text representation of results, with actual adjusted p-values
results[[1]]






