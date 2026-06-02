################# GO_MWU for Gulf Killifish #######################
# This script follows directly after the PWCtablemodify.R script


# Antrelle D. Clark, Auburn University Department of Biological Sciences, Bernal Lab, Fall 2024; adc0064@auburn.edu
# R-script version 2024.09.0+375


# This is a mock GO R script based off of the works of Ally Swank, PhD student @ Boston University: Abudefduf-Brains GO_MWU.R (), Adam H. De-novo Transcriptome Assembly (https://github.com/adhallaj/Denovo-Transcriptome-Assembly/blob/main/README.md) with help from Logan Turner, PhD student @ Auburn University
# https://github.com/z0on/GO_MWU
# Mikhail V. Matz, UT Austin, February 2015; matz@utexas.edu


# GO_MWU uses continuous measure of significance (such as fold-change or -log(p-value)) to identify GO categories that are
# significantly enriches with either up- or down- regulated genes. The advantage - no need to impose arbitrary signifcance cutoff.

# If the measure is binary (0 or 1) the script will perform a typical "GO enrichment" analysis based FIsher's exact test: it will show GO categories over-represented among the genes that have 1 as their measure.
# On the plot, different fonts are used to indicate significance and color indicates enrichment with either up (red) or down (blue) regulated genes. No colors are shown for binary measure analysis.
# The tree on the plot is hierarchical clustering of GO categories based on shared genes. Categories with no branch length between them are subsets of each other.
# The fraction next to GO category name indicates the fraction of "good" genes in it; "good" genes being the ones exceeding the arbitrary absValye cutoff (option in gomwwuPlot). FOr Fisher's based test, specify absValue=0.5. THis value does not affect statistics and is used for plotting only.
# Stretch the plot manually to match tree to text


#Load the library
library(clusterProfiler)
library(dplyr)

# Set working directory
setwd("/Users/antre/Box/Bernal_lab/Antrelle/GulfKillifish")

# For each pairwise comparison, remove columns that are not EnsemblID and Log2FoldChange


######################## FOR SIGNIFICANT GENES #########################


#### LongTide vs LongTideRunoff ####
LTvsLTR_pairwisefiltered <- read.csv("C:\\Users\\antre\\Box\\Bernal_lab\\Antrelle\\GulfKillifish\\LongTide_vs_LongTideRunoff_pairwisejoined_filtered.csv", header = TRUE, sep = ",", fill = TRUE)
colnames(LTvsLTR_pairwisefiltered)
  # [1] "X"              "EnsemblID"      "log2FoldChange" "padj"           "Gene_name"      "GO_ID"
#Subset the data to only keep the columns we want: EnsemblID and log2FoldChange
LTvsLTR_EnsemblLogFC <- LTvsLTR_pairwisefiltered[c("EnsemblID", "log2FoldChange")]
head(LTvsLTR_EnsemblLogFC) # "head" to ensure data is there; "colnames" just to see column names
  #EnsemblID log2FoldChange
  #1 ENSFHEG00000001645     -0.9454084
  #2 ENSFHEG00000001645     -0.9454084
  #3 ENSFHEG00000001645     -0.9454084
  #4 ENSFHEG00000001645     -0.9454084
  #5 ENSFHEG00000001645     -0.9454084
  #6 ENSFHEG00000001645     -0.9454084
#Filter data to remove repeated lines
LTvsLTR_Ensembl_LogFC <- LTvsLTR_EnsemblLogFC %>%
  distinct()
  #EnsemblID log2FoldChange
  #1 ENSFHEG00000001645     -0.9454084
  #2 ENSFHEG00000004018      0.7102799
  #3 ENSFHEG00000004541      0.7315790
  #4 ENSFHEG00000008149     -3.2127419
  #5 ENSFHEG00000008393     -0.3894052
  #6 ENSFHEG00000010807      0.8235978
# Extract table as a .csv file
write.table(LTvsLTR_Ensembl_LogFC, "LTvsLTR_Ensembl_LogFC.csv", sep = ",", quote = F)


#### LongTide vs WeakleyTide ####
LTvsWT_pairwisefiltered <- read.csv("C:\\Users\\antre\\Box\\Bernal_lab\\Antrelle\\GulfKillifish\\LongTide_vs_WeakleyTide_pairwisejoined_filtered.csv", header = TRUE, sep = ",", fill = TRUE)
colnames(LTvsWT_pairwisefiltered)
  # [1] "X"              "EnsemblID"      "log2FoldChange" "padj"           "Gene_name"      "GO_ID"  
#Subset the data to only keep the columns we want: EnsemblID and log2FoldChange
LTvsWT_EnsemblLogFC <- LTvsWT_pairwisefiltered[c("EnsemblID", "log2FoldChange")]
#Filter data to remove repeated lines
LTvsWT_Ensembl_LogFC <- LTvsWT_EnsemblLogFC %>%
  distinct()
# Extract table as a .csv file
write.table(LTvsWT_Ensembl_LogFC, "LTvsWT_Ensembl_LogFC.csv", sep = ",", quote = F)


#### LongTide vs WeakleyTideRunoff ####
LTvsWTR_pairwisefiltered <- read.csv("C:\\Users\\antre\\Box\\Bernal_lab\\Antrelle\\GulfKillifish\\LongTide_vs_WeakleyTideRunoff_pairwisejoined_filtered.csv", header = TRUE, sep = ",", fill = TRUE)
colnames(LTvsWTR_pairwisefiltered)
  # [1] "X"              "EnsemblID"      "log2FoldChange" "padj"           "Gene_name"      "GO_ID
#Subset the data to only keep the columns we want: EnsemblID and log2FoldChange
LTvsWTR_EnsemblLogFC <- LTvsWTR_pairwisefiltered[c("EnsemblID", "log2FoldChange")]
#Filter data to remove repeated lines
LTvsWTR_Ensembl_LogFC <- LTvsWTR_EnsemblLogFC %>%
  distinct()
# Extract table as a .csv file
write.table(LTvsWTR_Ensembl_LogFC, "LTvsWTR_Ensembl_LogFC.csv", sep = ",", quote = F)


#### Weakley Tide vs WeaklyTideRunoff ####
WTvsWTR_pairwisefiltered <- read.csv("C:\\Users\\antre\\Box\\Bernal_lab\\Antrelle\\GulfKillifish\\WeakleyTide_vs_WeakleyTideRunoff_pairwisejoined_filtered.csv", header = TRUE, sep = ",", fill = TRUE)
colnames(WTvsWTR_pairwisefiltered)
  # [1] "X"              "EnsemblID"      "log2FoldChange" "padj"           "Gene_name"      "GO_ID
#Subset the data to only keep the columns we want: EnsemblID and log2FoldChange
WTvsWTR_EnsemblLogFC <- WTvsWTR_pairwisefiltered[c("EnsemblID", "log2FoldChange")]
WTvsWTR_Ensembl_LogFC <- WTvsWTR_EnsemblLogFC %>%
  distinct()
# Extract table as a .csv file
write.table(WTvsWTR_Ensembl_LogFC, "WTvsWTR_Ensembl_LogFC.csv", sep = ",", quote = F)


#### LongTideRunoff vs WeakleyTideRunoff ####
LTRvsWTR_pairwisefiltered <- read.csv("C:\\Users\\antre\\Box\\Bernal_lab\\Antrelle\\GulfKillifish\\LongTideRunoff_vs_WeakleyTideRunoff_pairwisejoined_filtered.csv", header = TRUE, sep = ",", fill = TRUE)
colnames(LTRvsWTR_pairwisefiltered)
  # [1] "X"              "EnsemblID"      "log2FoldChange" "padj"           "Gene_name"      "GO_ID
#Subset the data to only keep the columns we want: EnsemblID and log2FoldChange
LTRvsWTR_EnsemblLogFC <- LTRvsWTR_pairwisefiltered[c("EnsemblID", "log2FoldChange")]
LTRvsWTR_Ensembl_LogFC <- LTRvsWTR_EnsemblLogFC %>%
  distinct()
# Extract table as a .csv file
write.table(LTRvsWTR_Ensembl_LogFC, "LTRvsWTR_Ensembl_LogFC.csv", sep = ",", quote = F)


#### LongTideRunoff vs WeakleyTide ####
LTRvsWT_pairwisefiltered <- read.csv("C:\\Users\\antre\\Box\\Bernal_lab\\Antrelle\\GulfKillifish\\LongTideRunoff_vs_WeakleyTide_pairwisejoined_filtered.csv", header = TRUE, sep = ",", fill = TRUE)
colnames(LTRvsWT_pairwisefiltered)
# [1] "X"              "EnsemblID"      "log2FoldChange" "padj"           "Gene_name"      "GO_ID
#Subset the data to only keep the columns we want: EnsemblID and log2FoldChange
LTRvsWT_EnsemblLogFC <- LTRvsWT_pairwisefiltered[c("EnsemblID", "log2FoldChange")]
LTRvsWT_Ensembl_LogFC <- LTRvsWT_EnsemblLogFC %>%
  distinct()
# Extract table as a .csv file
write.table(LTRvsWT_Ensembl_LogFC, "LTRvsWT_Ensembl_LogFC.csv", sep = ",", quote = F)




######################## FOR ALL GENES #########################

# For each pairwise comparison, remove columns that are not EnsemblID and Log2FoldChange

#### LongTide vs WeakleyTide ####
LTvsWT_pairwisefiltered_ALLGENES <- read.csv("C:\\Users\\antre\\Box\\Bernal_lab\\Antrelle\\GulfKillifish\\LongTide_vs_WeakleyTide_ALLGENES_pairwisejoined_filtered.csv", header = TRUE, sep = ",", fill = TRUE)
colnames(LTvsWT_pairwisefiltered_ALLGENES)
#   [1] "EnsemblID"      "log2FoldChange" "padj"           "Gene_name"      "GO_ID"
#Subset the data to only keep the columns we want: EnsemblID and log2FoldChange
LTvsWT_EnsemblLogFC_ALLGENES <- LTvsWT_pairwisefiltered_ALLGENES[c("EnsemblID", "log2FoldChange")]
#Filter data to remove repeated lines
LTvsWT_Ensembl_LogFC_ALLGENES <- LTvsWT_EnsemblLogFC_ALLGENES %>%
  distinct()
# Extract table as a .csv file
write.table(LTvsWT_Ensembl_LogFC_ALLGENES, "LTvsWT_Ensembl_LogFC_ALLGENES.csv", sep = ",", quote = F) #renamed to LTvsWT_ELFC_AG.csv for GO_MWU_GFKanalysis

#### LongTide vs. LongTideRunoff ####
LTvsLTR_pairwisefiltered_ALLGENES <- read.csv("C:\\Users\\antre\\Box\\Bernal_lab\\Antrelle\\GulfKillifish\\LongTide_vs_LongTideRunoff_ALLGENES_pairwisejoined_filtered.csv", header = TRUE, sep = ",", fill = TRUE)
colnames(LTvsLTR_pairwisefiltered_ALLGENES)
#   [1] "EnsemblID"      "log2FoldChange" "padj"           "Gene_name"      "GO_ID"
#Subset the data to only keep the columns we want: EnsemblID and log2FoldChange
LTvsLTR_EnsemblLogFC_ALLGENES <- LTvsLTR_pairwisefiltered_ALLGENES[c("EnsemblID", "log2FoldChange")]
#Filter data to remove repeated lines
LTvsLTR_Ensembl_LogFC_ALLGENES <- LTvsLTR_EnsemblLogFC_ALLGENES %>%
  distinct()
# Extract table as a .csv file
write.table(LTvsLTR_Ensembl_LogFC_ALLGENES, "LTvsLTR_Ensembl_LogFC_ALLGENES.csv", sep = ",", quote = F) #renamed to LTvsLTR_ELFC_AG.csv for GO_MWU_GFKanalysis


#### LongTide vs. WeakleyTideRunoff ####
LTvsWTR_pairwisefiltered_ALLGENES <- read.csv("C:\\Users\\antre\\Box\\Bernal_lab\\Antrelle\\GulfKillifish\\LongTide_vs_WeakleyTideRunoff_ALLGENES_pairwisejoined_filtered.csv", header = TRUE, sep = ",", fill = TRUE)
colnames(LTvsWTR_pairwisefiltered_ALLGENES)
#   [1] "EnsemblID"      "log2FoldChange" "padj"           "Gene_name"      "GO_ID"
#Subset the data to only keep the columns we want: EnsemblID and log2FoldChange
LTvsWTR_EnsemblLogFC_ALLGENES <- LTvsWTR_pairwisefiltered_ALLGENES[c("EnsemblID", "log2FoldChange")]
#Filter data to remove repeated lines
LTvsWTR_Ensembl_LogFC_ALLGENES <- LTvsWTR_EnsemblLogFC_ALLGENES %>%
  distinct()
# Extract table as a .csv file
write.table(LTvsWTR_Ensembl_LogFC_ALLGENES, "LTvsWTR_Ensembl_LogFC_ALLGENES.csv", sep = ",", quote = F) #renamed to LTvsWTR_ELFC_AG.csv for GO_MWU_GFKanalysis



#### WeakleyTide vs. WeakleyTideRunoff ####
WTvsWTR_pairwisefiltered_ALLGENES <- read.csv("C:\\Users\\antre\\Box\\Bernal_lab\\Antrelle\\GulfKillifish\\WeakleyTide_vs_WeakleyTideRunoff_ALLGENES_pairwisejoined_filtered.csv", header = TRUE, sep = ",", fill = TRUE)
colnames(WTvsWTR_pairwisefiltered_ALLGENES)
#   [1] "EnsemblID"      "log2FoldChange" "padj"           "Gene_name"      "GO_ID"
#Subset the data to only keep the columns we want: EnsemblID and log2FoldChange
WTvsWTR_EnsemblLogFC_ALLGENES <- WTvsWTR_pairwisefiltered_ALLGENES[c("EnsemblID", "log2FoldChange")]
#Filter data to remove repeated lines
WTvsWTR_Ensembl_LogFC_ALLGENES <- WTvsWTR_EnsemblLogFC_ALLGENES %>%
  distinct()
# Extract table as a .csv file
write.table(WTvsWTR_Ensembl_LogFC_ALLGENES, "WTvsWTR_Ensembl_LogFC_ALLGENES.csv", sep = ",", quote = F) #renamed to WTvsWTR_ELFC_AG.csv for GO_MWU_GFKanalysis



#### LongTideRunoff vs. WeakleyTideRunoff ####
LTRvsWTR_pairwisefiltered_ALLGENES <- read.csv("C:\\Users\\antre\\Box\\Bernal_lab\\Antrelle\\GulfKillifish\\LongTideRunoff_vs_WeakleyTideRunoff_ALLGENES_pairwisejoined_filtered.csv", header = TRUE, sep = ",", fill = TRUE)
colnames(LTRvsWTR_pairwisefiltered_ALLGENES)
#   [1] "EnsemblID"      "log2FoldChange" "padj"           "Gene_name"      "GO_ID"
#Subset the data to only keep the columns we want: EnsemblID and log2FoldChange
LTRvsWTR_EnsemblLogFC_ALLGENES <- LTRvsWTR_pairwisefiltered_ALLGENES[c("EnsemblID", "log2FoldChange")]
#Filter data to remove repeated lines
LTRvsWTR_Ensembl_LogFC_ALLGENES <- LTRvsWTR_EnsemblLogFC_ALLGENES %>%
  distinct()
# Extract table as a .csv file
write.table(LTRvsWTR_Ensembl_LogFC_ALLGENES, "LTRvsWTR_Ensembl_LogFC_ALLGENES.csv", sep = ",", quote = F) #renamed to LTRvsWTR_ELFC_AG.csv for GO_MWU_GFKanalysis



#### LongTideRunoff vs. WeakleyTide ####
LTRvsWT_pairwisefiltered_ALLGENES <- read.csv("C:\\Users\\antre\\Box\\Bernal_lab\\Antrelle\\GulfKillifish\\LongTideRunoff_vs_WeakleyTide_ALLGENES_pairwisejoined_filtered.csv", header = TRUE, sep = ",", fill = TRUE)
colnames(LTRvsWT_pairwisefiltered_ALLGENES)
#   [1] "EnsemblID"      "log2FoldChange" "padj"           "Gene_name"      "GO_ID"
#Subset the data to only keep the columns we want: EnsemblID and log2FoldChange
LTRvsWT_EnsemblLogFC_ALLGENES <- LTRvsWT_pairwisefiltered_ALLGENES[c("EnsemblID", "log2FoldChange")]
#Filter data to remove repeated lines
LTRvsWT_Ensembl_LogFC_ALLGENES <- LTRvsWT_EnsemblLogFC_ALLGENES %>%
  distinct()
# Extract table as a .csv file
write.table(LTRvsWT_Ensembl_LogFC_ALLGENES, "LTRvsWT_Ensembl_LogFC_ALLGENES.csv", sep = ",", quote = F) #renamed to WTvsLTR_ELFC_AG.csv for GO_MWU_GFKanalysis


################################################################################


# Next step is to create a table that has all the GO_IDs on one line separated by ";" that corresponds to the same gene: EnsemblID and GO_ID

# Go back to reference data table that was created: MummichogDataInfoforGOterm
MummichogGOreferencetable <- read.csv("C:\\Users\\antre\\Box\\Bernal_lab\\Antrelle\\GulfKillifish\\MummichogDataInfoforGOterm.csv", header = TRUE, sep = ",", fill = TRUE)
# Combine all GO IDs for the same EnsemblID
MummichogGOref_combinedGOIDs <- MummichogGOreferencetable %>%
  group_by(EnsemblID) %>% #grouping the data by the EnsemblID column
  summarize(GO_ID = paste(GO_ID, collapse = ";")) #concentrates all GO IDs for a single gene, separating them by semicolons
print(MummichogGOref_combinedGOIDs) #check to make sure the GO IDs for a single gene has been combined
# Extract table as a .csv file
write.csv(MummichogGOref_combinedGOIDs, "Mummichog_aggregated_GO_terms.csv", row.names = FALSE)
# Subset the table that genes from the reference that is blank or has commas in it are deleted (essentially those that do not have GO_IDs)
MummichogGOref_combinedGOIDs_removed <- subset(MummichogGOref_combinedGOIDs, !(GO_ID == "" | GO_ID == "," | GO_ID == ",," | is.na(GO_ID)))
print(MummichogGOref_combinedGOIDs_removed) #check to make sure rows with blank GO IDs are removed from the reference data
# Extract table as a .csv file
write.csv(MummichogGOref_combinedGOIDs_removed, "Mummichog_aggregated_GO_terms_filtered.csv", row.names = FALSE)



# Proceed to GO_MWU_GFKanalysis.R


