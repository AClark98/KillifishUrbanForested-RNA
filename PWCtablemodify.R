################# Pairwise Comparison Table Modify for Gulf Killifish #######################

# Antrelle D. Clark, Auburn University Department of Biological Sciences, Bernal Lab, Fall 2024; adc0064@auburn.edu
# R-script version 2024.09.0+375


#################################################################

#Install packages needed
library("BiocManager")
BiocManager::install("clusterProfiler") #the package that is doing the actual GO analysis
BiocManager::install("AnnotationDbi") #the package running in the background
install.packages("GO.db")
install.packages("GO")

#Load the library
library(clusterProfiler)
library(dplyr)
library(GO.db)

setwd("/Users/antre/Box/Bernal_lab/Antrelle/GulfKillifish")


# The GO terms were retrieved from http://useast.ensembl.org/biomart/martview/ via the following steps:
# In the drop down menu [-CHOOSE DATABASE -]  select [Ensembl Genes 113]
# In the [-CHOOSE DATASET MENU -] select the organism you want to use as a reference - my reference for this dataset is: Mummichog genes (Fundulus_heteroclitus-3.0.2)
# Under [Attributes] select [Features], under [Gene] check the [Gene Stable ID], [Gene description] and [Gene name] box; the [Transcript stable ID] might be checked and that is fine.
# Under [External] check the [GO term accession] abd [GO term name] boxes
# Click [Results], from the file type drop down menu select [CSV]
# Click [Export all results to COmpressed web file (notify by email)]
# Enter your email and hit [Go]
# You should get a link to download your database sequencese as a .txt.gz file in fasta format

#OUTPUT:
  #Your results are being compiled in the background.
  #Your reference is martquery_1108163932_584.txt.gz.
  #An email will be sent to you when they are ready.
  #^ this file will include the genes and what they code for

#Transfer to supercomputer:
# PS C:\Users\antre> scp -r "C:\Users\antre\Box\Bernal_lab\Antrelle\GulfKillifish\martquery_1108163932_584.txt" aubadc002@asax.asc.edu:"\home\aubadc002\killifish\gfkGOtermtest"
# Usually one would have to blast multiple fish species to get an accurate GO profile if their species isn't readily available; however, Mummichog is sister species and also the reference genome used for filtering reads so step skipped

# Download notepad ++ (for Windows) or 
# Open pairwise files and remove the word "gene" from in front of all genes, For example: change gene:ENSFHEG00000008393 to ENSFHEG00000008393 and save the file as the modified version of the original
# Give the gene column the name "EnsemblID" because the columns have to match in order to merge
# Also modify the martquery_1108163932_584 name and change it to: MummichogbiomaRtdata and then also change the "gene stable ID" column name to "Ensembl ID"


# Make two tables using geneIDs original table by removing columns: one that will contain EnsemblID and GO_ID and one that will contain EnsemblID and Gene_name
    # EnsemblID and GO_ID

#library(dplyr)

#### the following code is what you would use for manually deleting columns (this is the method used for this analysis) ####
#download the Mummichog data (will be downloaded as a .txt file)
# open excel and go to the "data" tab, and then "get & transform data: from text/csv" and import the file
# In the Mummichog data, change the following names:
# Gene Stable ID ---> EnsemblID
# Gene Name ---> Gene_name
# GO term accession ---> GO_ID
# then:
# manually delete the other columns because we don't want those
MummichogbiomaRtdata <- read.csv("C:\\Users\\antre\\Box\\Bernal_lab\\Antrelle\\GulfKillifish\\MummichogDataInfoforGOterm.csv", header = TRUE, sep = ",", fill = TRUE)
colnames(MummichogbiomaRtdata)
  # [1] "EnsemblID" "Gene_name" "GO_ID"
# Subset the data frame to select only the columns you want to keep: EnsemblID and GO_ID (this could also be done manually by saving the file under another name and deleting the row you don't want and saving)
EnsemblID_GO_ID_df <- MummichogbiomaRtdata[c("EnsemblID", "GO_ID")]
colnames(EnsemblID_GO_ID_df)
  # [1] "EnsemblID" "GO_ID"
# Extract the table as a .csv
write.table(EnsemblID_GO_ID_df, "EnsemblID_GO_ID_df.csv", sep = ",", quote = F)
# Subset the data frame to select only the columns you want to keep: EnsemblID and Gene_name
EnsemblID_Gene_name_df <- MummichogbiomaRtdata[c("EnsemblID", "Gene_name")]
colnames(EnsemblID_Gene_name_df)
  # [1] "EnsemblID" "Gene_name"
# Extract the table as a .csv
write.table(EnsemblID_Gene_name_df, "EnsemblID_Gene_name_df.csv", sep = ",", quote = F)
# save as "C:/Users/antre/Box/Bernal_lab/Antrelle/GulfKillifish/EnsemblID_GO_ID_modified.csv" and "C:/Users/antre/Box/Bernal_lab/Antrelle/GulfKillifish/EnsemblID_Gene_name_modified.csv"


######## Merge pairwise tables containing significant genes, EnsemblID, GO_ID, and Gene_name ########

# Use notepad++ to manually add a column before EnsemblID to shift all of the data over (and just number them as shown as the line of code; ex: code line 1 so just place a 1 there to fill it because you cannot have the same row names) so that R can read EnsemblID as a column

######################## FOR SIGNIFICANT GENES #########################

#### LongTide vs LongTideRunoff ####
LongTide_vs_LongTiderunoff_modified<-read.table("C:/Users/antre/Box/Bernal_lab/Antrelle/GulfKillifish/longtide_vs_longtiderunoff_pairwise_gkf_modified.txt", header = TRUE, row.names = 1)
colnames(LongTide_vs_LongTiderunoff_modified)
  # [1] "EnsemblID"      "baseMean"       "log2FoldChange" "lfcSE"          "stat"           "pvalue"        
  # [7] "padj"
colnames(EnsemblID_GO_ID_df)
  # [1] "EnsemblID" "GO_ID"
# First, merge the first two data frames: LongTide_vs_LongTiderunoff_modified + EnsemblID_GO_ID_df
LTvsLTR_merged_df_1 <- merge(LongTide_vs_LongTiderunoff_modified, EnsemblID_GO_ID_df, by = c("EnsemblID"))
# Then, merge the result with the third data frame: LTvsLTR_merged_df_1 + EnsemblID_Gene_name_df; !!!! check colnames first to make sure there is a column with the same name in each: EnsemblID! !!!!
colnames(LTvsLTR_merged_df_1)
  # [1] "EnsemblID"      "baseMean"       "log2FoldChange" "lfcSE"          "stat"           "pvalue"        
  # [7] "padj"           "GO_ID"  
colnames(EnsemblID_Gene_name_df)
  # [1] "EnsemblID" "Gene_name"
LongTide_vs_LongTiderunoff_modified_merged <- merge(LTvsLTR_merged_df_1, EnsemblID_Gene_name_df, by = c("EnsemblID"))
colnames(LongTide_vs_LongTiderunoff_modified_merged)
  # [1] "EnsemblID"      "baseMean"       "log2FoldChange" "lfcSE"          "stat"           "pvalue"        
  # [7] "padj"           "GO_ID"          "Gene_name"
# Subset the table to only include "EnsemblID" , "log2FoldChange", "padj", "Gene_name", "GO_ID"
LongTide_vs_LongTiderunoff_modified_merged_subset <- subset(LongTide_vs_LongTiderunoff_modified_merged, select = c(EnsemblID, log2FoldChange, padj, Gene_name, GO_ID)) #you are specifying the columns you want to keep
colnames(LongTide_vs_LongTiderunoff_modified_merged_subset)
  # [1] "EnsemblID"      "log2FoldChange" "padj"           "Gene_name"      "GO_ID"
# Extract the table as a .csv
write.table(LongTide_vs_LongTiderunoff_modified_merged_subset, "LongTide_vs_LongTideRunoff_pairwisejoined.csv", sep = ",", quote = F)
# Filter the data to make sure there are no duplicate rows (rows may contain the same gene, which is fine, but the GO_ID must be different)
LongTide_vs_LongTideRunoff_modified_merged_subset_filtered <- LongTide_vs_LongTiderunoff_modified_merged_subset %>%
  distinct()
# Extract the filtered table as a .csv
write.table(LongTide_vs_LongTideRunoff_modified_merged_subset_filtered, "LongTide_vs_LongTideRunoff_pairwisejoined_filtered.csv", sep = ",", quote = F)
# originally had 343 rows; now it has 41


#### LongTide vs WeakleyTide ####

# Use notepad++ to manually add a column before EnsemblID to shift all of the data over
LongTide_vs_WeakleyTide_modified<-read.table("C:/Users/antre/Box/Bernal_lab/Antrelle/GulfKillifish/longtide_vs_weakleytide_pairwise_gfk_modified.txt", header = TRUE, row.names = 1)
colnames(LongTide_vs_WeakleyTide_modified)
  # [1] "EnsemblID"      "baseMean"       "log2FoldChange" "lfcSE"          "stat"           "pvalue"        
  # [7] "padj"
colnames(EnsemblID_GO_ID_df)
  # [1] "EnsemblID" "GO_ID"
# First, merge the first two data frames: LongTide_vs_WeakleyTide_modified + EnsemblID_GO_ID_df
LTvsWT_merged_df_1 <- merge(LongTide_vs_WeakleyTide_modified, EnsemblID_GO_ID_df, by = c("EnsemblID"))
# Then, merge the result with the third data frame:
head(LTvsWT_merged_df_1)
# Then, merge the result with the third data frame:
colnames(LTvsWT_merged_df_1)
  # [1] "EnsemblID"      "baseMean"       "log2FoldChange" "lfcSE"          "stat"           "pvalue"        
  # [7] "padj"           "GO_ID" 
colnames(EnsemblID_Gene_name_df)
  # [1] "EnsemblID" "Gene_name"
LongTide_vs_WeakleyTide_modified_merged <- merge(LTvsWT_merged_df_1, EnsemblID_Gene_name_df, by = c("EnsemblID"))
colnames(LongTide_vs_WeakleyTide_modified_merged)
  # [1] "EnsemblID"      "baseMean"       "log2FoldChange" "lfcSE"          "stat"           "pvalue"        
  # [7] "padj"           "GO_ID"          "Gene_name" 
head(LongTide_vs_WeakleyTide_modified_merged) #checking to ensure data is present
# Subset the table to only include "EnsemblID" , "log2FoldChange", "padj", "Gene_name", "GO_ID"
LongTide_vs_WeakleyTide_modified_merged_subset <- subset(LongTide_vs_WeakleyTide_modified_merged, select = c(EnsemblID, log2FoldChange, padj, Gene_name, GO_ID)) #you are specifying the columns you want to keep
colnames(LongTide_vs_WeakleyTide_modified_merged_subset)
  # [1] "EnsemblID"      "log2FoldChange" "padj"           "Gene_name"      "GO_ID" 
# Extract the table as a .csv
write.table(LongTide_vs_WeakleyTide_modified_merged_subset, "LongTide_vs_WeakleyTide_pairwisejoined.csv", sep = ",", quote = F)
# Filter the data to make sure there are no duplicate rows (rows may contain the same gene, which is fine, but the GO_ID must be different)
LongTide_vs_WeakleyTide_modified_merged_subset_filtered <- LongTide_vs_WeakleyTide_modified_merged_subset %>%
  distinct()
# Extract the filtered table as a .csv
write.table(LongTide_vs_WeakleyTide_modified_merged_subset_filtered, "LongTide_vs_WeakleyTide_pairwisejoined_filtered.csv", sep = ",", quote = F)
# originally had 251 rows; now it has 30


#### LongTide vs WeakleyTideRunoff ####

# Use notepad++ to manually add a column before EnsemblID to shift all of the data over
LongTide_vs_WeakleyTideRunoff_modified<-read.table("C:/Users/antre/Box/Bernal_lab/Antrelle/GulfKillifish/longtide_vs_weakleytiderunoff_pairwise_gfk_modified.txt", header = TRUE, row.names = 1)
colnames(LongTide_vs_WeakleyTideRunoff_modified)
  # [1] "EnsemblID"      "baseMean"       "log2FoldChange" "lfcSE"          "stat"           "pvalue"        
  # [7] "padj"
# First, merge the first two data frames: LongTide_vs_WeakleyTideRunoff_modified + EnsemblID_GO_ID_df
LTvsWTR_merged_df_1 <- merge(LongTide_vs_WeakleyTideRunoff_modified, EnsemblID_GO_ID_df, by = c("EnsemblID"))
# Then, merge the result with the third data frame:
colnames(LTvsWTR_merged_df_1)
  # [1] "EnsemblID"      "baseMean"       "log2FoldChange" "lfcSE"          "stat"           "pvalue"        
  # [7] "padj"           "GO_ID" 
colnames(EnsemblID_Gene_name_df)
  # [1] "EnsemblID" "Gene_name"
LongTide_vs_WeakleyTideRunoff_modified_merged <- merge(LTvsWTR_merged_df_1, EnsemblID_Gene_name_df, by = c("EnsemblID"))
colnames(LongTide_vs_WeakleyTideRunoff_modified_merged)
  # [1] "EnsemblID"      "baseMean"       "log2FoldChange" "lfcSE"          "stat"           "pvalue"        
  # [7] "padj"           "GO_ID"          "Gene_name"
# Subset the table to only include "EnsemblID" , "log2FoldChange", "padj", "Gene_name", "GO_ID"
LongTide_vs_WeakleyTideRunoff_modified_merged_subset <- subset(LongTide_vs_WeakleyTideRunoff_modified_merged, select = c(EnsemblID, log2FoldChange, padj, Gene_name, GO_ID)) #you are specifying the columns you want to keep
colnames(LongTide_vs_WeakleyTideRunoff_modified_merged_subset)
  # [1] "EnsemblID"      "log2FoldChange" "padj"           "Gene_name"      "GO_ID" 
# Subset the table to only include "EnsemblID" , "log2FoldChange", "padj", "Gene_name", "GO_ID"
LongTide_vs_WeakleyTideRunoff_modified_merged_subset <- subset(LongTide_vs_WeakleyTideRunoff_modified_merged, select = c(EnsemblID, log2FoldChange, padj, Gene_name, GO_ID)) #you are specifying the columns you want to keep
colnames(LongTide_vs_WeakleyTideRunoff_modified_merged_subset)
  # [1] "EnsemblID"      "log2FoldChange" "padj"           "Gene_name"      "GO_ID" 
# Extract the table as a .csv
write.table(LongTide_vs_WeakleyTideRunoff_modified_merged_subset, "LongTide_vs_WeakleyTideRunoff_pairwisejoined.csv", sep = ",", quote = F)
# Filter the data to make sure there are no duplicate rows (rows may contain the same gene, which is fine, but the GO_ID must be different)
LongTide_vs_WeakleyTideRunoff_modified_merged_subset_filtered <- LongTide_vs_WeakleyTideRunoff_modified_merged_subset %>%
  distinct()
# Extract the filtered table as a .csv
write.table(LongTide_vs_WeakleyTideRunoff_modified_merged_subset_filtered, "LongTide_vs_WeakleyTideRunoff_pairwisejoined_filtered.csv", sep = ",", quote = F)
# originally had 40570 rows; now it has 1323 


#### Weakley Tide vs WeaklyTideRunoff ####

# Use notepad++ to manually add a column before EnsemblID to shift all of the data over
WeakleyTide_vs_WeakleyTideRunoff_modified<-read.table("C:/Users/antre/Box/Bernal_lab/Antrelle/GulfKillifish/weakleytide_vs_weakleytiderunoff_pairwise_gfk_modified.txt", header = TRUE, row.names = 1)
colnames(WeakleyTide_vs_WeakleyTideRunoff_modified)
  # [1] "EnsemblID"      "baseMean"       "log2FoldChange" "lfcSE"          "stat"           "pvalue"        
  # [7] "padj"
colnames(EnsemblID_GO_ID_df)
  # [1] "EnsemblID" "GO_ID" 
# First, merge the first two data frames: WeakleyTide_vs_WeakleyTideRunoff_modified + EnsemblID_GO_ID_df
WTvsWTR_merged_df_1 <- merge(WeakleyTide_vs_WeakleyTideRunoff_modified, EnsemblID_GO_ID_df, by = c("EnsemblID"))
# Then, merge the result with the third data frame:
colnames(WTvsWTR_merged_df_1)
  # [1] "EnsemblID"      "baseMean"       "log2FoldChange" "lfcSE"          "stat"           "pvalue"        
  # [7] "padj"           "GO_ID" 
colnames(EnsemblID_Gene_name_df)
  # [1] "EnsemblID" "Gene_name"
WeakleyTide_vs_WeakleyTideRunoff_modified_merged <- merge(WTvsWTR_merged_df_1, EnsemblID_Gene_name_df, by = c("EnsemblID"))
colnames(WeakleyTide_vs_WeakleyTideRunoff_modified_merged)
  # [1] "EnsemblID"      "baseMean"       "log2FoldChange" "lfcSE"          "stat"           "pvalue"        
  # [7] "padj"           "GO_ID"          "Gene_name"
# Subset the table to only include "EnsemblID" , "log2FoldChange", "padj", "Gene_name", "GO_ID"
WeakleyTide_vs_WeakleyTideRunoff_modified_merged_subset <- subset(WeakleyTide_vs_WeakleyTideRunoff_modified_merged, select = c(EnsemblID, log2FoldChange, padj, Gene_name, GO_ID)) #you are specifying the columns you want to keep
colnames(WeakleyTide_vs_WeakleyTideRunoff_modified_merged_subset)
  # [1] "EnsemblID"      "log2FoldChange" "padj"           "Gene_name"      "GO_ID"
# Extract the table as a .csv
write.table(WeakleyTide_vs_WeakleyTideRunoff_modified_merged_subset, "WeakleyTide_vs_WeakleyTideRunoff_pairwisejoined.csv", sep = ",", quote = F)
# Filter the data to make sure there are no duplicate rows (rows may contain the same gene, which is fine, but the GO_ID must be different)
WeakleyTide_vs_WeakleyTideRunoff_modified_merged_subset_filtered <- WeakleyTide_vs_WeakleyTideRunoff_modified_merged_subset %>%
  distinct()
# Extract the filtered table as a .csv
write.table(WeakleyTide_vs_WeakleyTideRunoff_modified_merged_subset_filtered, "WeakleyTide_vs_WeakleyTideRunoff_pairwisejoined_filtered.csv", sep = ",", quote = F)
# originally had 30771 rows; now it has 941 


#### LongTideRunoff vs WeakleyTideRunoff ####

# Use notepad++ to manually add a column before EnsemblID to shift all of the data over
LongTideRunoff_vs_WeakleyTideRunoff_modified<-read.table("C:/Users/antre/Box/Bernal_lab/Antrelle/GulfKillifish/longtiderunoff_vs_weakleytiderunoff_pairwise_gfk_modified.txt", header = TRUE, row.names = 1)
colnames(LongTideRunoff_vs_WeakleyTideRunoff_modified)
  # [1] "EnsemblID"      "baseMean"       "log2FoldChange" "lfcSE"          "stat"           "pvalue"        
  # [7] "padj"
colnames(EnsemblID_GO_ID_df)
  # [1] "EnsemblID" "GO_ID"
# First, merge the first two data frames: LongTideRunoff_vs_WeakleyTideRunoff_modified + EnsemblID_GO_ID_df
LTRvsWTR_merged_df_1 <- merge(LongTideRunoff_vs_WeakleyTideRunoff_modified, EnsemblID_GO_ID_df, by = c("EnsemblID"))
colnames(LTRvsWTR_merged_df_1)
  # [1] "EnsemblID"      "baseMean"       "log2FoldChange" "lfcSE"          "stat"           "pvalue"        
  # [7] "padj"           "GO_ID"
# Then, merge the result with the third data frame:
colnames(EnsemblID_Gene_name_df)
  # [1] "EnsemblID" "Gene_name"
LongTideRunoff_vs_WeakleyTideRunoff_modified_merged <- merge(LTRvsWTR_merged_df_1, EnsemblID_Gene_name_df, by = c("EnsemblID"))
colnames(LongTideRunoff_vs_WeakleyTideRunoff_modified_merged)
  # [1] "EnsemblID"      "baseMean"       "log2FoldChange" "lfcSE"          "stat"           "pvalue"        
  # [7] "padj"           "GO_ID"          "Gene_name"
# Subset the table to only include "EnsemblID" , "log2FoldChange", "padj", "Gene_name", "GO_ID"
LongTideRunoff_vs_WeakleyTideRunoff_modified_merged_subset <- subset(LongTideRunoff_vs_WeakleyTideRunoff_modified_merged, select = c(EnsemblID, log2FoldChange, padj, Gene_name, GO_ID)) #you are specifying the columns you want to keep
colnames(LongTideRunoff_vs_WeakleyTideRunoff_modified_merged_subset)
  # [1] "EnsemblID"      "log2FoldChange" "padj"           "Gene_name"      "GO_ID"  
# Extract the table as a .csv
write.table(LongTideRunoff_vs_WeakleyTideRunoff_modified_merged_subset, "LongTideRunoff_vs_WeakleyTideRunoff_pairwisejoined.csv", sep = ",", quote = F)
# Filter the data to make sure there are no duplicate rows (rows may contain the same gene, which is fine, but the GO_ID must be different)
LongTideRunoff_vs_WeakleyTideRunoff_modified_merged_subset_filtered <- LongTideRunoff_vs_WeakleyTideRunoff_modified_merged_subset %>%
  distinct()
# Extract the filtered table as a .csv
write.table(LongTideRunoff_vs_WeakleyTideRunoff_modified_merged_subset_filtered, "LongTideRunoff_vs_WeakleyTideRunoff_pairwisejoined_filtered.csv", sep = ",", quote = F)
# originally had 40750 rows; now it has 1323


#### LongTideRunoff vs WeakleyTide ####

# Use notepad++ to manually add a column before EnsemblID to shift all of the data over
LongTideRunoff_vs_WeakleyTide_modified<-read.table("C:/Users/antre/Box/Bernal_lab/Antrelle/GulfKillifish/longtiderunoff_vs_weakleytide_pairwise_gfk_modified.txt", header = TRUE, row.names = 1)
colnames(LongTideRunoff_vs_WeakleyTide_modified)
  # [1] "EnsemblID"      "baseMean"       "log2FoldChange" "lfcSE"          "stat"           "pvalue"        
  # [7] "padj"
colnames(EnsemblID_GO_ID_df)
  # [1] "EnsemblID" "GO_ID"
# First, merge the first two data frames: LongTideRunoff_vs_WeakleyTideRunoff_modified + EnsemblID_GO_ID_df
LTRvsWT_merged_df_1 <- merge(LongTideRunoff_vs_WeakleyTide_modified, EnsemblID_GO_ID_df, by = c("EnsemblID"))
colnames(LTRvsWT_merged_df_1)
  # [1] "EnsemblID"      "baseMean"       "log2FoldChange" "lfcSE"          "stat"           "pvalue"        
  # [7] "padj"           "GO_ID"
# Then, merge the result with the third data frame:
colnames(EnsemblID_Gene_name_df)
  # [1] "EnsemblID" "Gene_name"
LongTideRunoff_vs_WeakleyTide_modified_merged <- merge(LTRvsWT_merged_df_1, EnsemblID_Gene_name_df, by = c("EnsemblID"))
colnames(LongTideRunoff_vs_WeakleyTide_modified_merged)
  # [1] "EnsemblID"      "baseMean"       "log2FoldChange" "lfcSE"          "stat"           "pvalue"        
  # [7] "padj"           "GO_ID"          "Gene_name" 
# Subset the table to only include "EnsemblID" , "log2FoldChange", "padj", "Gene_name", "GO_ID"
LongTideRunoff_vs_WeakleyTide_modified_merged_subset <- subset(LongTideRunoff_vs_WeakleyTide_modified_merged, select = c(EnsemblID, log2FoldChange, padj, Gene_name, GO_ID)) #you are specifying the columns you want to keep
colnames(LongTideRunoff_vs_WeakleyTide_modified_merged_subset)
  # [1] "EnsemblID"      "log2FoldChange" "padj"           "Gene_name"      "GO_ID" 
# Extract the table as a .csv
write.table(LongTideRunoff_vs_WeakleyTide_modified_merged_subset, "LongTideRunoff_vs_WeakleyTide_pairwisejoined.csv", sep = ",", quote = F)
# Filter the data to make sure there are no duplicate rows (rows may contain the same gene, which is fine, but the GO_ID must be different)
LongTideRunoff_vs_WeakleyTide_modified_merged_subset_filtered <- LongTideRunoff_vs_WeakleyTide_modified_merged_subset %>%
  distinct()
# Extract the filtered table as a .csv
write.table(LongTideRunoff_vs_WeakleyTide_modified_merged_subset_filtered, "LongTideRunoff_vs_WeakleyTide_pairwisejoined_filtered.csv", sep = ",", quote = F)
# originally had 1286 rows; now it has 113




######################## FOR ALL GENES #########################

# Use notepad++ to manually replace the column name "basemean" with EnsemblID and \t to shift all of the data over #see Adding_columns_in_notepadplusplus script to easily perform this


#### LongTide vs WeakleyTide ####

# Use notepad++ to manually add a column before EnsemblID to shift all of the data over
LongTide_vs_WeakleyTide_modified_ALLGENES<-read.table("C:\\Users\\antre\\Box\\Bernal_lab\\Antrelle\\GulfKillifish\\longtide_vs_weakleytide_pairwise_gkf_ALLGENES.txt", header = TRUE, row.names = 1)
colnames(LongTide_vs_WeakleyTide_modified_ALLGENES)
  #[1] "X1"             "EnsemblID"      "baseMean"       "log2FoldChange" "lfcSE"          "stat"          
  #[7] "pvalue"         "padj" 
colnames(EnsemblID_GO_ID_df)
  # [1] "EnsemblID" "GO_ID"
# First, merge the first two data frames: LongTide_vs_WeakleyTide_modified_ALLGENES + EnsemblID_GO_ID_df
LTvsWT_merged_df_AG <- merge(LongTide_vs_WeakleyTide_modified_ALLGENES, EnsemblID_GO_ID_df, by = c("EnsemblID"))
# Then, merge the result with the third data frame:
LongTide_vs_WeakleyTide_modified_merged_ALLGENES <- merge(LTvsWT_merged_df_AG, EnsemblID_Gene_name_df, by = c("EnsemblID"))
colnames(LongTide_vs_WeakleyTide_modified_merged_ALLGENES)
  #[1] "EnsemblID"      "X1"             "baseMean"       "log2FoldChange" "lfcSE"          "stat"          
  #[7] "pvalue"         "padj"           "GO_ID"          "Gene_name"     
# Subset the table to only include "EnsemblID" , "log2FoldChange", "padj", "Gene_name", "GO_ID"
LongTide_vs_WeakleyTide_modified_merged_ALLGENES_subset <- subset(LongTide_vs_WeakleyTide_modified_merged_ALLGENES, select = c(EnsemblID, log2FoldChange, padj, Gene_name, GO_ID)) #you are specifying the columns you want to keep
colnames(LongTide_vs_WeakleyTide_modified_merged_ALLGENES_subset)
# [1] "EnsemblID"      "log2FoldChange" "padj"           "Gene_name"      "GO_ID" 
# Extract the table as a .csv
write.table(LongTide_vs_WeakleyTide_modified_merged_ALLGENES_subset, "LongTide_vs_WeakleyTide_ALLGENES_pairwisejoined.csv", sep = ",", quote = F)
# Filter the data to make sure there are no duplicate rows (rows may contain the same gene, which is fine, but the GO_ID must be different)
LongTide_vs_WeakleyTide_modified_merged_ALLGENES_subset_filtered <- LongTide_vs_WeakleyTide_modified_merged_ALLGENES_subset %>%
  distinct()
# Extract the filtered table as a .csv
write.table(LongTide_vs_WeakleyTide_modified_merged_ALLGENES_subset_filtered, "LongTide_vs_WeakleyTide_ALLGENES_pairwisejoined_filtered.csv", sep = ",", quote = F)
# originally had 1048575 rows; now it has 98539


#### LongTide vs. LongTideRunoff ####

# Use notepad++ to manually add a column before basemean, named "EnsemblID" to shift all of the data over and save the original .tsv as a .txt file
LongTide_vs_LongTideRunoff_modified_ALLGENES<-read.table("C:\\Users\\antre\\Box\\Bernal_lab\\Antrelle\\GulfKillifish\\longtide_vs_longtiderunoff_pairwise_gkf_ALLGENES.txt", header = TRUE, row.names = 1)
colnames(LongTide_vs_LongTideRunoff_modified_ALLGENES)
# First, merge the first two data frames: LongTide_vs_LongTideRunoff_modified_ALLGENES + EnsemblID_GO_ID_df
LTvsLTR_merged_df_AG <- merge(LongTide_vs_LongTideRunoff_modified_ALLGENES, EnsemblID_GO_ID_df, by = c("EnsemblID"))
# Then, merge the result with the third data frame:
LongTide_vs_LongTideRunoff_modified_merged_ALLGENES <- merge(LTvsLTR_merged_df_AG, EnsemblID_Gene_name_df, by = c("EnsemblID"))
colnames(LongTide_vs_LongTideRunoff_modified_merged_ALLGENES)
# Subset the table to only include "EnsemblID" , "log2FoldChange", "padj", "Gene_name", "GO_ID"
LongTide_vs_LongTideRunoff_modified_merged_ALLGENES_subset <- subset(LongTide_vs_LongTideRunoff_modified_merged_ALLGENES, select = c(EnsemblID, log2FoldChange, padj, Gene_name, GO_ID)) #you are specifying the columns you want to keep
colnames(LongTide_vs_LongTideRunoff_modified_merged_ALLGENES_subset)
# Extract the table as a .csv
write.table(LongTide_vs_LongTideRunoff_modified_merged_ALLGENES_subset, "LongTide_vs_LongTideRunoff_ALLGENES_pairwisejoined.csv", sep = ",", quote = F)
# Filter the data to make sure there are no duplicate rows (rows may contain the same gene, which is fine, but the GO_ID must be different)
LongTide_vs_LongTideRunoff_modified_merged_ALLGENES_subset_filtered <- LongTide_vs_LongTideRunoff_modified_merged_ALLGENES_subset %>%
  distinct()
# Extract the filtered table as a .csv
write.table(LongTide_vs_LongTideRunoff_modified_merged_ALLGENES_subset_filtered, "LongTide_vs_LongTideRunoff_ALLGENES_pairwisejoined_filtered.csv", sep = ",", quote = F)
# originally had 2182834 rows; now it has 98539




#### LongTide vs WeakleyTideRunoff ####
# Use notepad++ to manually add a column before EnsemblID to shift all of the data over and save the original .tsv as a .txt file
LongTide_vs_WeakleyTideRunoff_modified_ALLGENES<-read.table("C:\\Users\\antre\\Box\\Bernal_lab\\Antrelle\\GulfKillifish\\longtide_vs_weakleytiderunoff_pairwise_gkf_ALLGENES.txt", header = TRUE, row.names = 1)
colnames(LongTide_vs_WeakleyTideRunoff_modified_ALLGENES)
# First, merge the first two data frames: LongTide_vs_WeakleyTideRunoff_modified_ALLGENES + EnsemblID_GO_ID_df
LTvsWTR_merged_df_AG <- merge(LongTide_vs_WeakleyTideRunoff_modified_ALLGENES, EnsemblID_GO_ID_df, by = c("EnsemblID"))
# Then, merge the result with the third data frame:
LongTide_vs_WeakleyTideRunoff_modified_merged_ALLGENES <- merge(LTvsWTR_merged_df_AG, EnsemblID_Gene_name_df, by = c("EnsemblID"))
colnames(LongTide_vs_WeakleyTideRunoff_modified_merged_ALLGENES)
# Subset the table to only include "EnsemblID" , "log2FoldChange", "padj", "Gene_name", "GO_ID"
LongTide_vs_WeakleyTideRunoff_modified_merged_ALLGENES_subset <- subset(LongTide_vs_WeakleyTideRunoff_modified_merged_ALLGENES, select = c(EnsemblID, log2FoldChange, padj, Gene_name, GO_ID)) #you are specifying the columns you want to keep
colnames(LongTide_vs_WeakleyTideRunoff_modified_merged_ALLGENES_subset)
# Extract the table as a .csv
write.table(LongTide_vs_WeakleyTideRunoff_modified_merged_ALLGENES_subset, "LongTide_vs_WeakleyTideRunoff_ALLGENES_pairwisejoined.csv", sep = ",", quote = F)
# Filter the data to make sure there are no duplicate rows (rows may contain the same gene, which is fine, but the GO_ID must be different)
LongTide_vs_WeakleyTideRunoff_modified_merged_ALLGENES_subset_filtered <- LongTide_vs_WeakleyTideRunoff_modified_merged_ALLGENES_subset %>%
  distinct()
# Extract the filtered table as a .csv
write.table(LongTide_vs_WeakleyTideRunoff_modified_merged_ALLGENES_subset_filtered, "LongTide_vs_WeakleyTideRunoff_ALLGENES_pairwisejoined_filtered.csv", sep = ",", quote = F)
# originally had 2182834 rows; now it has 98539 



#### WeakleyTide vs. WeakleyTideRunoff ####

# Use notepad++ to manually add a column before EnsemblID to shift all of the data over and save the original .tsv as a .txt file
WeakleyTide_vs_WeakleyTideRunoff_modified_ALLGENES<-read.table("C:\\Users\\antre\\Box\\Bernal_lab\\Antrelle\\GulfKillifish\\weakleytide_vs_weakleytiderunoff_pairwise_gkf_ALLGENES.txt", header = TRUE, row.names = 1)
colnames(WeakleyTide_vs_WeakleyTideRunoff_modified_ALLGENES)
# First, merge the first two data frames: WeakleyTide_vs_WeakleyTideRunoff_modified_ALLGENES + EnsemblID_GO_ID_df
WTvsWTR_merged_df_AG <- merge(WeakleyTide_vs_WeakleyTideRunoff_modified_ALLGENES, EnsemblID_GO_ID_df, by = c("EnsemblID"))
# Then, merge the result with the third data frame:
WeakleyTide_vs_WeakleyTideRunoff_modified_merged_ALLGENES <- merge(WTvsWTR_merged_df_AG, EnsemblID_Gene_name_df, by = c("EnsemblID"))
colnames(WeakleyTide_vs_WeakleyTideRunoff_modified_merged_ALLGENES)
# Subset the table to only include "EnsemblID" , "log2FoldChange", "padj", "Gene_name", "GO_ID"
WeakleyTide_vs_WeakleyTideRunoff_modified_merged_ALLGENES_subset <- subset(WeakleyTide_vs_WeakleyTideRunoff_modified_merged_ALLGENES, select = c(EnsemblID, log2FoldChange, padj, Gene_name, GO_ID)) #you are specifying the columns you want to keep
colnames(WeakleyTide_vs_WeakleyTideRunoff_modified_merged_ALLGENES_subset)
# Extract the table as a .csv
write.table(WeakleyTide_vs_WeakleyTideRunoff_modified_merged_ALLGENES_subset, "WeakleyTide_vs_WeakleyTideRunoff_ALLGENES_pairwisejoined.csv", sep = ",", quote = F)
# Filter the data to make sure there are no duplicate rows (rows may contain the same gene, which is fine, but the GO_ID must be different)
WeakleyTide_vs_WeakleyTideRunoff_modified_merged_ALLGENES_subset_filtered <- WeakleyTide_vs_WeakleyTideRunoff_modified_merged_ALLGENES_subset %>%
  distinct()
# Extract the filtered table as a .csv
write.table(WeakleyTide_vs_WeakleyTideRunoff_modified_merged_ALLGENES_subset_filtered, "WeakleyTide_vs_WeakleyTideRunoff_ALLGENES_pairwisejoined_filtered.csv", sep = ",", quote = F)
# originally had 2182833 rows; now it has 98539



#### LongTideRunoff vs. WeakleyTideRunoff ####

# Use notepad++ to manually add a column before EnsemblID to shift all of the data over and save the original .tsv as a .txt file
LongTideRunoff_vs_WeakleyTideRunoff_modified_ALLGENES<-read.table("C:\\Users\\antre\\Box\\Bernal_lab\\Antrelle\\GulfKillifish\\longtiderunoff_vs_weakleytiderunoff_pairwise_gkf_ALLGENES.txt", header = TRUE, row.names = 1)
colnames(LongTideRunoff_vs_WeakleyTideRunoff_modified_ALLGENES)
# First, merge the first two data frames: LongTideRunoff_vs_WeakleyTideRunoff_modified_ALLGENES + EnsemblID_GO_ID_df
LTRvsWTR_merged_df_AG <- merge(LongTideRunoff_vs_WeakleyTideRunoff_modified_ALLGENES, EnsemblID_GO_ID_df, by = c("EnsemblID"))
# Then, merge the result with the third data frame:
LongTideRunoff_vs_WeakleyTideRunoff_modified_merged_ALLGENES <- merge(LTRvsWTR_merged_df_AG, EnsemblID_Gene_name_df, by = c("EnsemblID"))
colnames(LongTideRunoff_vs_WeakleyTideRunoff_modified_merged_ALLGENES)
# Subset the table to only include "EnsemblID" , "log2FoldChange", "padj", "Gene_name", "GO_ID"
LongTideRunoff_vs_WeakleyTideRunoff_modified_merged_ALLGENES_subset <- subset(LongTideRunoff_vs_WeakleyTideRunoff_modified_merged_ALLGENES, select = c(EnsemblID, log2FoldChange, padj, Gene_name, GO_ID)) #you are specifying the columns you want to keep
colnames(LongTideRunoff_vs_WeakleyTideRunoff_modified_merged_ALLGENES_subset)
# Extract the table as a .csv
write.table(LongTideRunoff_vs_WeakleyTideRunoff_modified_merged_ALLGENES_subset, "LongTideRunoff_vs_WeakleyTideRunoff_ALLGENES_pairwisejoined.csv", sep = ",", quote = F)
# Filter the data to make sure there are no duplicate rows (rows may contain the same gene, which is fine, but the GO_ID must be different)
LongTideRunoff_vs_WeakleyTideRunoff_modified_merged_ALLGENES_subset_filtered <- LongTideRunoff_vs_WeakleyTideRunoff_modified_merged_ALLGENES_subset %>%
  distinct()
# Extract the filtered table as a .csv
write.table(LongTideRunoff_vs_WeakleyTideRunoff_modified_merged_ALLGENES_subset_filtered, "LongTideRunoff_vs_WeakleyTideRunoff_ALLGENES_pairwisejoined_filtered.csv", sep = ",", quote = F)
# originally had 2182833 rows; now it has 98539



#### LongTideRunoff vs. WeakleyTide ####

# Use notepad++ to manually add a column before EnsemblID to shift all of the data over and save the original .tsv as a .txt file
LongTideRunoff_vs_WeakleyTide_modified_ALLGENES<-read.table("C:\\Users\\antre\\Box\\Bernal_lab\\Antrelle\\GulfKillifish\\longtiderunoff_vs_weakleytide_pairwise_gkf_ALLGENES.txt", header = TRUE, row.names = 1)
colnames(LongTideRunoff_vs_WeakleyTide_modified_ALLGENES)
# First, merge the first two data frames: LongTideRunoff_vs_WeakleyTide_modified_ALLGENES + EnsemblID_GO_ID_df
LTRvsWT_merged_df_AG <- merge(LongTideRunoff_vs_WeakleyTide_modified_ALLGENES, EnsemblID_GO_ID_df, by = c("EnsemblID"))
# Then, merge the result with the third data frame:
LongTideRunoff_vs_WeakleyTide_modified_merged_ALLGENES <- merge(LTRvsWT_merged_df_AG, EnsemblID_Gene_name_df, by = c("EnsemblID"))
colnames(LongTideRunoff_vs_WeakleyTide_modified_merged_ALLGENES)
# Subset the table to only include "EnsemblID" , "log2FoldChange", "padj", "Gene_name", "GO_ID"
LongTideRunoff_vs_WeakleyTide_modified_merged_ALLGENES_subset <- subset(LongTideRunoff_vs_WeakleyTide_modified_merged_ALLGENES, select = c(EnsemblID, log2FoldChange, padj, Gene_name, GO_ID)) #you are specifying the columns you want to keep
colnames(LongTideRunoff_vs_WeakleyTide_modified_merged_ALLGENES_subset)
# Extract the table as a .csv
write.table(LongTideRunoff_vs_WeakleyTide_modified_merged_ALLGENES_subset, "LongTideRunoff_vs_WeakleyTide_ALLGENES_pairwisejoined.csv", sep = ",", quote = F)
# Filter the data to make sure there are no duplicate rows (rows may contain the same gene, which is fine, but the GO_ID must be different)
LongTideRunoff_vs_WeakleyTide_modified_merged_ALLGENES_subset_filtered <- LongTideRunoff_vs_WeakleyTide_modified_merged_ALLGENES_subset %>%
  distinct()
# Extract the filtered table as a .csv
write.table(LongTideRunoff_vs_WeakleyTide_modified_merged_ALLGENES_subset_filtered, "LongTideRunoff_vs_WeakleyTide_ALLGENES_pairwisejoined_filtered.csv", sep = ",", quote = F)
# originally had 282833 rows; now it has 98539











# Proceed to GO_MWU_preparation_GFK.R







