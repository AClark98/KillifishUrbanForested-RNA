#################DESeq Analysis for Gulf Killifish - RNA-seq #####################

# Antrelle D. Clark - Auburn University Department of Biological Sciences 2024
# R-script version 4.4.1 (2024.09.0+375)


setwd("/Users/antre/Box/Bernal_lab/Antrelle/GulfKillifish")

library(BiocManager)
BiocManager::install("DEGreport")
library(DEGreport)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(stringr)
install.packages("svglite")
library(svglite)


# Input files containing gene data
gulfkillifishcountdata<-as.matrix(read.csv("gene_count_matrix - copy.csv", row.names="gene_id"))
gulfkillifishcoldata<-read.csv("sample_info_GKF.csv", header = TRUE, row.names = 1)

# Use the "coldata$___" command to to convert the treatment and origin columns to a **factor**.
# These columns are most likely stored as a character vector (text) or numeric vector, depending on how it was
# originally loaded or created.
# **Assigning back**:
# By assigning `as.factor(adultallcoldata$treatment)` back to `adultallcoldata$treatment`, this code overwrites the 
# original `treatment` column with its factor-converted version.
# Why Convert to a Factor?
# Factors are especially useful in statistical analyses and modeling in R. Many R packages (including DESeq2) use 
# factors to interpret categorical variables correctly. For example, DESeq2 uses factors to recognize and handle 
# experimental conditions like `treatment` groups, which are often needed for differential expression analysis.
# Without converting `treatment` to a factor, the program might not recognize it as a categorical variable, leading 
# to potential errors or misinterpretations in the analysis.
gulfkillifishcoldata$treatment<-as.factor(gulfkillifishcoldata$treatment)
gulfkillifishcoldata$origin<-as.factor(gulfkillifishcoldata$origin)

# Check that the columns and rows have the same name
all(rownames(gulfkillifishcoldata) == colnames(gulfkillifishcountdata))
[1] TRUE

# Run a likelihood ratio test to determine if the effect of origin:treatment, origin, or treatment is significant
GFKdds<-DESeqDataSetFromMatrix(countData = gulfkillifishcountdata, colData = gulfkillifishcoldata, design =~ origin + treatment)
GFKdds<-DESeq(GFKdds)

# This gives you the number of genes due to origin of population: long_bayou_F (forested) vs weakley_bayou_U (urbanized)
originLRT<-DESeq(GFKdds, test = "LRT", reduced = ~treatment)
originLRT.res<-results(originLRT, alpha = 0.05)
summary(originLRT.res)

out of 22430 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 12, 0.053%
LFC < 0 (down)     : 15, 0.067%
outliers [1]       : 0, 0%
low counts [2]     : 5695, 25%
(mean count < 17)

# This gives you the number of genes due to the type of treatment: tide vs tide_runoff
treatmentLRT<-DESeq(GFKdds, test = "LRT", reduced = ~origin)
treatmentLRT.res<-results(treatmentLRT, alpha = 0.05)
summary(treatmentLRT.res)

out of 22430 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 37, 0.16%
LFC < 0 (down)     : 21, 0.094%
outliers [1]       : 0, 0%
low counts [2]     : 6128, 27%
(mean count < 21)

# This gives you the number of genes due to both origin and treatment
originandtreatmentLRT <- DESeq(GFKdds, test = "LRT", reduced = ~ 1)
originandtreatmentLRT.res <-results(originandtreatmentLRT, alpha = 0.05)
summary(originandtreatmentLRT.res)

out of 22430 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 58, 0.26%
LFC < 0 (down)     : 54, 0.24%
outliers [1]       : 0, 0%
low counts [2]     : 3092, 14%
(mean count < 2)

# Output the LFC for (1) origin, (2) treatment, (3) origin + treatment
DESeq2::plotMA(originLRT.res, main = "Log2-Fold Change vs. Base Mean", ylim=c(-5,5))
DESeq2::plotMA(treatmentLRT.res, main = "Log2-Fold Change vs. Base Mean", ylim=c(-5,5))
DESeq2::plotMA(originandtreatmentLRT.res, main = "Log2-Fold Change vs. Base Mean", ylim=c(-5,5))
#^ saved as "LFCorigin", "LFCtreatment", and "LFCoriginandtreatment"

# Save the significantly differentially expressed genes results for (1) origin, (2) treatment, (3) origin + treatment
#########write.csv(originres_sigs, file = "origin_gene_expression_full.csv")
originres_sigs<-subset(originLRT.res, padj < 0.05)
summary(originres_sigs)

  out of 27 with nonzero total read count
  adjusted p-value < 0.05
  LFC > 0 (up)       : 12, 44%
  LFC < 0 (down)     : 15, 56%
  outliers [1]       : 0, 0%
  low counts [2]     : 0, 0%
  (mean count < 17)
  
write.csv(originLRT.res, file = "origin_gene_expression_results_sigs.csv")


treatmentres_sigs<-subset(treatmentLRT.res, padj < 0.05)
summary(treatmentres_sigs)

  out of 58 with nonzero total read count
  adjusted p-value < 0.05
  LFC > 0 (up)       : 37, 64%
  LFC < 0 (down)     : 21, 36%
  outliers [1]       : 0, 0%
  low counts [2]     : 0, 0%
  (mean count < 21)
  
write.csv(treatmentres_sigs, file = "treatment_gene_expression_results_sigs.csv")


originandtreatmentres_sigs<-subset(originandtreatmentLRT.res, padj < 0.05)
summary(originandtreatmentres_sigs)

  out of 112 with nonzero total read count
  adjusted p-value < 0.05
  LFC > 0 (up)       : 58, 52%
  LFC < 0 (down)     : 54, 48%
  outliers [1]       : 0, 0%
  low counts [2]     : 0, 0%
  (mean count < 2)
  
write.csv(originandtreatmentres_sigs, file = "originandtreatment_gene_expression_results_sigs.csv")

# LRT and PCA: origin
originLRTvsd<-vst(originLRT)
# We want just the significant genes
originLRTsigs<-subset(originLRT.res, padj < 0.05)
origin_sigs<-rownames(originLRTsigs)
VSD.origin.subset<- originLRTvsd[rownames(originLRTvsd) %in% origin_sigs, ]
summary(VSD.origin.subset)
[1] "DESeqTransform object of length 27 with 24 metadata columns"
plotPCA(VSD.origin.subset, intgroup = "origin")
gulfkillifish_PCA_dim<-plotPCA(VSD.origin.subset, intgroup = "origin", returnData = T) #save
ggplot(gulfkillifish_PCA_dim, aes(x = PC1, y = PC2, color = origin)) +
geom_point(size = 3) +
labs(x = "PC1: 47% variance", y = "PC2: 18% variance", title = "Sample Distance Matrix by Origin of Population") +
theme_minimal()
############### saved as a "svg" and "png" under SDMorigin" ################### can click "export and save as image" or "ggsave"

# LRT and PCA: treatment
treatmentLRTvsd<-vst(treatmentLRT)
# We want just the significant genes
treatmentLRTsigs<-subset(treatmentLRT.res, padj < 0.05)
treatment_sigs<-rownames(treatmentLRTsigs)
VSD.treatment.subset<- treatmentLRTvsd[rownames(treatmentLRTvsd) %in% treatment_sigs, ]
summary(VSD.treatment.subset)
[1] "DESeqTransform object of length 58 with 24 metadata columns"
plotPCA(VSD.treatment.subset, intgroup = "treatment")
gulfkillifish_PCA_t_dim<-plotPCA(VSD.treatment.subset, intgroup = "treatment", returnData = T) #save
ggplot(gulfkillifish_PCA_t_dim, aes(x = PC1, y = PC2, color = treatment)) +
geom_point(size =3) +
labs(x = "PC1: 51% variance", y = "PC2: 11% variance", title = "Sample Distance Matrix by Treatment") +
theme_minimal()
############### saved as a "svg" and "png" under SDMtreatment" ###################

# LRT and PCA: origin and treatment
originandtreatmentLRTvsd<-vst(originandtreatmentLRT)
# We want just the significant genes
originandtreatmentLRTsigs<-subset(originandtreatmentLRT.res, padj < 0.05)
originandtreatment_sigs<-rownames(originandtreatmentLRTsigs)
VSD.originandtreatment.subset<- originandtreatmentLRTvsd[rownames(originandtreatmentLRTvsd) %in% originandtreatment_sigs, ]
summary(VSD.originandtreatment.subset)
[1] "DESeqTransform object of length 112 with 24 metadata columns"
plotPCA(VSD.originandtreatment.subset, intgroup = c("origin", "treatment"))
gulfkillifish_PCA_o_t_dim<-plotPCA(VSD.originandtreatment.subset, intgroup = c("origin", "treatment"), returnData = T) #save
install.packages("cluster")
library(cluster)
install.packages("factoextra")
library(factoextra)
# Perform k-means clustering (let's assume 4 clusters)
set.seed(123)  # For reproducibility
kmeans_result <- kmeans(gulfkillifish_PCA_o_t_dim[, c("PC1", "PC2")], centers = 4)

# Add cluster information to the PCA data
gulfkillifish_PCA_o_t_dim$Cluster <- as.factor(kmeans_result$cluster)
# create PCA
ggplot(gulfkillifish_PCA_o_t_dim, aes(x = PC1, y = PC2, color = origin, shape = treatment)) +
geom_point(size = 3) +
stat_ellipse(aes(color = cluster), level = 0.95) +
labs(x = "PC1: 31% variance", y = "PC2: 13% variance", title = "Sample Distance Matrix by Origin and Treatment") +
theme_minimal()


# this is to add the ellipses to show the clustering of the groups

# LRT and PCA: origin and treatment
originandtreatmentLRTvsd <- vst(originandtreatmentLRT)

# We want just the significant genes
originandtreatmentLRTsigs <- subset(originandtreatmentLRT.res, padj < 0.05)
originandtreatment_sigs <- rownames(originandtreatmentLRTsigs)

VSD.originandtreatment.subset <- originandtreatmentLRTvsd[rownames(originandtreatmentLRTvsd) %in% originandtreatment_sigs, ]
summary(VSD.originandtreatment.subset)

# PCA plot using plotPCA to get the PCA data for plotting
gulfkillifish_PCA_o_t_dim <- plotPCA(VSD.originandtreatment.subset, intgroup = c("origin", "treatment"), returnData = TRUE) # save PCA data

# Install and load necessary packages
install.packages("cluster")
library(cluster)
install.packages("factoextra")
library(factoextra)

# Combine origin and treatment into a single variable for color coding
gulfkillifish_PCA_o_t_dim$origin_treatment <- paste(gulfkillifish_PCA_o_t_dim$origin, gulfkillifish_PCA_o_t_dim$treatment, sep = "_")

# Now plot with circles as points, 4 colors for origin+treatment combinations, and add ellipses
ggplot(gulfkillifish_PCA_o_t_dim, aes(x = PC1, y = PC2, color = origin_treatment)) +
  geom_point(size = 3, shape = 16) +  # Keep points as circles (shape = 16)
  stat_ellipse(aes(color = origin_treatment), level = 0.95) +  # Add ellipses based on origin+treatment
  labs(x = "PC1: 31% variance", y = "PC2: 13% variance", title = "Sample Distance Matrix by Origin and Treatment with Clusters") +
  theme_minimal()


# color-blind friendly palette
library(RColorBrewer)

# Choose a color-blind friendly palette
palette <- brewer.pal(4, "Set1")  # Set1 is a good option with 4 distinct colors

# Create the PCA plot using the new palette
ggplot(gulfkillifish_PCA_o_t_dim, aes(x = PC1, y = PC2, color = origin_treatment)) +
  geom_point(size = 3, shape = 16) +  # Keep points as circles (shape = 16)
  stat_ellipse(aes(color = origin_treatment), level = 0.95) +  # Add ellipses based on origin+treatment
  scale_color_manual(values = palette) +  # Apply the color palette
  labs(x = "PC1: 31% variance", y = "PC2: 13% variance", title = "Sample Distance Matrix by Origin and Treatment with Clusters") +
  theme_minimal()

######## Make PCA with different shapes which is easier to differentiate for color-blind friendliness ##########


# LRT and PCA: origin and treatment
originandtreatmentLRTvsd <- vst(originandtreatmentLRT)

# We want just the significant genes
originandtreatmentLRTsigs <- subset(originandtreatmentLRT.res, padj < 0.05)
originandtreatment_sigs <- rownames(originandtreatmentLRTsigs)

VSD.originandtreatment.subset <- originandtreatmentLRTvsd[rownames(originandtreatmentLRTvsd) %in% originandtreatment_sigs, ]
summary(VSD.originandtreatment.subset)

# PCA plot using plotPCA to get the PCA data for plotting
gulfkillifish_PCA_o_t_dim <- plotPCA(VSD.originandtreatment.subset, intgroup = c("origin", "treatment"), returnData = TRUE) # save PCA data

# Create a new column for origin_treatment combinations (to color by origin+treatment)
gulfkillifish_PCA_o_t_dim$origin_treatment <- paste(gulfkillifish_PCA_o_t_dim$origin, gulfkillifish_PCA_o_t_dim$treatment, sep = "_")

# Assign numeric shape values for each origin_treatment combination
# Numeric values for shapes: 16=Circle, 17=Triangle, 15=Square, 8=Star (will be replaced with text stars)
gulfkillifish_PCA_o_t_dim$shape_group <- factor(gulfkillifish_PCA_o_t_dim$origin_treatment,
                                                levels = unique(gulfkillifish_PCA_o_t_dim$origin_treatment),
                                                labels = c(16, 17, 15, 8))  # Assign numeric shape codes here

# Now, we replace shape = 8 (Star) with the actual Unicode star symbol using geom_text()
# Plot without geom_point() for shape = 8, and only show stars as text
ggplot(gulfkillifish_PCA_o_t_dim, aes(x = PC1, y = PC2, color = origin_treatment, shape = shape_group)) +
  # Use geom_point for other shapes (circle, triangle, square)
  geom_point(data = subset(gulfkillifish_PCA_o_t_dim, shape_group != 8), aes(shape = shape_group), size = 3) +  # Regular shapes for points
  # Add ellipses based on origin+treatment
  stat_ellipse(aes(color = origin_treatment), level = 0.95) +  
  labs(x = "PC1: 31% variance", y = "PC2: 13% variance", title = "Sample Distance Matrix by Origin and Treatment with Clusters") +
  theme_gray() +  # Set the gray background theme
  coord_fixed(ratio = 1) +  # Ensure aspect ratio is fixed to keep shapes proportional
  scale_shape_manual(values = c(16, 17, 15, 8)) +  # Ensure manual assignment of shapes to levels
  # Replace the star shape with actual Unicode stars (and remove original shape)
  geom_text(aes(label = ifelse(shape_group == 8, "★", "")), size = 10, color = "purple")  # Add Unicode stars for shape 8




############### saved as a "svg" and "png" under SDMoriginandtreatment" ###################

# Produce a HEATMAP to identify any outliers using the combination: origin and treatment
library(pheatmap)
library(RColorBrewer)
library(ComplexHeatmap)
sampleDists <- dist(t(assay(VSD.originandtreatment.subset)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- c("L1C1", "L1C2", "L1C3", "L1C4", "L1C5", "L2U1", "L2U2", "L2U3", "L2U4", "L2U5", "L3C1", "L3C2", "L3C3", "L3C4", "L3C5", "L4U1", "L4U2", "L4U3", "L4U4", "L4U5", "W1C1", "W1C2", "W1C3", "W1C4", "W2U1", "W2U2", "W2U3", "W2U4", "W2U5", "W3C1", "W3C2", "W3C3", "W3C4", "W4U1", "W4U2", "W4U3", "W4U4")
colnames(sampleDistMatrix) <- c("tide_F", "tide_F", "tide_F", "tide_F", "tide_F", "tide_runoff_F", "tide_runoff_F", "tide_runoff_F", "tide_runoff_F", "tide_runoff_F", "tide_F", "tide_F", "tide_F", "tide_F", "tide_F", "tide_runoff_F", "tide_runoff_F", "tide_runoff_F", "tide_runoff_F", "tide_runoff_F", "tide_U", "tide_U", "tide_U", "tide_U", "tide_runoff_U", "tide_runoff_U", "tide_runoff_U", "tide_runoff_U", "tide_runoff_U", "tide_U", "tide_U", "tide_U", "tide_U", "tide_runoff_U", "tide_runoff_U", "tide_runoff_U", "tide_runoff_U")
pheatmap(sampleDistMatrix,
cluster_rows = TRUE,
cluster_cols = TRUE,
display_numbers = FALSE, ##### or you can do "display_numbers = TRUE" to add numbers to the heatmap #####
main = "Heatmap crossing origin and treatment")

### saved as .svg and .png under "pHM_originandtreatment_displaynumbers" and "pHM_originandtreatment" ###

# Pairwise-contrast comparisons (comparisons between two individual treatments:
# (1) Long + Tide VS. Weakley + Tide
# (2) Long + Tide_runoff VS. Weakley + Tide_runoff
# (3) Long + Tide VS. Long + Tide_runoff
# (4) Weakley + Tide VS. Weakley + Tide_runoff
# (5) Long + Tide VS. Weakley + Tide_runoff
# (6) Long + Tide_runoff VS. Weakley + Tide


# (1) Long + Tide VS. Weakley + Tide
LongTide_WeakleyTidedds<-DESeqDataSetFromMatrix(countData = gulfkillifishcountdata, colData = gulfkillifishcoldata, design =~ combined_for_pairwise)
LongTide_WeakleyTidedds<-DESeq(LongTide_WeakleyTidedds)
# To get all genes
LongTide_vs_WeakleyTide_ALL <- results(LongTide_WeakleyTidedds, contrast = c("combined_for_pairwise", "long_tide", "weakley_tide"))
LongTide_vs_WeakleyTide_ALLGENES <- as.data.frame(LongTide_vs_WeakleyTide_ALL)
# Write output files to be used in GO term tests
write.table(LongTide_vs_WeakleyTide_ALLGENES, "longtide_vs_weakleytide_pairwise_gkf_ALLGENES.tsv", sep = "\t", quote = F)
# To get significant genes
LongTide_vs_WeakleyTide <- results(LongTide_WeakleyTidedds, alpha = 0.05, contrast = c("combined_for_pairwise", "long_tide", "weakley_tide"))
LongTide_vs_WeakleyTide
## ^log fold change and Wald test p-value will show up; next we must adjust just to get significant ##
LongTide_vs_WeakleyTide_sigs<- subset(LongTide_vs_WeakleyTide, padj < 0.05)
summary(LongTide_vs_WeakleyTide_sigs)
#OUTPUT:
  out of 6 with nonzero total read count
  adjusted p-value < 0.05
  LFC > 0 (up)       : 4, 67%
  LFC < 0 (down)     : 2, 33%
  outliers [1]       : 0, 0%
  low counts [2]     : 0, 0%
  (mean count < 0)
# Write output files to be used in GO term tests
LongTide_vs_WeakleyTide_sigs<-as.data.frame(LongTide_vs_WeakleyTide_sigs)
write.table(LongTide_vs_WeakleyTide_sigs, "longtide_vs_weakleytide_pairwise_gkf.tsv", sep = "\t", quote = F) #the "sig genes" will be used to identify shared genes between groups, while the "all genes" will be used for GO term analysis

# Create a PCA for the LT vs WT PWC
LongTide_WeakleyTidedds<-DESeqDataSetFromMatrix(countData = gulfkillifishcountdata, colData = gulfkillifishcoldata, design =~ combined_for_pairwise)
LongTide_WeakleyTidedds<-DESeq(LongTide_WeakleyTidedds)
LTvsWTdds_vst <- vst(LongTide_WeakleyTidedds, blind = FALSE)
# Subset the DESeqDataSet to include only the "long_tide" and "weakley_tide" treatments
LTvsWT_PWC_subsetdata <- subset(colData(LongTide_WeakleyTidedds), combined_for_pairwise %in% c("long_tide", "weakley_tide"))
# Check the distribution of conditions in the subset
table(LongTide_WeakleyTidedds_subset$combined_for_pairwise)
#Drop the unused levels
LongTide_WeakleyTidedds$combined_for_pairwise <- droplevels(LongTide_WeakleyTidedds$combined_for_pairwise)
LongTide_WeakleyTidedds_subset <- LongTide_WeakleyTidedds[, rownames(LTvsWT_PWC_subsetdata)]
# Check the levels after dropping unused ones
levels(LongTide_WeakleyTidedds$combined_for_pairwise)
# Subset the data to only include the two treatments
subset_data <- subset(colData(LongTide_WeakleyTidedds), combined_for_pairwise %in% c("long_tide", "weakley_tide"))
# Subset the DESeqDataSet using the rownames of the subset data
LongTide_WeakleyTidedds_subset <- LongTide_WeakleyTidedds[, rownames(subset_data)]
# Re-check the distribution of conditions after subsetting
table(subset_data$combined_for_pairwise)
# Run DESeq2 with a simple design formula (no additional covariates)
LongTide_WeakleyTidedds_subset <- DESeqDataSetFromMatrix(
countData = gulfkillifishcountdata[, rownames(subset_data)],
colData = subset_data,
design = ~combined_for_pairwise
)
# Run DESeq2 analysis
LongTide_WeakleyTidedds_subset <- DESeq(LongTide_WeakleyTidedds_subset)
# Extract results and filter by padj < 0.05
LTvsWT_results <- results(LongTide_WeakleyTidedds_subset, alpha = 0.05, contrast = c("combined_for_pairwise", "long_tide", "weakley_tide"))
# Filter for genes with padj < 0.05
LTvsWT_sig_results <- LTvsWT_results[LTvsWT_results$padj < 0.05, ]
# Apply VST to the subsetted DESeqDataSet
LTvsWTdds_vst_subset <- vst(LongTide_WeakleyTidedds_subset, blind = FALSE)
# Plot PCA for the subsetted VST data
plotPCA(LTvsWTdds_vst_subset, intgroup = "combined_for_pairwise")




# (2) Long + Tide_runoff VS. Weakley + Tide_runoff
LongTideRunoff_WeakleyTideRunoffdds<-DESeqDataSetFromMatrix(countData = gulfkillifishcountdata, colData = gulfkillifishcoldata, design =~ combined_for_pairwise)
LongTideRunoff_WeakleyTideRunoffdds<-DESeq(LongTideRunoff_WeakleyTideRunoffdds)
# To get all genes
LongTideRunoff_vs_WeakleyTideRunoff_ALL <- results(LongTideRunoff_WeakleyTideRunoffdds, contrast = c("combined_for_pairwise", "long_tiderunoff", "weakley_tiderunoff"))
LongTideRunoff_vs_WeakleyTideRunoff_ALLGENES <- as.data.frame(LongTideRunoff_vs_WeakleyTideRunoff_ALL)
# Write output files to be used in GO term tests
write.table(LongTideRunoff_vs_WeakleyTideRunoff_ALLGENES, "longtiderunoff_vs_weakleytiderunoff_pairwise_gkf_ALLGENES.tsv", sep = "\t", quote = F)
# To get significant genes
LongTideRunoff_WeakleyTideRunoff <- results(LongTideRunoff_WeakleyTideRunoffdds, alpha = 0.05, contrast = c("combined_for_pairwise", "long_tiderunoff", "weakley_tiderunoff"))
LongTideRunoff_WeakleyTideRunoff
## ^log fold change and Wald test p-value will show up; next we must adjust just to get significant ##
LongTideRunoff_WeakleyTideRunoff_sigs<- subset(LongTideRunoff_WeakleyTideRunoff, padj < 0.05)
summary(LongTideRunoff_WeakleyTideRunoff_sigs)
#OUTPUT:
  out of 294 with nonzero total read count
  adjusted p-value < 0.05
  LFC > 0 (up)       : 170, 58%
  LFC < 0 (down)     : 124, 42%
  outliers [1]       : 0, 0%
  low counts [2]     : 0, 0%
  (mean count < 21)
LongTideRunoff_WeakleyTideRunoff_sigs<-as.data.frame(LongTideRunoff_WeakleyTideRunoff_sigs)
write.table(LongTideRunoff_WeakleyTideRunoff_sigs, "longtiderunoff_vs_weakleytiderunoff_pairwise_gkf.tsv", sep = "\t", quote = F) #the "sig genes" will be used to identify shared genes between groups, while the "all genes" will be used for GO term analysis
  
#Create PCA for PWC LTR vs WTR
LongTR_WeakleyTRdds<-DESeqDataSetFromMatrix(countData = gulfkillifishcountdata, colData = gulfkillifishcoldata, design =~ combined_for_pairwise)
LongTR_WeakleyTRdds<-DESeq(LongTR_WeakleyTRdds)
LTRvsWTRdds_vst <- vst(LongTR_WeakleyTRdds, blind = FALSE)
# Subset the DESeqDataSet to include only the "long_tiderunoff" and "weakley_tiderunoff" treatments
LTRvsWTR_PWC_subsetdata <- subset(colData(LongTR_WeakleyTRdds), combined_for_pairwise %in% c("long_tiderunoff", "weakley_tiderunoff"))
# Check the distribution of conditions in the subset
table(LTRvsWTR_PWC_subsetdata$combined_for_pairwise)
#Drop the unused levels
LongTR_WeakleyTRdds$combined_for_pairwise <- droplevels(LongTR_WeakleyTRdds$combined_for_pairwise)
LongTR_WeakleyTRdds_subset <- LongTR_WeakleyTRdds[, rownames(LTRvsWTR_PWC_subsetdata)]
# Check the levels after dropping unused ones
levels(LongTR_WeakleyTRdds$combined_for_pairwise)
# Subset the data to only include the two treatments
LTRvsWTRsubset_data <- subset(colData(LongTR_WeakleyTRdds), combined_for_pairwise %in% c("long_tiderunoff", "weakley_tiderunoff"))
# Subset the DESeqDataSet using the rownames of the subset data
LongTide_WeakleyTidedds_subset <- LongTide_WeakleyTidedds[, rownames(LTRvsWTRsubset_data)]
# Re-check the distribution of conditions after subsetting
table(subset_data$combined_for_pairwise)
# Run DESeq2 with a simple design formula (no additional covariates)
LongTR_WeakleyTRdds_subset <- DESeqDataSetFromMatrix(
  countData = gulfkillifishcountdata[, rownames(LTRvsWTRsubset_data)],
  colData = LTRvsWTRsubset_data,
  design = ~combined_for_pairwise
)
# Run DESeq2 analysis
LongTR_WeakleyTRdds_subset <- DESeq(LongTR_WeakleyTRdds_subset)
# Apply VST to the subsetted DESeqDataSet
LTRvsWTRdds_vst_subset <- vst(LongTR_WeakleyTRdds_subset, blind = FALSE)
# Plotting PCA with customized colors based on the 'combined_for_pairwise' condition
library(ggplot2)
# Step 1: Create a PCA plot and specify the color variable
pca_plot <- plotPCA(LTRvsWTRdds_vst_subset, intgroup = "combined_for_pairwise", returnData = TRUE)
# Step 2: Modify the plot with custom colors using ggplot2
pca_plot <- ggplot(pca_plot, aes(x = PC1, y = PC2, color = combined_for_pairwise)) +
  geom_point(size = 3) + # Adjust point size
  scale_color_manual(values = c("long_tiderunoff" = "green", "weakley_tiderunoff" = "purple")) + # Set colors for each condition
  labs(title = "PCA of Significant Genes") +
  theme_minimal()
# Step 3: Display the PCA plot
print(pca_plot)



# Load necessary libraries
library(DESeq2)
library(pheatmap)
library(ggplot2)

# Step 1: Subset the significant genes (those with adjusted p-value < 0.05)
LongTideRunoff_WeakleyTideRunoff_sigs <- subset(LongTideRunoff_WeakleyTideRunoff, padj < 0.05)

# Step 2: Sort the significant genes by log2FoldChange (top 20 and bottom 20)
# Top 20 genes with the highest log2FoldChange (upregulated)
top_20_genes <- LongTideRunoff_WeakleyTideRunoff_sigs[order(LongTideRunoff_WeakleyTideRunoff_sigs$log2FoldChange, decreasing = TRUE), ][1:20, ]

# Bottom 20 genes with the lowest log2FoldChange (downregulated)
bottom_20_genes <- LongTideRunoff_WeakleyTideRunoff_sigs[order(LongTideRunoff_WeakleyTideRunoff_sigs$log2FoldChange, decreasing = FALSE), ][1:20, ]

# Step 3: Combine the top and bottom 20 genes
top_bottom_40_genes <- rbind(top_20_genes, bottom_20_genes)

# Step 4: Get the normalized expression data for the selected genes (VST-transformed)
expression_data <- assay(LTRvsWTRdds_vst_subset)[rownames(top_bottom_40_genes), ]

# Step 5: Generate the heatmap
pheatmap(expression_data,
         scale = "row",  # Normalize by rows (genes) to highlight expression differences
         clustering_distance_rows = "euclidean",  # Cluster rows (genes) based on similarity
         clustering_distance_cols = "euclidean",  # Cluster columns (samples) based on similarity
         main = "Heatmap of Top 20 and Bottom 20 Significant Genes (log2FoldChange)",
         show_rownames = TRUE,  # Show gene names on rows
         show_colnames = TRUE,  # Show sample names on columns
         color = colorRampPalette(c("blue", "white", "red"))(50))  # Customize color scale




#for top 40 and bottom 40 genes
# Step 2: Sort the significant genes by log2FoldChange (top 20 and bottom 20)
# Top 20 genes with the highest log2FoldChange (upregulated)
#top_40_genes <- LongTideRunoff_WeakleyTideRunoff_sigs[order(LongTideRunoff_WeakleyTideRunoff_sigs$log2FoldChange, decreasing = TRUE), ][1:40, ]
# Bottom 20 genes with the lowest log2FoldChange (downregulated)
#bottom_40_genes <- LongTideRunoff_WeakleyTideRunoff_sigs[order(LongTideRunoff_WeakleyTideRunoff_sigs$log2FoldChange, decreasing = FALSE), ][1:40, ]
# Step 3: Combine the top and bottom 20 genes
#top_bottom_80_genes <- rbind(top_40_genes, bottom_40_genes)
# Step 4: Get the normalized expression data for the selected genes (VST-transformed)
#expression_data <- assay(LTRvsWTRdds_vst_subset)[rownames(top_bottom_80_genes), ]
# Step 5: Generate the heatmap
#pheatmap(expression_data,
 #        scale = "row",  # Normalize by rows (genes) to highlight expression differences
#       clustering_distance_rows = "euclidean",  # Cluster rows (genes) based on similarity
  #       clustering_distance_cols = "euclidean",  # Cluster columns (samples) based on similarity
    #     main = "Heatmap of Top 40 and Bottom 40 Significant Genes (log2FoldChange)",
    #     show_rownames = TRUE,  # Show gene names on rows
    #     show_colnames = TRUE,  # Show sample names on columns
    #     color = colorRampPalette(c("blue", "white", "red"))(50))  # Customize color scale




# (3) Long + Tide VS. Long + Tide_runoff
LongTide_LongTideRunoffdds<-DESeqDataSetFromMatrix(countData = gulfkillifishcountdata, colData = gulfkillifishcoldata, design =~ combined_for_pairwise)
LongTide_LongTideRunoffdds<-DESeq(LongTide_LongTideRunoffdds)
# To get all genes
LongTide_vs_LongTideRunoff_ALL <- results(LongTide_LongTideRunoffdds, contrast = c("combined_for_pairwise", "long_tide", "long_tiderunoff"))
LongTide_vs_LongTideRunoff_ALLGENES <- as.data.frame(LongTide_vs_LongTideRunoff_ALL)
# Write output files to be used in GO term tests
write.table(LongTide_vs_LongTideRunoff_ALLGENES, "longtide_vs_longtiderunoff_pairwise_gkf_ALLGENES.tsv", sep = "\t", quote = F)
# To get significant genes
LongTide_LongTideRunoff <- results(LongTide_LongTideRunoffdds, alpha = 0.05, contrast = c("combined_for_pairwise", "long_tide", "long_tiderunoff"))
LongTide_LongTideRunoff
## ^log fold change and Wald test p-value will show up; next we must adjust just to get significant ##
LongTide_LongTideRunoff_sigs<- subset(LongTide_LongTideRunoff, padj < 0.05)
summary(LongTide_LongTideRunoff_sigs)
#OUTPUT:
  out of 9 with nonzero total read count
  adjusted p-value < 0.05
  LFC > 0 (up)       : 5, 56%
  LFC < 0 (down)     : 4, 44%
  outliers [1]       : 0, 0%
  low counts [2]     : 0, 0%
  (mean count < 0)
LongTide_LongTideRunoff_sigs<-as.data.frame(LongTide_LongTideRunoff_sigs)
write.table(LongTide_LongTideRunoff_sigs, "longtide_vs_longtiderunoff_pairwise_gkf.tsv", sep = "\t", quote = F) #the "sig genes" will be used to identify shared genes between groups, while the "all genes" will be used for GO term analysis

  
  
    
# (4) Weakley + Tide VS. Weakley + Tide_runoff
WeakleyTide_WeakleyTideRunoffdds<-DESeqDataSetFromMatrix(countData = gulfkillifishcountdata, colData = gulfkillifishcoldata, design =~ combined_for_pairwise)
WeakleyTide_WeakleyTideRunoffdds<-DESeq(WeakleyTide_WeakleyTideRunoffdds)
# To get all genes
WeakleyTide_vs_WeakleyTideRunoff_ALL <- results(WeakleyTide_WeakleyTideRunoffdds, contrast = c("combined_for_pairwise", "weakley_tide", "weakley_tiderunoff"))
WeakleyTide_vs_WeakleyTideRunoff_ALLGENES <- as.data.frame(WeakleyTide_vs_WeakleyTideRunoff_ALL)
# Write output files to be used in GO term tests
write.table(WeakleyTide_vs_WeakleyTideRunoff_ALLGENES, "weakleytide_vs_weakleytiderunoff_pairwise_gkf_ALLGENES.tsv", sep = "\t", quote = F)
# To get significant genes
WeakleyTide_WeakleyTideRunoff <- results(WeakleyTide_WeakleyTideRunoffdds, alpha = 0.05, contrast = c("combined_for_pairwise", "weakley_tide", "weakley_tiderunoff"))
WeakleyTide_WeakleyTideRunoff
## ^log fold change and Wald test p-value will show up; next we must adjust just to get significant ##
WeakleyTide_WeakleyTideRunoff_sigs<- subset(WeakleyTide_WeakleyTideRunoff, padj < 0.05)
summary(WeakleyTide_WeakleyTideRunoff_sigs)
#OUTPUT:
  out of 199 with nonzero total read count
  adjusted p-value < 0.05
  LFC > 0 (up)       : 96, 48%
  LFC < 0 (down)     : 103, 52%
  outliers [1]       : 0, 0%
  low counts [2]     : 0, 0%
  (mean count < 17)
WeakleyTide_WeakleyTideRunoff_sigs<-as.data.frame(WeakleyTide_WeakleyTideRunoff_sigs)
write.table(WeakleyTide_WeakleyTideRunoff_sigs, "weakleytide_vs_weakleytiderunoff_pairwise_gfk.tsv", sep = "\t", quote = F) #the "sig genes" will be used to identify shared genes between groups, while the "all genes" will be used for GO term analysis
  
  

# (5) Long + Tide VS. Weakley + Tide_runoff
LongTide_WeakleyTideRunoffdds<- DESeqDataSetFromMatrix(countData = gulfkillifishcountdata, colData = gulfkillifishcoldata, design =~ combined_for_pairwise)
LongTide_WeakleyTideRunoffdds<-DESeq(LongTide_WeakleyTideRunoffdds)
# To get all genes
LongTide_vs_WeakleyTideRunoff_ALL <- results(LongTide_WeakleyTideRunoffdds, contrast = c("combined_for_pairwise", "long_tide", "weakley_tiderunoff"))
LongTide_vs_WeakleyTideRunoff_ALLGENES <- as.data.frame(LongTide_vs_WeakleyTideRunoff_ALL)
# Write output files to be used in GO term tests
write.table(LongTide_vs_WeakleyTideRunoff_ALLGENES, "longtide_vs_weakleytiderunoff_pairwise_gkf_ALLGENES.tsv", sep = "\t", quote = F)
# To get significant genes
LongTide_WeakleyTideRunoff <- results(LongTide_WeakleyTideRunoffdds, alpha = 0.05, contrast = c("combined_for_pairwise", "long_tide", "weakley_tiderunoff"))
LongTide_WeakleyTideRunoff
## ^log fold change and Wald test p-value will show up; next we must adjust just to get significant ##
LongTide_WeakleyTideRunoff_sigs<- subset(LongTide_WeakleyTideRunoff, padj < 0.05)
summary(LongTide_WeakleyTideRunoff_sigs)
#OUTPUT:
  out of 69 with nonzero total read count
  adjusted p-value < 0.05
  LFC > 0 (up)       : 29, 42%
  LFC < 0 (down)     : 40, 58%
  outliers [1]       : 0, 0%
  low counts [2]     : 0, 0%
  (mean count < 2)
LongTide_WeakleyTideRunoff_sigs<-as.data.frame(LongTide_WeakleyTideRunoff_sigs)
write.table(LongTide_WeakleyTideRunoff_sigs, "longtide_vs_weakleytiderunoff_pairwise_gfk.tsv", sep = "\t", quote = F) #the "sig genes" will be used to identify shared genes between groups, while the "all genes" will be used for GO term analysis

  
  
  
# (6) Long + Tide_runoff VS. Weakley + Tide
LongTideRunoff_WeakleyTidedds<- DESeqDataSetFromMatrix(countData = gulfkillifishcountdata, colData = gulfkillifishcoldata, design =~ combined_for_pairwise)
LongTideRunoff_WeakleyTidedds<-DESeq(LongTideRunoff_WeakleyTidedds)
# To get all genes
LongTideRunoff_vs_WeakleyTide_ALL <- results(LongTideRunoff_WeakleyTidedds, contrast = c("combined_for_pairwise", "long_tiderunoff", "weakley_tide"))
LongTideRunoff_vs_WeakleyTide_ALLGENES <- as.data.frame(LongTideRunoff_vs_WeakleyTide_ALL)
# Write output files to be used in GO term tests
write.table(LongTideRunoff_vs_WeakleyTide_ALLGENES, "longtiderunoff_vs_weakleytide_pairwise_gkf_ALLGENES.tsv", sep = "\t", quote = F)
# To get significant genes
LongTideRunoff_WeakleyTide <- results(LongTide_WeakleyTideRunoffdds, alpha = 0.05, contrast = c("combined_for_pairwise", "long_tiderunoff", "weakley_tide"))
LongTideRunoff_WeakleyTide  
## ^log fold change and Wald test p-value will show up; next we must adjust just to get significant ##
LongTideRunoff_WeakleyTide_sigs<-subset(LongTideRunoff_WeakleyTide, padj < 0.05)
summary(LongTideRunoff_WeakleyTide_sigs)
#OUTPUT:
  out of 29 with nonzero total read count
  adjusted p-value < 0.05
  LFC > 0 (up)       : 7, 24%
  LFC < 0 (down)     : 22, 76%
  outliers [1]       : 0, 0%
  low counts [2]     : 0, 0%
  (mean count < 0)
LongTideRunoff_WeakleyTide_sigs<-as.data.frame(LongTideRunoff_WeakleyTide_sigs)
write.table(LongTideRunoff_WeakleyTide_sigs, "longtiderunoff_vs_weakleytide_pairwise_gfk.tsv", sep = "\t", quote = F) #the "sig genes" will be used to identify shared genes between groups, while the "all genes" will be used for GO term analysis
  
  


# Create PCA and Heatmap for treatments: (2) Long + Tide_Runoff vs Weakley + Tide_Runoff [294 sigs], (3) Long + Tide vs Long + Tide_Runoff [9 sigs], (4) Weakley + Tide vs Weakley + Tide_Runoff [199 sigs]


# (2) Long + Tide_Runoff vs Weakley + Tide_Runoff [294 sigs]
LTRvsWTRvsd <-vst(LongTideRunoff_WeakleyTideRunoffdds)
# We just want the significant genes
LTRvsWTRsigs<-subset(LongTideRunoff_WeakleyTideRunoff_sigs, padj <0.05)
LTRvsWTR_sigs<-rownames(LTRvsWTRsigs)
VSD.LTRvsWTR.subset<- LTRvsWTRvsd[rownames(LTRvsWTRvsd) %in% LTRvsWTR_sigs, ]
summary(VSD.LTRvsWTR.subset)
[1] "DESeqTransform object of length 294 with 31 metadata columns"

LTRvsWTR_groups_to_keep <- c("long_tiderunoff", "weakley_tiderunoff")
print(LTRvsWTR_groups_to_keep)
[1] "long_tiderunoff"    "weakley_tiderunoff"

LTRvsWTRsubset_samples <- colData(VSD.LTRvsWTR.subset)$combined_for_pairwise %in% LTRvsWTR_groups_to_keep
LTRvsWTRVSD_subset <- VSD.LTRvsWTR.subset[, LTRvsWTRsubset_samples]
# plotPCA and save as .svg and .png - PCALTRvsWTR294sigs_pairwise
plotPCA(LTRvsWTRVSD_subset, intgroup = "combined_for_pairwise")
#ggplot and save as .svg and .png - SDM_LTRvsWTR294sigs_pairwise
LTRvsWTRgulfkillifish_PCA_t_dim<-plotPCA(LTRvsWTRVSD_subset, intgroup = "combined_for_pairwise", returnData = T) #save
ggplot(LTRvsWTRgulfkillifish_PCA_t_dim, aes(x = PC1, y = PC2, color = combined_for_pairwise)) +
  geom_point(size =3) +
  labs(x = "PC1: 46% variance", y = "PC2: 10% variance", title = "Sample Distance Matrix - LTR vs WTR") +
  theme_minimal()
# pHeatmap and save as .svg and .png - pHM_LTRvsWTR294sigs_pairwise
library(pheatmap)
library(RColorBrewer)
library(ComplexHeatmap)
sampleDists <- dist(t(assay(LTRvsWTRVSD_subset)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- c("L2U1", "L2U2", "L2U3", "L2U4", "L2U5", "L4U1", "L4U2", "L4U3", "L4U4", "L4U5", "W2U1", "W2U2", "W2U3", "W2U4", "W2U5", "W4U1", "W4U2", "W4U3", "W4U4")
colnames(sampleDistMatrix) <- c("tide_runoff_F", "tide_runoff_F", "tide_runoff_F", "tide_runoff_F", "tide_runoff_F", "tide_runoff_F", "tide_runoff_F", "tide_runoff_F", "tide_runoff_F", "tide_runoff_F", "tide_runoff_U", "tide_runoff_U", "tide_runoff_U", "tide_runoff_U", "tide_runoff_U", "tide_runoff_U", "tide_runoff_U", "tide_runoff_U", "tide_runoff_U")
pheatmap(sampleDistMatrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = FALSE, 
         main = "Heatmap LTR vs WTR")


# (3) Long + Tide vs Long + Tide_Runoff [9 sigs]
LTvsLTRvsd <-vst(LongTide_LongTideRunoffdds)
# We just want the significant genes
LTvsLTRsigs<-subset(LongTide_LongTideRunoff_sigs, padj <0.05)
LTvsLTR_sigs<-rownames(LTvsLTRsigs)
VSD.LTvsLTR.subset<- LTvsLTRvsd[rownames(LTvsLTRvsd) %in% LTvsLTR_sigs, ]
summary(VSD.LTvsLTR.subset)
[1] "DESeqTransform object of length 9 with 31 metadata columns"

LTvsLTR_groups_to_keep <- c("long_tide", "long_tiderunoff")
print(LTvsLTR_groups_to_keep)
[1] "long_tide"       "long_tiderunoff"

LTvsLTRsubset_samples <- colData(VSD.LTvsLTR.subset)$combined_for_pairwise %in% LTvsLTR_groups_to_keep
LTvsLTRVSD_subset <- VSD.LTvsLTR.subset[, LTvsLTRsubset_samples]
# plotPCA and save as .svg and .png - PCALTvsLTR9sigs_pairwise
plotPCA(LTvsLTRVSD_subset, intgroup = "combined_for_pairwise")
#ggplot and save as .svg and .png - SDM_LTvsLTR9sigs_pairwise
LTvsLTRgulfkillifish_PCA_t_dim<-plotPCA(LTvsLTRVSD_subset, intgroup = "combined_for_pairwise", returnData = T) #save
ggplot(LTvsLTRgulfkillifish_PCA_t_dim, aes(x = PC1, y = PC2, color = combined_for_pairwise)) +
  geom_point(size =3) +
  labs(x = "PC1: 53% variance", y = "PC2: 23% variance", title = "Sample Distance Matrix - LT vs LTR") +
  theme_minimal()
# pHeatmap and save as .svg and .png - pHM_LTvsLTR9sigs_pairwise
library(pheatmap)
library(RColorBrewer)
library(ComplexHeatmap)
sampleDists <- dist(t(assay(LTvsLTRVSD_subset)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- c("L1C1", "L1C2", "L1C3", "L1C4", "L1C5", "L2U1", "L2U2", "L2U3", "L2U4", "L2U5", "L3C1", "L3C2", "L3C3", "L3C4", "L3C5", "L4U1", "L4U2", "L4U3", "L4U4", "L4U5")
colnames(sampleDistMatrix) <- c("tide_F", "tide_F", "tide_F", "tide_F", "tide_F", "tide_runoff_F", "tide_runoff_F", "tide_runoff_F", "tide_runoff_F", "tide_runoff_F", "tide_F", "tide_F", "tide_F", "tide_F", "tide_F", "tide_runoff_F", "tide_runoff_F", "tide_runoff_F", "tide_runoff_F", "tide_runoff_F")
pheatmap(sampleDistMatrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = FALSE, 
         main = "Heatmap LT vs LTR")


# (4) Weakley + Tide vs Weakley + Tide_Runoff [199 sigs]
WTvsWTRvsd <-vst(WeakleyTide_WeakleyTideRunoffdds)
# We just want the significant genes
WTvsWTRsigs<-subset(WeakleyTide_WeakleyTideRunoff_sigs, padj <0.05)
WTvsWTR_sigs<-rownames(WTvsWTRsigs)
VSD.WTvsWTR.subset<- WTvsWTRvsd[rownames(WTvsWTRvsd) %in% WTvsWTR_sigs, ]
summary(VSD.WTvsWTR.subset)
[1] "DESeqTransform object of length 199 with 31 metadata columns"

WTvsWTR_groups_to_keep <- c("weakley_tide", "weakley_tiderunoff")
print(WTvsWTR_groups_to_keep)
[1] "weakley_tide"       "weakley_tiderunoff"

WTvsWTRsubset_samples <- colData(VSD.WTvsWTR.subset)$combined_for_pairwise %in% WTvsWTR_groups_to_keep
WTvsWTRVSD_subset <- VSD.WTvsWTR.subset[, WTvsWTRsubset_samples]
# plotPCA and save as .svg and .png - PCAWTvsWTR199sigs_pairwise
plotPCA(WTvsWTRVSD_subset, intgroup = "combined_for_pairwise")
#ggplot and save as .svg and .png - SDM_WTvsWTR199sigs_pairwise
WTvsWTRgulfkillifish_PCA_t_dim<-plotPCA(WTvsWTRVSD_subset, intgroup = "combined_for_pairwise", returnData = T) #save
ggplot(WTvsWTRgulfkillifish_PCA_t_dim, aes(x = PC1, y = PC2, color = combined_for_pairwise)) +
  geom_point(size =3) +
  labs(x = "PC1: 45% variance", y = "PC2: 10% variance", title = "Sample Distance Matrix - WT vs WTR") +
  theme_minimal()
# pHeatmap and save as .svg and .png - pHM_WTvsWTR199sigs_pairwise
library(pheatmap)
library(RColorBrewer)
library(ComplexHeatmap)
sampleDists <- dist(t(assay(WTvsWTRVSD_subset)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- c("W1C1", "W1C2", "W1C3", "W1C4", "W2U1", "W2U2", "W2U3", "W2U4", "W2U5", "W3C1", "W3C2", "W3C3", "W3C4", "W4U1", "W4U2", "W4U3", "W4U4")
colnames(sampleDistMatrix) <- c("tide_U", "tide_U", "tide_U", "tide_U", "tide_runoff_U", "tide_runoff_U", "tide_runoff_U", "tide_runoff_U", "tide_runoff_U", "tide_U", "tide_U", "tide_U", "tide_U", "tide_runoff_U", "tide_runoff_U", "tide_runoff_U", "tide_runoff_U")
pheatmap(sampleDistMatrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = FALSE, 
         main = "Heatmap WT vs WTR")



# Proceed to PWCtablemodify.R

#############################################################################################


