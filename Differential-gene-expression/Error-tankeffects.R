# Input files containing gene data
gulfkillifishcountdata<-as.matrix(read.csv("/Users/adc0064/Library/CloudStorage/Box-Box/Bernal_lab/Antrelle/GulfKillifish/gene_count_matrix-Copycopy.csv", row.names="gene_id"))
gulfkillifishcoldatatank<-read.csv("/Users/adc0064/Library/CloudStorage/Box-Box/Bernal_lab/Antrelle/GulfKillifish/sample_info_GKF_tank.csv", header = TRUE, row.names = 1)

# Use the "coldata$___" command to to convert the tank and origin columns to a **factor**.
# These columns are most likely stored as a character vector (text) or numeric vector, depending on how it was
# originally loaded or created.
# **Assigning back**:
# By assigning `as.factor(gulfkillifishcoldata$tank)` back to `gulfkillifishcoldata$tank`, this code overwrites the 
# original `tank` column with its factor-converted version.
# Why Convert to a Factor?
# Factors are especially useful in statistical analyses and modeling in R. Many R packages (including DESeq2) use 
# factors to interpret categorical variables correctly. For example, DESeq2 uses factors to recognize and handle 
# experimental conditions like `tank` groups, which are often needed for differential expression analysis.
# Without converting `tank` to a factor, the program might not recognize it as a categorical variable, leading 
# to potential errors or misinterpretations in the analysis.
gulfkillifishcoldatatank$origin<-as.factor(gulfkillifishcoldatatank$origin)
gulfkillifishcoldatatank$tank<-as.factor(gulfkillifishcoldatatank$tank)

# Check that the columns and rows have the same name
all(rownames(gulfkillifishcoldatatank) == colnames(gulfkillifishcountdata))
[1] TRUE


# Run a likelihood ratio test to determine if the effect of origin:tank, origin, or tank is significant
GFKdds_tank<-DESeqDataSetFromMatrix(countData = gulfkillifishcountdata, colData = gulfkillifishcoldatatank, design =~ tank)
GFKdds_tank<-DESeq(GFKdds_tank)


tank_LRT <- DESeq(GFKdds_tank, test = "LRT", reduced = ~1)
tank_LRT.res <- results(tank_LRT, alpha = 0.05)
summary(tank_LRT.res)
out of 22485 with nonzero total read count
              # adjusted p-value < 0.05
              # LFC > 0 (up)       : 626, 2.8%
              # LFC < 0 (down)     : 996, 4.4%
              # outliers [1]       : 116, 0.52%
              # low counts [2]     : 4291, 19%
              # (mean count < 7)

##This is not helpful because this is too many. This model is only considering "tank" which doesn't help because
##the model will not run with origin + tank or treatment + tank because its saying its a linear combination as in
##tank is already embedded into each main variable origin and treatment.

###Going to run pairwise comparison on just the replicate tanks to see how many genes there are

          ###Metadata:
              ##Tanks 1 and 2 = long bayou + tide
              ##Tanks 3 and 4 = long bayou + tide/runoff
              ##Tanks 5 and 6 = weakley bayou + tide
              ##Tanks 7 and 8 = weakley bayou + tide/runoff


##Pairwise comparison for Tanks 1 and 2 (long bayou + tide)
res_1v2 <- results(GFKdds_tank,
                   contrast = c("tank", "1", "2"),
                   alpha = 0.05)

summary(res_1v2)
              # out of 22485 with nonzero total read count
              # adjusted p-value < 0.05
              # LFC > 0 (up)       : 198, 0.88%
              # LFC < 0 (down)     : 97, 0.43%
              # outliers [1]       : 116, 0.52%
              # low counts [2]     : 6009, 27%
              # (mean count < 22)

res_1v2_sig <- subset(res_1v2, padj < 0.05)


##Pairwise comparison for Tanks 3 and 4 (long bayou + tide/runoff)
res_3v4 <- results(GFKdds_tank,
                   contrast = c("tank", "3", "4"),
                   alpha = 0.05)

summary(res_3v4)
              # out of 22485 with nonzero total read count
              # adjusted p-value < 0.05
              # LFC > 0 (up)       : 11, 0.049%
              # LFC < 0 (down)     : 14, 0.062%
              # outliers [1]       : 116, 0.52%
              # low counts [2]     : 433, 1.9%
              # (mean count < 0)

res_3v4_sig <- subset(res_3v4, padj < 0.05)

##Pairwise comparison for Tanks 5 and 6 (weakley bayou + tide)
res_5v6 <- results(GFKdds_tank,
                   contrast = c("tank", "5", "6"),
                   alpha = 0.05)

summary(res_5v6)
              # out of 22485 with nonzero total read count
              # adjusted p-value < 0.05
              # LFC > 0 (up)       : 28, 0.12%
              # LFC < 0 (down)     : 34, 0.15%
              # outliers [1]       : 116, 0.52%
              # low counts [2]     : 433, 1.9%
(mean count < 0)
res_5v6_sig <- subset(res_5v6, padj < 0.05)


##Pairwise comparison for Tanks 7 and 8 (weakley bayou + tide/runoff)
res_7v8 <- results(GFKdds_tank,
                   contrast = c("tank", "7", "8"),
                   alpha = 0.05)

summary(res_7v8)
              # out of 22485 with nonzero total read count
              # adjusted p-value < 0.05
              # LFC > 0 (up)       : 34, 0.15%
              # LFC < 0 (down)     : 24, 0.11%
              # outliers [1]       : 116, 0.52%
              # low counts [2]     : 2175, 9.7%
              # (mean count < 1)
res_7v8_sig <- subset(res_7v8, padj < 0.05)


DESeq2::plotMA(tank_LRT.res, main = "Log2-Fold Change vs. Base Mean", ylim=c(-5,5))

          # Error in checkFullRank(modelMatrix) : 
          #   the model matrix is not full rank, so the model cannot be fit as specified.
          # One or more variables or interaction terms in the design formula are linear
          # combinations of the others and must be removed.
          # 
          # Please read the vignette section 'Model matrix not full rank':
          #   
          #   vignette('DESeq2')






