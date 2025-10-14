################# WGCNA  - RNA-seq #####################
######### Gulf Killifish - Fundulus grandis - Gill ########

# Antrelle D. Clark - Auburn University Department of Biological Sciences - Fall 2025
# Gratitude to Ally Swank, MS and Katie Eaton, PhD for the WGCNA script
# R-script version 4.4.1 (2024-06-14)

setwd("~/Library/CloudStorage/Box-Box/Bernal_lab/Antrelle/GulfKillifish")

library(DESeq2)
BiocManager::install("GO.db")
BiocManager::install("WGCNA")
BiocManager::install("impute")
BiocManager::install("preprocessCore")
install.packages("WGCNA")
library('WGCNA')
#WGCNA v1.73
library(dynamicTreeCut)
library(fastcluster)
library(flashClust)
library(Hmisc)
library(matrixStats)
library(ggforce)
library(dplyr)
library(tidyr)
library(ggplot2)
library(foreach)
library(doParallel)
library(ComplexHeatmap)
library(circlize)
library(splines, lib.loc = "/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library")
library(survival, lib.loc = "/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library")


options(stringAsFactors = FALSE)

file.exists("gene_count_matrix.csv")
#[1] TRUE

x<-read.csv("/Users/adc0064/Library/CloudStorage/Box-Box/Bernal_lab/Antrelle/GulfKillifish/gene_count_matrix.csv", row.names=1, header=TRUE, stringsAsFactors=FALSE )
x <- x[,colnames(x)!="X"]

colnames(x)
#[1] "L1C1" "L1C2" "L1C3" "L1C4" "L1C5" "L2U1" "L2U2" "L2U3" "L2U4" "L2U5" "L3C1" "L3C2" "L3C3" "L3C4" "L3C5" "L4U1"
#[17] "L4U2" "L4U3" "L4U4" "L4U5" "W1C1" "W1C2" "W1C3" "W1C4" "W2U1" "W2U2" "W2U3" "W2U4" "W2U5" "W3C1" "W3C2" "W3C3"
#[33] "W3C4" "W4U1" "W4U2" "W4U3" "W4U4"

head(x)
#L1C1 L1C2 L1C3 L1C4 L1C5 L2U1 L2U2 L2U3 L2U4 L2U5 L3C1 L3C2 L3C3 L3C4 L3C5 L4U1 L4U2 L4U3
#gene:ENSFHEG00000014326  531  368  667  558  415  470  616  429  383  551  552  903  477  542  519  329  482  458
#gene:ENSFHEG00000000334    0    0    0    4    2    0    0    0    0    0    2    0    0    0    0    0    0    0
#gene:ENSFHEG00000013847    0    4    0    0    0    0    0    0    0    0    0    0    0    0    0    0    3    0
#gene:ENSFHEG00000013795    6    3   10    7   13    8    3    2    5   18   12    7    3    6    0    6    9    6
#gene:ENSFHEG00000013942  853  692  740  739  851  708  821  674 1269 1312  779  955  659  762  479  502  809  473
#gene:ENSFHEG00000013756   52   45   55   72   52   24   52   35   54   50   30   75   43   34   64   54   29   53
#L4U4 L4U5 W1C1 W1C2 W1C3 W1C4 W2U1 W2U2 W2U3 W2U4 W2U5 W3C1 W3C2 W3C3 W3C4 W4U1 W4U2 W4U3
#gene:ENSFHEG00000014326  586  315  413  589  539  395  680  540  806  676  584  362  691  326  373  648  497  410
#gene:ENSFHEG00000000334    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    2    0    0
#gene:ENSFHEG00000013847    0    0    0    0    0    0    0    0    0    0    3    0    4    0    0    0    0    0
#gene:ENSFHEG00000013795    8    1    6    6    4   10   12    5    8   10    8   10    2    0    5    3    7    5
#gene:ENSFHEG00000013942  868  913  636  780  609  887  784  831  409  737  676  596  626  519  651  870  509  711
#gene:ENSFHEG00000013756   54   30   75   32   47   83   60  100   43   97   69   41   76   48   50   63   26   27
#W4U4
#gene:ENSFHEG00000014326  662
#gene:ENSFHEG00000000334    0
#gene:ENSFHEG00000013847    0
#gene:ENSFHEG00000013795    6
#gene:ENSFHEG00000013942  695
#gene:ENSFHEG00000013756   61

nrow(x)
#[1] 23469


#filter low counts
x$filtr<-apply(x, 1, function(k) mean(k > 10)) > 0.5
x<-x[!(x$filtr=="FALSE"),]
nrow(x)
#[1] 17108

#This is the code to only analyze the 50% with the highest variation
#counts$variance = apply(counts, 1, var)
#counts2 = counts[counts$variance >= quantile(counts$variance, c(.50)), ] #50% most variable genes
#counts2$variance <- NULL
#counts=counts2

x$filtr=NULL
#outlier removal if necessary

metaData <- read.csv('f.grandis_attributes.csv', header = TRUE, sep = ",") 
head(metaData)
nrow(metaData)
#sample_ID   treatment       origin combined_for_pairwise control.salinity treatment.salinity length..cm.
#1      L1C1        tide long_bayou_F             long_tide               15                 10         6.2
#2      L1C2        tide long_bayou_F             long_tide               15                 10         6.7
#3      L1C3        tide long_bayou_F             long_tide               15                 10        10.6
#4      L1C4        tide long_bayou_F             long_tide               15                 10         8.5
#5      L1C5        tide long_bayou_F             long_tide               15                 10         8.2
#6      L2U1 tide_runoff long_bayou_F       long_tiderunoff               15                  5         8.5
#weight..g. sex
#1        6.0 unk
#2        5.8 unk
#3       20.8   M
#4       10.1 unk
#5        9.3 unk
#6       11.2 unk

totalCounts=colSums(x) 
min(totalCounts)
#[1] 14684991
max(totalCounts)
#[1] 23312236
mean(totalCounts) 
#[1] 17490210


dds <- DESeqDataSetFromMatrix(
  countData = x, # Our prepped data frame with counts
  colData = metaData, # Data frame with annotation for our samples
  design = ~1 # Here we are not specifying a model
)
dds<-dds[rowSums(counts(dds))>0,]
summary(dds)
#[1] "DESeqDataSet object of length 17108 with 0 metadata columns"


#to obtain variance stabilized data follow: 
head(dds)
#class: DESeqDataSet 
#dim: 6 37 
#metadata(1): version
#assays(1): counts
#rownames(6): gene:ENSFHEG00000014326 gene:ENSFHEG00000013942 ... gene:ENSFHEG00000013549
#gene:ENSFHEG00000013406
#rowData names(0):
#colnames(37): L1C1 L1C2 ... W4U3 W4U4
#colData names(9): sample_ID treatment ... weight..g. sex


vsd=vst(dds)
head(vsd)


#########START HERE!!!
# Check that the data has the correct format for many functions operating on multiple sets:
#WGCNA THE DATA NEEDS TO BE TRANSPOSED to make dendrogram

gsg = goodSamplesGenes(x, verbose = 3);
gsg$allOK
#[1] TRUE

xx <- assay(vsd) %>%
  t() # Transpose this data

dim(vsd) 
#[1] 17108    37
dim(xx) 
#[1]    37 17108
rownames(xx) #make sure the rownames are actually 
#[1] "L1C1" "L1C2" "L1C3" "L1C4" "L1C5" "L2U1" "L2U2" "L2U3" "L2U4" "L2U5" "L3C1" "L3C2" "L3C3" "L3C4" "L3C5" "L4U1"
#[17] "L4U2" "L4U3" "L4U4" "L4U5" "W1C1" "W1C2" "W1C3" "W1C4" "W2U1" "W2U2" "W2U3" "W2U4" "W2U5" "W3C1" "W3C2" "W3C3"
#[33] "W3C4" "W4U1" "W4U2" "W4U3" "W4U4"


sampleTree = hclust(dist(xx), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
#plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")

######Potential Outliers but not removed: L1C2, L4U1, and W2U2

#saved as WGCNA_detectoutliers.pdf


#trait table
# population: 0 = LB 1 = WB
# treatment: 0 = control 1 = stress
# salinity: 0 = 10 1 = 5

traitData = read.csv("f.grandis_traits.csv",header=TRUE);
dim(traitData)
#[1] 37  4
names(traitData)
#[1] "sample_ID"  "population" "treatment"  "salinity"


sample_ID = rownames(xx);
traitRows = match(sample_ID, traitData$sample_ID);
datTraits = traitData[traitRows, -1];
rownames(datTraits) = traitData[traitRows, 1];
collectGarbage();

#store the data in an R function
save(xx, datTraits, file = "Fgrand-WGCNA-step1.RData")


##To Find blocks of modules, use the data that is not transposed 
# Choose a set of soft-thresholding powers
lnames= load(file="Fgrand-WGCNA-step1.RData")
lnames
#[1] "xx"        "datTraits"
exprSize = checkSets(xx, checkStructure = TRUE);
nSets = exprSize$nSets;
nSets = checkSets(xx, checkStructure = TRUE)$nSets

powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(xx, powerVector = powers, verbose = 6)
#pickSoftThreshold: will use block size 2615
# Plot the results:
sizeGrWindow(9, 5) #this will return NULL because you are not producing a plot yet
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


###Look at the "Scale Independence" and "Mean Connectivity" to see where the curves start to stabilize
###Based on just "Scale Independence": Typically you want to go with the lowest power that reaches ~0.9 
###and is stable but mine takes a dip and goes back up and reaches 0.9 at soft threshold power of 10 
###and then continues to increase. Could choose 4-7 or 10, but going with 10 to be safe and avoid dips.
###Based on both (which is why it is important and helps make the decision clear): both curves stabilize
###at the soft threshold (power) of 10 so going with 10 for downstream analyses
#build networks
net = blockwiseModules(xx, power = 10, corType="bicor", networkType = "signed",
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = FALSE,
                       deepSplit = 2,
                       pamRespectsDendro = F,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "Fgrand-net",
                       verbose = 5, maxBlockSize = 5000)
names(net)
#[1] "colors"         "unmergedColors" "MEs"            "goodSamples"    "goodGenes"      "dendrograms"   
#[7] "TOMFiles"       "blockGenes"     "blocks"         "MEsOK


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
consMEs = net$MEs;
consTree = net$dendrograms[[1]]

sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

##Graph is titled Cluster Dendrogram

save(consMEs, moduleLabels, moduleColors, consTree, file = "Fgrand-NetworkConstruction-auto.RData")

############################

#Clusters to All Samples (not very clear pattern)
module_df <- data.frame(gene_id = names(net$colors),
                        colors = labels2colors(net$colors))
module_df
module_df[1:5,]
write.table(module_df, file ="Fgrand-gene_modules.txt",sep = "\t", quote=FALSE)

write.csv(module_df, file="gene_modules.csv", quote = FALSE)


MEs0 <- moduleEigengenes(xx, mergedColors)$eigengenes
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)
MEs0$treatment = row.names(MEs0)


mME = MEs0 %>% pivot_longer(-treatment) %>% mutate(
  name = gsub("ME", "", name),
  name = factor(name, levels = module_order))

#-treatment take all columns except treatment and pivot them into long format

mME %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")


############################

#Cluster Dendrogram Dissimilarity and Merged Dynamic

softPower = 10; #10 with outliers included
adjacency = adjacency(xx, power = softPower, type='signed');
TOM = TOMsimilarity(adjacency, TOMType='signed');
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
#sizeGrWindow(12,9)
#pdf(file=Dendrogram_signed_BM10.pdf, width=20, height=20)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
#dev.off()
minModuleSize = 30;

# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
#dynamicMods


# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
#dynamicColors

# Plot the dendrogram and colors underneath
#sizeGrWindow(8,6)
#pdf(file=Dendrogram_signed_BM10_colors.pdf, width=20, height=20)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
#dev.off()

# Calculate eigengenes
MEList = moduleEigengenes(xx, colors = dynamicColors)
MEList$eigengenes #gives you the eigenes by sample 
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
#pdf(file=ClusteringEigengenes.pdf, width=20, height=20)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
#dev.off()



MEDissThres = 0.25 
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(xx, dynamicColors, cutHeight = MEDissThres, verbose = 3)

# The merged module colors
mergedColors = merge$colors; #Change "net" to "merge" to do the rest with the merged colors!!!!
table(mergedColors)
mergedColors
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

sizeGrWindow(12, 9)
#pdf(file = "DendroAndColors_sft6_bm10.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off() #this is to save files

# MERGING: Rename to moduleColors that are similar
moduleColors = mergedColors
table(moduleColors)
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "Fgrand-RNA_WGCNA_networkConstruct_signed.RData")
# Define numbers of genes and samples
nGenes = ncol(xx);
nSamples = nrow(xx);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(xx, moduleColors)$eigengenes 
MEs = orderMEs(MEs0)


######################
#correlations of traits with eigengenes

moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
Colors=sub("ME","",names(MEs))


## Will display their p-values
#textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "");
textMatrix =  paste(signif(moduleTraitPvalue, 1))
dim(textMatrix) = dim(moduleTraitCor)
par(mfrow = c(1,1))
par(mar=c(1,1,1,1))
# Display the correlation values within a heatmap plot
pdf(file="trait-cor-heatmap.pdf", width=7, height=7)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               xLabelsAngle = 0,
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = TRUE,
               cex.text = 0.75,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

####for this heatmap, box color = correlation and number in box is P-value.
###You can edit the text to include both if you want where the correlation number will be shown
###followed by a number in parentheses which would be the P-value. Here, we only plotted the
###P-value inside the boxes with the color of the boxes being associated with correlation

               
###AND/OR###

## Will display their correlation values
textMatrix =  paste(signif(moduleTraitCor, 2))
#textMatrix =  paste(signif(moduleTraitPvalue, 1))
dim(textMatrix) = dim(moduleTraitCor)
par(mfrow = c(1,1))
par(mar=c(1,1,1,1))
# Display the correlation values within a heatmap plot
pdf(file="trait-cor-heatmap-included.pdf", width=7, height=7)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               xLabelsAngle = 0,
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = TRUE,
               cex.text = 0.75,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()
               


MEDiss= 1-cor(MEs)
METree= flashClust(as.dist(MEDiss), method= "average")
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")


#################################################################################
### make heatmap of specific modules of interest
head(MEs)

#Which modules have biggest differences across treatment groups?
#First we double check that our samples are still in order.
all.equal(metaData$sample_ID, rownames(MEs))

# Create the design matrix from the `treatment` variable
des_mat <- model.matrix(~ metaData$treatment)

#Run linear model on each module. 
#Limma wants our tests to be per row, so we also need to transpose so the eigengenes are rows
# lmFit() needs a transposed version of the matrix
fit <- limma::lmFit(t(MEs), design = des_mat)
# Apply empirical Bayes to smooth standard errors
fit <- limma::eBayes(fit)

# Apply multiple testing correction and obtain stats
stats_df <- limma::topTable(fit, number = ncol(MEs)) %>%
  tibble::rownames_to_column("module")

head(stats_df)
#module       logFC       AveExpr         t     P.Value  adj.P.Val         B
#1   MEdarkgreen  0.15605197  1.522806e-16  2.874688 0.004146416 0.09951399 -2.048523
#2     MEdarkred  0.12403071  1.181487e-16  2.284813 0.022572323 0.27086788 -3.442233
#3      MEsalmon -0.10472049 -3.375678e-17 -1.929093 0.054055927 0.28759867 -4.128666
#4   MEsteelblue -0.10376058  1.912884e-16 -1.911410 0.056291965 0.28759867 -4.159765
#5 MElightyellow -0.09970102 -8.326673e-17 -1.836627 0.066618123 0.28759867 -4.288119
#6      MEgrey60  0.09510429 -9.639436e-17  1.751950 0.080147376 0.28759867 -4.427274



#retrieve gene list for the modules, can be modified in excel to run GO analyses
annot=read.table(file="~/Library/CloudStorage/Box-Box/Bernal_lab/Antrelle/GulfKillifish/Mummichog_aggregated_GO_terms_filtered.tab", header=T, sep='\t')
probes = colnames(xx)
head(probes)
#[1] NA NA NA NA NA NA

# Remove the "gene:" prefix
probes_clean <- sub("^gene:", "", probes)
#Match
probes2annot <- match(probes_clean, annot$V1)  # or annot$gene if your column is named 'gene'

head(probes2annot)
#[1] 12742 12402 12232 12248 12043 11918


datGS.Traits=data.frame(cor(xx,datTraits,use="p"))
names(datGS.Traits)=paste("cor",names(datGS.Traits),sep=".")
datME=moduleEigengenes(xx,moduleColors)$eigengenes
datKME=signedKME(xx, datME, outputColumnName="MM.")
datOutput=list(gene=colnames(xx),colors = moduleColors)
datOutput=as.vector(datOutput)
datOutput$gene = as.character(datOutput$gene)

datOutput0 = as.data.frame(datOutput)

#As a sanity check, let's use ggplot to see what module XX's eigengene looks like between treatment groups.
#First we need to set up the module eigengene for this module with the sample metadata labels we need.
module_darkolivegreen_df <- MEs %>%
  tibble::rownames_to_column("FishName") %>%
  # Here we are performing an inner join with a subset of metadata
  dplyr::inner_join(metaData %>%
                      dplyr::select(sample_ID, treatment),
                    by = c("FishName" = "sample_ID")
  )


ggplot(
  module_darkolivegreen_df,
  aes(
    x = treatment,
    y = MEdarkolivegreen,
    color = treatment
  )
) +
  # a boxplot with outlier points hidden (they will be in the sina plot)
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  # A single plot to show all of the individual data points
  ggforce::geom_sina(maxwidth = 0.3) +
  theme_classic()

#If you want to know which of your genes make up a modules, you can look at the $colors slot.
gene_module_key <- as.data.frame(datOutput) %>%
  # Let's add the `ME` part so its more clear what these numbers are and it matches elsewhere
  dplyr::mutate(colors = paste0("ME", colors))
gene_module_key <- gene_module_key %>%
  dplyr::filter(colors == "MEdarkolivegreen")

head(gene_module_key)
gene_module_key
#gene           colors
#1  gene:ENSFHEG00000006376 MEdarkolivegreen
#2  gene:ENSFHEG00000005960 MEdarkolivegreen
#3  gene:ENSFHEG00000019450 MEdarkolivegreen
#4  gene:ENSFHEG00000016616 MEdarkolivegreen
#5  gene:ENSFHEG00000015482 MEdarkolivegreen
#6  gene:ENSFHEG00000008330 MEdarkolivegreen
#7  gene:ENSFHEG00000015138 MEdarkolivegreen
#8  gene:ENSFHEG00000016678 MEdarkolivegreen
#9  gene:ENSFHEG00000020935 MEdarkolivegreen
#10 gene:ENSFHEG00000003380 MEdarkolivegreen
#11 gene:ENSFHEG00000018346 MEdarkolivegreen
#12 gene:ENSFHEG00000009096 MEdarkolivegreen
#13 gene:ENSFHEG00000002951 MEdarkolivegreen
#14 gene:ENSFHEG00000001966 MEdarkolivegreen
#15 gene:ENSFHEG00000016045 MEdarkolivegreen
#16 gene:ENSFHEG00000018638 MEdarkolivegreen
#17 gene:ENSFHEG00000018922 MEdarkolivegreen
#18 gene:ENSFHEG00000001609 MEdarkolivegreen
#19 gene:ENSFHEG00000023006 MEdarkolivegreen
#20 gene:ENSFHEG00000018653 MEdarkolivegreen
#21 gene:ENSFHEG00000012270 MEdarkolivegreen
#22 gene:ENSFHEG00000022916 MEdarkolivegreen
#23 gene:ENSFHEG00000001789 MEdarkolivegreen
#24 gene:ENSFHEG00000020575 MEdarkolivegreen
#25 gene:ENSFHEG00000009627 MEdarkolivegreen
#26 gene:ENSFHEG00000003206 MEdarkolivegreen
#27 gene:ENSFHEG00000017446 MEdarkolivegreen
#28 gene:ENSFHEG00000023023 MEdarkolivegreen
#29 gene:ENSFHEG00000007017 MEdarkolivegreen
#30 gene:ENSFHEG00000016153 MEdarkolivegreen
#31 gene:ENSFHEG00000022139 MEdarkolivegreen
#32 gene:ENSFHEG00000002494 MEdarkolivegreen
#33 gene:ENSFHEG00000005258 MEdarkolivegreen
#34 gene:ENSFHEG00000015620 MEdarkolivegreen
#35 gene:ENSFHEG00000022566 MEdarkolivegreen
#36 gene:ENSFHEG00000017925 MEdarkolivegreen


#Save the gene to module key to a TSV file
#this list was used to identify annotated genes on the command line
readr::write_csv(gene_module_key, "gill_wgcna_module_darkolivegreen.csv")

write.csv(net$colors, "module-genes_gill_WGCNA.csv")

#############################################################################


#Make a custom heatmap function for specific modules

make_module_heatmap <- function(module_name,
                                expression_mat = xx,
                                metadata_df = metaData,
                                gene_module_key_df = gene_module_key,
                                module_eigengenes_df = module_eigengenes) {
  
  
  # Create a summary heatmap of a given module.
  
  #
  # Args:
  # module_name: a character indicating what module should be plotted, e.g. "ME0"
  # expression_mat: The full gene expression matrix. Default is `normalized_counts`.
  # metadata_df: a data frame with Treatment and time_point
  #              as columns. Default is `metadata`.
  # gene_module_key: a data.frame indicating what genes are a part of what modules. Default is `gene_module_key`.
  # module_eigengenes: a sample x eigengene data.frame with samples as row names. Default is `module_eigengenes`.
  #
  # Returns:
  # A heatmap of expression matrix for a module's genes, with a barplot of the
  # eigengene expression for that module.
  
  # Set up the module eigengene with its Treatment
  #module_eigengene <- module_eigengenes_df %>%
  #dplyr::select(all_of(module_name)) %>%
  #tibble::rownames_to_column("FishID")
  
  
  # Set up column annotation from metadata
  col_annot_df <- metadata_df %>%
    # Only select the combined for pairwise and sample ID columns
    dplyr::select(combined_for_pairwise, sample_ID) %>%
    # Add on the eigengene expression by joining with sample IDs
    #dplyr::inner_join(module_eigengene, by = "sample_ID") %>%
    # Arrange by patient and time point
    dplyr::arrange(combined_for_pairwise, sample_ID) %>%
    # Store sample
    tibble::column_to_rownames("sample_ID")
  
  
  # Create the ComplexHeatmap column annotation object
  col_annot <- ComplexHeatmap::HeatmapAnnotation(
    # Supply treatment labels
    combined_for_pairwise = col_annot_df$combined_for_pairwise,
    # Add annotation barplot
    #module_eigengene = ComplexHeatmap::anno_barplot(dplyr::select(col_annot_df, module_name)),
    # Pick colors for each experimental group in treatment
    col = list(combined_for_pairwise = c("long_tide" = "#FF0000", "long_tiderunoff" = "#009E73", "weakley_tide" = "#56B4E9", "weakley_tiderunoff" = "#800080"))
  )
  
  
  # Get a vector of the Ensembl gene IDs that correspond to this module
  module_genes <- gene_module_key_df %>%
    dplyr::filter(colors == module_name) %>%
    dplyr::pull(gene)
  
  # Set up the gene expression data frame
  mod_mat <- expression_mat %>%
    t() %>%
    as.data.frame() %>%
    # Only keep genes from this module
    dplyr::filter(rownames(.) %in% module_genes) %>%
    # Order the samples to match col_annot_df
    dplyr::select(rownames(col_annot_df)) %>%
    # Data needs to be a matrix
    as.matrix()
  
  
  # Normalize the gene expression values
  mod_mat <- mod_mat %>%
    # Scale can work on matrices, but it does it by column so we will need to
    # transpose first
    t() %>%
    scale() %>%
    # And now we need to transpose back
    t()
  
  
  # Create a color function based on standardized scale
  color_func <- circlize::colorRamp2(
    c(-2, 0, 2),
    c("blue", "white", "red"))
  
  
  # Plot on a heatmap
  heatmap <- ComplexHeatmap::Heatmap(mod_mat,
                                     name = module_name,
                                     # Supply color function
                                     col = color_func,
                                     # Supply column annotation
                                     bottom_annotation = col_annot,
                                     # We don't want to cluster samples
                                     cluster_columns = FALSE,
                                     # We don't need to show sample or gene labels
                                     show_row_names = FALSE,
                                     show_column_names = FALSE
  )
  
  # Return heatmap
  return(heatmap)
} 



module_darkgreen_heatmap <- make_module_heatmap(module_name = "MEdarkgreen")
ComplexHeatmap::draw(module_darkgreen_heatmap)
# Print out the plot
module_darkgreen_heatmap

#Save plot to png
pdf("darkgreen_heatmap.pdf")
module_darkgreen_heatmap
dev.off()



# Make the Heatmap object
module_darkgreen_heatmap <- make_module_heatmap(module_name = "MEdarkgreen")

# Plot it on screen
ComplexHeatmap::draw(module_darkgreen_heatmap)

# Save to file
pdf("darkgreen_heatmap.pdf")
ComplexHeatmap::draw(module_darkgreen_heatmap)
dev.off()


head(annot)


#As a sanity check, let's use ggplot to see what module XX's eigengene looks like between treatment groups.
#First we need to set up the module darkgreen eigengene for this module with the sample metadata labels we need.
module_darkgreen_df <- MEs %>%
  tibble::rownames_to_column("FishName") %>%
  # Here we are performing an inner join with a subset of metadata
  dplyr::inner_join(metaData %>%
                      dplyr::select(sample_ID, treatment),
                    by = c("FishName" = "sample_ID")
  )


ggplot(
  module_darkgreen_df,
  aes(
    x = treatment,
    y = MEdarkgreen,
    color = treatment
  )
) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  ggforce::geom_sina(maxwidth = 0.3) +
  scale_color_manual(
    values = c("tide" = "blue", "tide_runoff" = "red")
  ) +
  theme_classic()

#If you want to know which of your genes make up a modules, you can look at the $colors slot.
gene_module_key <- as.data.frame(datOutput) %>%
  # Let's add the `ME` part so its more clear what these numbers are and it matches elsewhere
  dplyr::mutate(colors = paste0("ME", colors))
gene_module_key <- gene_module_key %>%
  dplyr::filter(colors == "MEdarkgreen")

head(gene_module_key)
gene_module_key


table(moduleColors) #this shows you all of the modules and the number of genes for each
#moduleColors
#blue           cyan      darkgreen       darkgrey darkolivegreen        darkred  darkturquoise 
#2156           1597            806            301             36            201            318 
#green    greenyellow         grey60     lightgreen    lightyellow   midnightblue  paleturquoise 
#1235            484           1908           1054            220            296             76 
#pink         purple            red    saddlebrown         salmon        skyblue      steelblue 
#626            488           1406            118            353            124             87 
#tan      turquoise          white 
#412           2668            138 


#Save the gene to module key to a TSV file
#this list was used to identify annotated genes on the command line
readr::write_csv(gene_module_key, "gill_wgcna_module_darkgreen.csv")



#extract the genes from the MEdarkgreen for heatmap

# Get the genes in darkgreen module
# Example: "darkgreen" module
genes_MEdarkgreen <- names(xx)[moduleColors == "MEdarkgreen"]


#Subset the data
expr_MEdarkgreen <- xx[, genes_MEdarkgreen]

# Scale by gene (row-wise scaling)
expr_MEdarkgreen_scaled <- t(scale(t(expr_MEdarkgreen)))


annotation_col <- data.frame(Groups = metaData$combined_for_pairwise)
rownames(annotation_col) <- rownames(metaData)   # sample IDs must match colnames(xx)

# Optional custom colors for treatments
ann_colors <- list(
  Groups = c("long_tide" = "#FF0000", 
             "long_tiderunoff" = "#009E73", 
             "weakley_tide" = "#56B4E9", 
             "weakley_tiderunoff" = "#800080")
)



#Plot heatmap with samples grouped by treatment
pheatmap::pheatmap(expr_MEdarkgreen_scaled,
                   show_rownames = FALSE,
                   annotation_row = annotation_col,  # annotation along rows since samples are rows
                   annotation_colors = ann_colors,
                   cluster_cols = TRUE,   # cluster genes
                   cluster_rows = FALSE)  # don’t cluster samples


################################################################################

#######THIS IS FOR CREATING THE MASTER GENE MODULE FILE NEEDED FOR THE GO ANALYSIS DOWNSTREAM
# Assume gene_module_key has two columns: gene, colors


geneIDs <- colnames(xx)   # assuming your expression matrix has genes as columns
gene_module_key <- data.frame(
  gene = geneIDs,
  colors = moduleColors
)

# Get all module names
all_modules <- sort(unique(gene_module_key$colors))

# Build the binary membership table
GMM <- gene_module_key %>%
  dplyr::mutate(dummy = 1) %>% 
  tidyr::pivot_wider(
    id_cols = gene,            # use only gene ID as the row identifier
    names_from = colors,       # make one column per module color
    values_from = dummy,
    values_fill = list(dummy = 0)
  )

# Add moduleColor column (original assignment)
GMM <- dplyr::left_join(
  gene_module_key,
  GMM,
  by = "gene"
)

# Reorder columns to match your example
GMM <- GMM %>%
  dplyr::rename(gene_ID = gene, moduleColor = colors) %>%
  dplyr::select(gene_ID, moduleColor, everything())

# Save to CSV for GO-MWU
write.csv(GMM, file = "gfk-master-gene-module-membership-binary.csv", quote = FALSE, row.names = FALSE)

#this table had "gene:" in front of the gene ensemblID so I went in on excel and replaced "gene:" with ____ (leave blank) to remove them all 



##########Boxplots for all significant modules (8 total)

#####MODULE SALMON
#As a sanity check, let's use ggplot to see what module salmon eigengene looks like between treatment groups.
#First we need to set up the module eigengene for this module with the sample metadata labels we need.
module_salmon_df <- MEs %>%
  tibble::rownames_to_column("FishName") %>%
  # Here we are performing an inner join with a subset of metadata
  dplyr::inner_join(metaData %>%
                      dplyr::select(sample_ID, combined_for_pairwise),
                    by = c("FishName" = "sample_ID")
  )


ggplot(
  module_salmon_df,
  aes(
    x = combined_for_pairwise,
    y = MEsalmon,
    color = combined_for_pairwise
  )
) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  ggforce::geom_sina(maxwidth = 0.3) +
  scale_color_manual(
    values = c(
      "long_tide" = "#F9776E",
      "long_tiderunoff" = "#7CAF04",
      "weakley_tide" = "#0BBFC4",
      "weakley_tiderunoff" = "#A020F0"
    )
  ) +
  theme_classic()


#####MODULE RED
#As a sanity check, let's use ggplot to see what module red eigengene looks like between treatment groups.
#First we need to set up the module eigengene for this module with the sample metadata labels we need.
module_red_df <- MEs %>%
  tibble::rownames_to_column("FishName") %>%
  # Here we are performing an inner join with a subset of metadata
  dplyr::inner_join(metaData %>%
                      dplyr::select(sample_ID, combined_for_pairwise),
                    by = c("FishName" = "sample_ID")
  )


ggplot(
  module_red_df,
  aes(
    x = combined_for_pairwise,
    y = MEred,
    color = combined_for_pairwise
  )
) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  ggforce::geom_sina(maxwidth = 0.3) +
  scale_color_manual(
    values = c(
      "long_tide" = "#F9776E",
      "long_tiderunoff" = "#7CAF04",
      "weakley_tide" = "#0BBFC4",
      "weakley_tiderunoff" = "#A020F0"
    )
  ) +
  theme_classic()


#####MODULE BLUE
#As a sanity check, let's use ggplot to see what module blue eigengene looks like between treatment groups.
#First we need to set up the module eigengene for this module with the sample metadata labels we need.
module_blue_df <- MEs %>%
  tibble::rownames_to_column("FishName") %>%
  # Here we are performing an inner join with a subset of metadata
  dplyr::inner_join(metaData %>%
                      dplyr::select(sample_ID, combined_for_pairwise),
                    by = c("FishName" = "sample_ID")
  )


ggplot(
  module_blue_df,
  aes(
    x = combined_for_pairwise,
    y = MEblue,
    color = combined_for_pairwise
  )
) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  ggforce::geom_sina(maxwidth = 0.3) +
  scale_color_manual(
    values = c(
      "long_tide" = "#F9776E",
      "long_tiderunoff" = "#7CAF04",
      "weakley_tide" = "#0BBFC4",
      "weakley_tiderunoff" = "#A020F0"
    )
  ) +
  theme_classic()


#####MODULE PINK
#As a sanity check, let's use ggplot to see what module pink eigengene looks like between treatment groups.
#First we need to set up the module eigengene for this module with the sample metadata labels we need.
module_pink_df <- MEs %>%
  tibble::rownames_to_column("FishName") %>%
  # Here we are performing an inner join with a subset of metadata
  dplyr::inner_join(metaData %>%
                      dplyr::select(sample_ID, combined_for_pairwise),
                    by = c("FishName" = "sample_ID")
  )


ggplot(
  module_pink_df,
  aes(
    x = combined_for_pairwise,
    y = MEpink,
    color = combined_for_pairwise
  )
) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  ggforce::geom_sina(maxwidth = 0.3) +
  scale_color_manual(
    values = c(
      "long_tide" = "#F9776E",
      "long_tiderunoff" = "#7CAF04",
      "weakley_tide" = "#0BBFC4",
      "weakley_tiderunoff" = "#A020F0"
    )
  ) +
  theme_classic()


#####MODULE WHITE
#As a sanity check, let's use ggplot to see what module white eigengene looks like between treatment groups.
#First we need to set up the module eigengene for this module with the sample metadata labels we need.
module_white_df <- MEs %>%
  tibble::rownames_to_column("FishName") %>%
  # Here we are performing an inner join with a subset of metadata
  dplyr::inner_join(metaData %>%
                      dplyr::select(sample_ID, combined_for_pairwise),
                    by = c("FishName" = "sample_ID")
  )


ggplot(
  module_white_df,
  aes(
    x = combined_for_pairwise,
    y = MEwhite,
    color = combined_for_pairwise
  )
) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  ggforce::geom_sina(maxwidth = 0.3) +
  scale_color_manual(
    values = c(
      "long_tide" = "#F9776E",
      "long_tiderunoff" = "#7CAF04",
      "weakley_tide" = "#0BBFC4",
      "weakley_tiderunoff" = "#A020F0"
    )
  ) +
  theme_classic()


#####MODULE TURQUOISE
#As a sanity check, let's use ggplot to see what module turquoise eigengene looks like between treatment groups.
#First we need to set up the module eigengene for this module with the sample metadata labels we need.
module_turquoise_df <- MEs %>%
  tibble::rownames_to_column("FishName") %>%
  # Here we are performing an inner join with a subset of metadata
  dplyr::inner_join(metaData %>%
                      dplyr::select(sample_ID, combined_for_pairwise),
                    by = c("FishName" = "sample_ID")
  )


ggplot(
  module_turquoise_df,
  aes(
    x = combined_for_pairwise,
    y = MEturquoise,
    color = combined_for_pairwise
  )
) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  ggforce::geom_sina(maxwidth = 0.3) +
  scale_color_manual(
    values = c(
      "long_tide" = "#F9776E",
      "long_tiderunoff" = "#7CAF04",
      "weakley_tide" = "#0BBFC4",
      "weakley_tiderunoff" = "#A020F0"
    )
  ) +
  theme_classic()


#####MODULE DARKRED
#As a sanity check, let's use ggplot to see what module darkred eigengene looks like between treatment groups.
#First we need to set up the module eigengene for this module with the sample metadata labels we need.
module_darkred_df <- MEs %>%
  tibble::rownames_to_column("FishName") %>%
  # Here we are performing an inner join with a subset of metadata
  dplyr::inner_join(metaData %>%
                      dplyr::select(sample_ID, combined_for_pairwise),
                    by = c("FishName" = "sample_ID")
  )


ggplot(
  module_darkred_df,
  aes(
    x = treatment,
    y = MEdarkred,
    color = treatment
  )
) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  ggforce::geom_sina(maxwidth = 0.3) +
  scale_color_manual(
    values = c("tide" = "blue", "tide_runoff" = "red")
  ) +
  theme_classic()





#################################################################



####PROCEED WITH GO ANALYSIS

###GO Analysis for salmon, red, darkred, darkgreen, blue, pink, white, and turquoise


### GO Enrichment Analysis ### 
# We have identified several interesting modules that are associated with treatment. 
# I'd like to know more about the function of genes in these modules, so I'll perform 
# gene ontology enrichment. I'm going to do this through a modification of code available 
# from M. Matz's GO_MWU pipeline (https://github.com/z0on/GO_MWU)
# Read in large gene module membership file
GMM = read.csv("gfk-master-gene-module-membership-binary.csv")


#######################MODULE DARKGREEN
color = GMM[, c("gene_ID", "darkgreen")]
write.csv(color, file = "gill_darkgreen.csv", quote = F, row.names = F)




######DARKGREEN - MF
# Look for GO categories significantly over-represented among module members, using Fisher's Exact Test
input = "gill_darkgreen.csv"
goAnnotations = "Mummichog_aggregated_GO_terms_filtered.txt"
goDatabase = "go.obo"
goDivision = "MF"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
           perlPath = "perl", 
           largest = 0.1, 
           smallest = 5,
           clusterCutHeight = 0.25)

#4 GO terms at 10% FDR
quartz()
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #absValue=0.001, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.05, # FDR cutoff to print in regular (not italic) font.
                  level3=0.01, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

# text representation of results, with actual adjusted p-values
results[[1]]

    #go_results <- results[[1]]
    #write.csv(go_results, "darkgreen_GO_MF_results.csv", row.names = FALSE)



######DARKGREEN - BP
# Look for GO categories significantly over-represented among module members, using Fisher's Exact Test
input = "gill_darkgreen.csv"
goAnnotations = "Mummichog_aggregated_GO_terms_filtered.txt"
goDatabase = "go.obo"
goDivision = "BP"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
           perlPath = "perl", 
           largest = 0.1, 
           smallest = 5,
           clusterCutHeight = 0.25)

#6 GO terms at 10% FDR

                #file for the enriched GO terms looks like: MWU_BP_gill_darkgreen.csv
quartz()
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #absValue=0.001, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.05, # FDR cutoff to print in regular (not italic) font.
                  level3=0.01, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

# text representation of results, with actual adjusted p-values
results[[1]]



######DARKGREEN - CC
# Look for GO categories significantly over-represented among module members, using Fisher's Exact Test
input = "gill_darkgreen.csv"
goAnnotations = "Mummichog_aggregated_GO_terms_filtered.txt"
goDatabase = "go.obo"
goDivision = "CC"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
           perlPath = "perl", 
           largest = 0.1, 
           smallest = 5,
           clusterCutHeight = 0.25)

#2 GO terms at 10% FDR

#file for the enriched GO terms looks like: MWU_BP_gill_darkgreen.csv
quartz()
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #absValue=0.001, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.05, # FDR cutoff to print in regular (not italic) font.
                  level3=0.01, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

# text representation of results, with actual adjusted p-values
results[[1]]



#######################MODULE DARKRED

color = GMM[, c("gene_ID", "darkred")]
write.csv(color, file = "gill_darkred.csv", quote = F, row.names = F)




######DARKRED - MF
# Look for GO categories significantly over-represented among module members, using Fisher's Exact Test
input = "gill_darkred.csv"
goAnnotations = "Mummichog_aggregated_GO_terms_filtered.txt"
goDatabase = "go.obo"
goDivision = "MF"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
           perlPath = "perl", 
           largest = 0.1, 
           smallest = 5,
           clusterCutHeight = 0.25)

#4 GO terms at 10% FDR
quartz()
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #absValue=0.001, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.05, # FDR cutoff to print in regular (not italic) font.
                  level3=0.01, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

# text representation of results, with actual adjusted p-values
results[[1]]



######DARKRED - BP
# Look for GO categories significantly over-represented among module members, using Fisher's Exact Test
input = "gill_darkred.csv"
goAnnotations = "Mummichog_aggregated_GO_terms_filtered.txt"
goDatabase = "go.obo"
goDivision = "BP"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
           perlPath = "perl", 
           largest = 0.1, 
           smallest = 5,
           clusterCutHeight = 0.25)

#10 GO terms at 10% FDR

#file for the enriched GO terms looks like: MWU_BP_gill_darkgreen.csv
quartz()
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #absValue=0.001, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.05, # FDR cutoff to print in regular (not italic) font.
                  level3=0.01, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

# text representation of results, with actual adjusted p-values
results[[1]]



######DARKRED - CC
# Look for GO categories significantly over-represented among module members, using Fisher's Exact Test
input = "gill_darkred.csv"
goAnnotations = "Mummichog_aggregated_GO_terms_filtered.txt"
goDatabase = "go.obo"
goDivision = "CC"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
           perlPath = "perl", 
           largest = 0.1, 
           smallest = 5,
           clusterCutHeight = 0.25)

#0 GO terms at 10% FDR
    #none passed the enrichment test under the parameters set

quartz()
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #absValue=0.001, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.05, # FDR cutoff to print in regular (not italic) font.
                  level3=0.01, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

# text representation of results, with actual adjusted p-values
results[[1]]





# We have identified several (six) interesting modules that are associated with the interaction between population and treatment. 
# I'd like to know more about the function of genes in these modules, so I'll perform 
# gene ontology enrichment. I'm going to do this through a modification of code available 
# from M. Matz's GO_MWU pipeline (https://github.com/z0on/GO_MWU)
# Read in large gene module membership file
GMM = read.csv("gfk-master-gene-module-membership-binary.csv")

#######################MODULE SALMON

color = GMM[, c("gene_ID", "salmon")]
write.csv(color, file = "gill_salmon.csv", quote = F, row.names = F)



######SALMON - MF
# Look for GO categories significantly over-represented among module members, using Fisher's Exact Test
input = "gill_salmon.csv"
goAnnotations = "Mummichog_aggregated_GO_terms_filtered.txt"
goDatabase = "go.obo"
goDivision = "MF"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
           perlPath = "perl", 
           largest = 0.1, 
           smallest = 5,
           clusterCutHeight = 0.25)

#14 GO terms at 10% FDR

quartz()
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #absValue=0.001, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.05, # FDR cutoff to print in regular (not italic) font.
                  level3=0.01, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

# text representation of results, with actual adjusted p-values
results[[1]]



######SALMON - BP
# Look for GO categories significantly over-represented among module members, using Fisher's Exact Test
input = "gill_salmon.csv"
goAnnotations = "Mummichog_aggregated_GO_terms_filtered.txt"
goDatabase = "go.obo"
goDivision = "BP"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
           perlPath = "perl", 
           largest = 0.1, 
           smallest = 5,
           clusterCutHeight = 0.25)

#48 GO terms at 10% FDR

quartz()
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #absValue=0.001, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.05, # FDR cutoff to print in regular (not italic) font.
                  level3=0.01, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

# text representation of results, with actual adjusted p-values
results[[1]]


######SALMON - CC
# Look for GO categories significantly over-represented among module members, using Fisher's Exact Test
input = "gill_salmon.csv"
goAnnotations = "Mummichog_aggregated_GO_terms_filtered.txt"
goDatabase = "go.obo"
goDivision = "CC"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
           perlPath = "perl", 
           largest = 0.1, 
           smallest = 5,
           clusterCutHeight = 0.25)

#18 GO terms at 10% FDR

quartz()
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #absValue=0.001, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.05, # FDR cutoff to print in regular (not italic) font.
                  level3=0.01, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

# text representation of results, with actual adjusted p-values
results[[1]]



#######################MODULE RED

color = GMM[, c("gene_ID", "red")]
write.csv(color, file = "gill_red.csv", quote = F, row.names = F)



######RED - MF
# Look for GO categories significantly over-represented among module members, using Fisher's Exact Test
input = "gill_red.csv"
goAnnotations = "Mummichog_aggregated_GO_terms_filtered.txt"
goDatabase = "go.obo"
goDivision = "MF"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
           perlPath = "perl", 
           largest = 0.1, 
           smallest = 5,
           clusterCutHeight = 0.25)

#37 GO terms at 10% FDR

quartz()
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #absValue=0.001, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.05, # FDR cutoff to print in regular (not italic) font.
                  level3=0.01, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

# text representation of results, with actual adjusted p-values
results[[1]]



######RED - BP
# Look for GO categories significantly over-represented among module members, using Fisher's Exact Test
input = "gill_red.csv"
goAnnotations = "Mummichog_aggregated_GO_terms_filtered.txt"
goDatabase = "go.obo"
goDivision = "BP"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
           perlPath = "perl", 
           largest = 0.1, 
           smallest = 5,
           clusterCutHeight = 0.25)

#118 GO terms at 10% FDR

quartz()
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #absValue=0.001, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.05, # FDR cutoff to print in regular (not italic) font.
                  level3=0.01, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

# text representation of results, with actual adjusted p-values
results[[1]]



######RED - CC
# Look for GO categories significantly over-represented among module members, using Fisher's Exact Test
input = "gill_red.csv"
goAnnotations = "Mummichog_aggregated_GO_terms_filtered.txt"
goDatabase = "go.obo"
goDivision = "CC"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
           perlPath = "perl", 
           largest = 0.1, 
           smallest = 5,
           clusterCutHeight = 0.25)

#21 GO terms at 10% FDR

quartz()
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #absValue=0.001, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.05, # FDR cutoff to print in regular (not italic) font.
                  level3=0.01, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

# text representation of results, with actual adjusted p-values
results[[1]]



#######################MODULE BLUE

color = GMM[, c("gene_ID", "blue")]
write.csv(color, file = "gill_blue.csv", quote = F, row.names = F)



######BLUE - MF
# Look for GO categories significantly over-represented among module members, using Fisher's Exact Test
input = "gill_blue.csv"
goAnnotations = "Mummichog_aggregated_GO_terms_filtered.txt"
goDatabase = "go.obo"
goDivision = "MF"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
           perlPath = "perl", 
           largest = 0.1, 
           smallest = 5,
           clusterCutHeight = 0.25)

#10 GO terms at 10% FDR

quartz()
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #absValue=0.001, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.05, # FDR cutoff to print in regular (not italic) font.
                  level3=0.01, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

# text representation of results, with actual adjusted p-values
results[[1]]



######BLUE - BP
# Look for GO categories significantly over-represented among module members, using Fisher's Exact Test
input = "gill_blue.csv"
goAnnotations = "Mummichog_aggregated_GO_terms_filtered.txt"
goDatabase = "go.obo"
goDivision = "BP"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
           perlPath = "perl", 
           largest = 0.1, 
           smallest = 5,
           clusterCutHeight = 0.25)

#12 GO terms at 10% FDR

quartz()
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #absValue=0.001, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.05, # FDR cutoff to print in regular (not italic) font.
                  level3=0.01, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

# text representation of results, with actual adjusted p-values
results[[1]]



######BLUE - CC
# Look for GO categories significantly over-represented among module members, using Fisher's Exact Test
input = "gill_blue.csv"
goAnnotations = "Mummichog_aggregated_GO_terms_filtered.txt"
goDatabase = "go.obo"
goDivision = "CC"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
           perlPath = "perl", 
           largest = 0.1, 
           smallest = 5,
           clusterCutHeight = 0.25)

#15 GO terms at 10% FDR

quartz()
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #absValue=0.001, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.05, # FDR cutoff to print in regular (not italic) font.
                  level3=0.01, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

# text representation of results, with actual adjusted p-values
results[[1]]




#######################MODULE PINK

color = GMM[, c("gene_ID", "pink")]
write.csv(color, file = "gill_pink.csv", quote = F, row.names = F)



######PINK - MF
# Look for GO categories significantly over-represented among module members, using Fisher's Exact Test
input = "gill_pink.csv"
goAnnotations = "Mummichog_aggregated_GO_terms_filtered.txt"
goDatabase = "go.obo"
goDivision = "MF"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
           perlPath = "perl", 
           largest = 0.1, 
           smallest = 5,
           clusterCutHeight = 0.25)

#8 GO terms at 10% FDR

quartz()
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #absValue=0.001, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.05, # FDR cutoff to print in regular (not italic) font.
                  level3=0.01, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

# text representation of results, with actual adjusted p-values
results[[1]]



######PINK - BP
# Look for GO categories significantly over-represented among module members, using Fisher's Exact Test
input = "gill_pink.csv"
goAnnotations = "Mummichog_aggregated_GO_terms_filtered.txt"
goDatabase = "go.obo"
goDivision = "BP"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
           perlPath = "perl", 
           largest = 0.1, 
           smallest = 5,
           clusterCutHeight = 0.25)

#1 GO terms at 10% FDR

quartz()
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #absValue=0.001, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.05, # FDR cutoff to print in regular (not italic) font.
                  level3=0.01, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

# text representation of results, with actual adjusted p-values
results[[1]]



######PINK - CC
# Look for GO categories significantly over-represented among module members, using Fisher's Exact Test
input = "gill_pink.csv"
goAnnotations = "Mummichog_aggregated_GO_terms_filtered.txt"
goDatabase = "go.obo"
goDivision = "CC"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
           perlPath = "perl", 
           largest = 0.1, 
           smallest = 5,
           clusterCutHeight = 0.25)

#0 GO terms at 10% FDR

quartz()
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #absValue=0.001, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.05, # FDR cutoff to print in regular (not italic) font.
                  level3=0.01, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

# text representation of results, with actual adjusted p-values
results[[1]]



#######################MODULE WHITE

color = GMM[, c("gene_ID", "white")]
write.csv(color, file = "gill_white.csv", quote = F, row.names = F)



######WHITE - MF
# Look for GO categories significantly over-represented among module members, using Fisher's Exact Test
input = "gill_white.csv"
goAnnotations = "Mummichog_aggregated_GO_terms_filtered.txt"
goDatabase = "go.obo"
goDivision = "MF"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
           perlPath = "perl", 
           largest = 0.1, 
           smallest = 5,
           clusterCutHeight = 0.25)

#0 GO terms at 10% FDR

quartz()
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #absValue=0.001, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.05, # FDR cutoff to print in regular (not italic) font.
                  level3=0.01, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

# text representation of results, with actual adjusted p-values
results[[1]]



######WHITE - BP
# Look for GO categories significantly over-represented among module members, using Fisher's Exact Test
input = "gill_white.csv"
goAnnotations = "Mummichog_aggregated_GO_terms_filtered.txt"
goDatabase = "go.obo"
goDivision = "BP"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
           perlPath = "perl", 
           largest = 0.1, 
           smallest = 5,
           clusterCutHeight = 0.25)

#2 GO terms at 10% FDR

quartz()
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #absValue=0.001, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.05, # FDR cutoff to print in regular (not italic) font.
                  level3=0.01, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

# text representation of results, with actual adjusted p-values
results[[1]]



######WHITE - CC
# Look for GO categories significantly over-represented among module members, using Fisher's Exact Test
input = "gill_white.csv"
goAnnotations = "Mummichog_aggregated_GO_terms_filtered.txt"
goDatabase = "go.obo"
goDivision = "CC"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
           perlPath = "perl", 
           largest = 0.1, 
           smallest = 5,
           clusterCutHeight = 0.25)

#1 GO terms at 10% FDR

quartz()
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #absValue=0.001, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.05, # FDR cutoff to print in regular (not italic) font.
                  level3=0.01, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

# text representation of results, with actual adjusted p-values
results[[1]]



#######################MODULE TURQUOISE

color = GMM[, c("gene_ID", "turquoise")]
write.csv(color, file = "gill_turquoise.csv", quote = F, row.names = F)



######TURQUOISE - MF
# Look for GO categories significantly over-represented among module members, using Fisher's Exact Test
input = "gill_turquoise.csv"
goAnnotations = "Mummichog_aggregated_GO_terms_filtered.txt"
goDatabase = "go.obo"
goDivision = "MF"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
           perlPath = "perl", 
           largest = 0.1, 
           smallest = 5,
           clusterCutHeight = 0.25)

#25 GO terms at 10% FDR

quartz()
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #absValue=0.001, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.05, # FDR cutoff to print in regular (not italic) font.
                  level3=0.01, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

# text representation of results, with actual adjusted p-values
results[[1]]



######TURQUOISE - BP
# Look for GO categories significantly over-represented among module members, using Fisher's Exact Test
input = "gill_turquoise.csv"
goAnnotations = "Mummichog_aggregated_GO_terms_filtered.txt"
goDatabase = "go.obo"
goDivision = "BP"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
           perlPath = "perl", 
           largest = 0.1, 
           smallest = 5,
           clusterCutHeight = 0.25)

#50 GO terms at 10% FDR

quartz()
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #absValue=0.001, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.05, # FDR cutoff to print in regular (not italic) font.
                  level3=0.01, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

# text representation of results, with actual adjusted p-values
results[[1]]



######TURQUOISE - CC
# Look for GO categories significantly over-represented among module members, using Fisher's Exact Test
input = "gill_turquoise.csv"
goAnnotations = "Mummichog_aggregated_GO_terms_filtered.txt"
goDatabase = "go.obo"
goDivision = "CC"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
           perlPath = "perl", 
           largest = 0.1, 
           smallest = 5,
           clusterCutHeight = 0.25)

#^this generates the file you need with the GO terms

#28 GO terms at 10% FDR

quartz()
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #absValue=0.001, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.05, # FDR cutoff to print in regular (not italic) font.
                  level3=0.01, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

# text representation of results, with actual adjusted p-values
results[[1]]


#######################MODULE TAN (almost significantly correlated with "origin of population" trait)

color = GMM[, c("gene_ID", "tan")]
write.csv(color, file = "gill_tan.csv", quote = F, row.names = F)



######TAN - MF
# Look for GO categories significantly over-represented among module members, using Fisher's Exact Test
input = "gill_tan.csv"
goAnnotations = "Mummichog_aggregated_GO_terms_filtered.txt"
goDatabase = "go.obo"
goDivision = "MF"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
           perlPath = "perl", 
           largest = 0.1, 
           smallest = 5,
           clusterCutHeight = 0.25)

#6 GO terms at 10% FDR

quartz()
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #absValue=0.001, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.05, # FDR cutoff to print in regular (not italic) font.
                  level3=0.01, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

# text representation of results, with actual adjusted p-values
results[[1]]



######TAN - BP
# Look for GO categories significantly over-represented among module members, using Fisher's Exact Test
input = "gill_tan.csv"
goAnnotations = "Mummichog_aggregated_GO_terms_filtered.txt"
goDatabase = "go.obo"
goDivision = "BP"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
           perlPath = "perl", 
           largest = 0.1, 
           smallest = 5,
           clusterCutHeight = 0.25)

#4 GO terms at 10% FDR

quartz()
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #absValue=0.001, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.05, # FDR cutoff to print in regular (not italic) font.
                  level3=0.01, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

# text representation of results, with actual adjusted p-values
results[[1]]



######TAN - CC
# Look for GO categories significantly over-represented among module members, using Fisher's Exact Test
input = "gill_tan.csv"
goAnnotations = "Mummichog_aggregated_GO_terms_filtered.txt"
goDatabase = "go.obo"
goDivision = "CC"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision, 
           perlPath = "perl", 
           largest = 0.1, 
           smallest = 5,
           clusterCutHeight = 0.25)


#0 GO terms at 10% FDR

quartz()
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #absValue=0.001, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.05, # FDR cutoff to print in regular (not italic) font.
                  level3=0.01, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

# text representation of results, with actual adjusted p-values
results[[1]]


               
###NOTE: compile all GO terms for a certain module into an excel sheet
###NOTE: you will have to reformat the table because the contents will be bunched together
###NOTE: once that is done, you will need to look at the padj and look for signficantly
#########enriched GO terms. They will be significant at the raw p-value but also significant
########after applying the FDR correction (GO-MWU); these are the terms you need
########If it has a 1 in the padj slot that means they did not survive the FDR correction of 0.1(10%)
####NOTE: it can say 14 GO terms passed 10% FDR but when you put them all in excel, you get 8.
#########This is normal. It removes redundancy during filtering so you can report how many
#########passed the 10% FDR but this many were unique to reduce redundancy.

################################################################################PROCEED BEYOND THIS POINT AT YOUR OWN RISK

### Intramodular analysis: identifying genes with high GS and MM ###
# Using the GS and MM measures, we can identify genes that have a high significance for the trait of interest
# as well as high module membership in interesting modules. This must be done module-by-module and trait-by-trait. 
# Plot shows that genes which are highly significantly associated with a trait (GS) 
# are often also the most important (central) elements of modules associated with that trait. 

dim(datExpr)   # should be samples x genes
    #[1] 17108    37
dim(MEs)       # should be samples x modules
    #[1] 37 24

  ####17108 genes x 37 samples
  ####37 samples x 24 modules
          #this is a mismatch because datExpr has genes in rows and samples in columns, but WGCNA functions expect samples to be in rows
          #and genes in to be in columns. We need to transpose datExpr before calculating correlation:


#Transpose expression data: make it samples x genes
datExpr_t = t(as.matrix(datExpr))
dim(datExpr_t)
    #[1]    37 17108
    #now this reads 37 samples as rows and genes in columns. Proceed with kME calculation

#Calculate module membership (kME)
geneModuleMembership = as.data.frame(cor(datExpr_t, MEs, use = "p"))

#Get p-values for the correlations
nSamples = nrow(datExpr_t)  # should be 37
    #37L shows in the working environment
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

#Calculate gene significance (GS) for all traits
# traitData: data frame with samples in rows, columns = traits
# Check rownames of expression and trait data (these must match or you will get error)
head(rownames(datExpr_t))
    #[1] "L1C1" "L1C2" "L1C3" "L1C4" "L1C5" "L2U1"
head(rownames(traitData))
    #[1] "1" "2" "3" "4" "5" "6"
#Assuming your data is in order (mine is), you can just make the traitData = to datExpr_t
rownames(traitData) = rownames(datExpr_t)
head(rownames(datExpr_t))
    #[1] "L1C1" "L1C2" "L1C3" "L1C4" "L1C5" "L2U1"
head(rownames(traitData))
    #[1] "L1C1" "L1C2" "L1C3" "L1C4" "L1C5" "L2U1"

#this is good
      
      #this will work if all of your columns are numeric; however, it considers my sample column so follow code underneath
      #geneSignificance = as.data.frame(cor(datExpr_t, traitData, use = "p"))
      #GS_pvalue = as.data.frame(corPvalueStudent(as.matrix(geneSignificance), nSamples))

      ###error/warning here again because all columns are not numeric. It is considering my sample_ID column
      ###fix:
        # Exclude sample_ID column
        traitData_numeric = traitData[, c("population", "treatment", "salinity")]

        # Now calculate gene significance (GS)
        geneSignificance = as.data.frame(cor(datExpr_t, traitData_numeric, use = "p"))
        GS_pvalue = as.data.frame(corPvalueStudent(as.matrix(geneSignificance), nSamples))
               
#########MODULE PINK (significantly associated with treatment variable in my dataset which is the population x treatment interaction)
        
# Choose module and trait
module = "pink"
trait = "treatment" #this is the trait that corresponds to population x treatment interaction in my dataset

#Find the ME column for this module
column = match(paste0("ME", module), colnames(MEs))

#Get genes that belong to this module
moduleGenes = moduleColors == module

# Build a table for genes in this module
intramodular = data.frame(
  gene = colnames(datExpr_t)[moduleGenes],     # gene IDs
  MM   = geneModuleMembership[moduleGenes, column],   # module membership
  GS   = geneSignificance[moduleGenes, trait],         # gene significance
  MM_pvalue = MMPvalue[moduleGenes, column],
  GS_pvalue = GS_pvalue[moduleGenes, trait]
)

#Add adjusted p-values
intramodular$MM_padjust = p.adjust(intramodular$MM_pvalue, method = "BH") #Benjamini-Hochberg correction applied
intramodular$GS_padjust = p.adjust(intramodular$GS_pvalue, method = "BH") #Benjamini-Hochberg correction applied


# Plot MM vs GS
sizeGrWindow(7,7)
par(mfrow = c(1,1))
verboseScatterplot(
  abs(geneModuleMembership[moduleGenes, column]),
  abs(geneControlSignificance[moduleGenes, 1]),
  xlab = paste("Module membership in", module, "module"),
  ylab = paste("Gene significance for", trait),
  main = paste("Module membership vs. gene significance\n"),
  cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2,
  col = module
)


####extract top hub genes for the MEpink module


# Rank by high MM and GS
topHubGenes = intramodular[order(-abs(intramodular$MM), -abs(intramodular$GS)), ]

# Top 10 hub genes
head(topHubGenes, 10)
                       #gene        MM         GS    MM_pvalue    GS_pvalue   MM_padjust GS_padjust
#371 gene:ENSFHEG00000020417 0.9085793 -0.3777578 7.883118e-15 0.0211613151 4.934832e-12 0.14398895
#172 gene:ENSFHEG00000023379 0.8771235 -0.3723898 1.076765e-12 0.0232234842 3.370276e-10 0.15000847
#163 gene:ENSFHEG00000002020 0.8662670 -0.3836307 4.330560e-12 0.0190822416 9.036435e-10 0.13272759
#321 gene:ENSFHEG00000001783 0.8612490 -0.3950253 7.913847e-12 0.0155333323 1.238517e-09 0.12185890
#352 gene:ENSFHEG00000001105 0.8399881 -0.4039462 8.029432e-11 0.0131582701 9.665396e-09 0.11202356
#146 gene:ENSFHEG00000011828 0.8385647 -0.3224544 9.263958e-11 0.0516021618 9.665396e-09 0.20575129
#468 gene:ENSFHEG00000020278 0.8290693 -0.2995729 2.324870e-10 0.0716539172 2.079098e-08 0.23962276
#528 gene:ENSFHEG00000017282 0.8256888 -0.5589604 3.182913e-10 0.0003231956 2.356861e-08 0.04261611
#548 gene:ENSFHEG00000008858 0.8250066 -0.3491964 3.388458e-10 0.0341437213 2.356861e-08 0.16986913
#203 gene:ENSFHEG00000013010 0.8221956 -0.2603709 4.372921e-10 0.1196136583 2.737448e-08 0.32288711


# Save all hub genes for pink module related to treatment (population x treatment interaction)
write.csv(topHubGenes, "HubGenes_pink_treatment.csv", row.names = FALSE)

#Gives you a list of genes that are highly correlated for reporting
hubGenes_filtered = subset(topHubGenes, abs(MM) > 0.8 & abs(GS) > 0.2)
head(hubGenes_filtered)
                       #gene        MM         GS    MM_pvalue  GS_pvalue   MM_padjust GS_padjust
#371 gene:ENSFHEG00000020417 0.9085793 -0.3777578 7.883118e-15 0.02116132 4.934832e-12  0.1439889
#172 gene:ENSFHEG00000023379 0.8771235 -0.3723898 1.076765e-12 0.02322348 3.370276e-10  0.1500085
#163 gene:ENSFHEG00000002020 0.8662670 -0.3836307 4.330560e-12 0.01908224 9.036435e-10  0.1327276
#321 gene:ENSFHEG00000001783 0.8612490 -0.3950253 7.913847e-12 0.01553333 1.238517e-09  0.1218589
#352 gene:ENSFHEG00000001105 0.8399881 -0.4039462 8.029432e-11 0.01315827 9.665396e-09  0.1120236
#146 gene:ENSFHEG00000011828 0.8385647 -0.3224544 9.263958e-11 0.05160216 9.665396e-09  0.2057513

write.csv(hubGenes_filtered, "HubGenes_filtered_pink_treatment.csv", row.names = FALSE)

               
               
#########MODULE DARKGREEN (most significantly associated with salinity trait (trait corresponds to treatment in my dataset)

# Choose module and trait
module = "darkgreen"
trait = "salinity" #this is the trait that corresponds to the type of treatment in my dataset

#Find the ME column for this module
column = match(paste0("ME", module), colnames(MEs))

#Get genes that belong to this module
moduleGenes = moduleColors == module

# Build a table for genes in this module
intramodular = data.frame(
  gene = colnames(datExpr_t)[moduleGenes],     # gene IDs
  MM   = geneModuleMembership[moduleGenes, column],   # module membership
  GS   = geneSignificance[moduleGenes, trait],         # gene significance
  MM_pvalue = MMPvalue[moduleGenes, column],
  GS_pvalue = GS_pvalue[moduleGenes, trait]
)

#Add adjusted p-values
intramodular$MM_padjust = p.adjust(intramodular$MM_pvalue, method = "BH") #Benjamini-Hochberg correction applied
intramodular$GS_padjust = p.adjust(intramodular$GS_pvalue, method = "BH") #Benjamini-Hochberg correction applied


# Plot MM vs GS
sizeGrWindow(7,7)
par(mfrow = c(1,1))
verboseScatterplot(
  abs(geneModuleMembership[moduleGenes, column]),
  abs(geneControlSignificance[moduleGenes, 1]),
  xlab = paste("Module membership in", module, "module"),
  ylab = paste("Gene significance for", trait),
  main = paste("Module membership vs. gene significance\n"),
  cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2,
  col = module
)

####extract top hub genes for the MEdarkgreen module


# Rank by high MM and GS
topHubGenes = intramodular[order(-abs(intramodular$MM), -abs(intramodular$GS)), ]

# Top 10 hub genes
head(topHubGenes, 10)
                      #gene        MM        GS    MM_pvalue   GS_pvalue   MM_padjust GS_padjust
#654 gene:ENSFHEG00000019007 0.8733233 0.3445413 1.778118e-12 0.036775338 1.433163e-09 0.14319286
#7   gene:ENSFHEG00000006701 0.8652118 0.4348735 4.925714e-12 0.007150419 1.985063e-09 0.06003373
#80  gene:ENSFHEG00000020075 0.8603864 0.4111710 8.757363e-12 0.011466916 2.352812e-09 0.07453495
#612 gene:ENSFHEG00000008645 0.8505435 0.4186400 2.657491e-11 0.009915994 5.354844e-09 0.06889906
#143 gene:ENSFHEG00000002699 0.8378052 0.3084076 9.992912e-11 0.063291057 1.516744e-08 0.18550034
#336 gene:ENSFHEG00000002411 0.8365725 0.4427553 1.129090e-10 0.006066185 1.516744e-08 0.05671409
#505 gene:ENSFHEG00000005692 0.8324995 0.2553721 1.678452e-10 0.127141283 1.912105e-08 0.28160163
#242 gene:ENSFHEG00000013401 0.8309436 0.3223702 1.947506e-10 0.051666709 1.912105e-08 0.16791680
#497 gene:ENSFHEG00000014664 0.8299732 0.4173371 2.135105e-10 0.010172897 1.912105e-08 0.06890214
#726 gene:ENSFHEG00000011919 0.8272840 0.3585879 2.746701e-10 0.029302145 2.038839e-08 0.12766232


# Save all hub genes for darkgreen module related to salinity (treatment exposure)
write.csv(topHubGenes, "HubGenes_darkgreen_salinity.csv", row.names = FALSE)

#gives you a list of genes that are highly correlated for reporting
hubGenes_filtered = subset(topHubGenes, abs(MM) > 0.8 & abs(GS) > 0.2)
head(hubGenes_filtered)
                      #gene        MM        GS
#654 gene:ENSFHEG00000019007 0.8733233 0.3445413
#7   gene:ENSFHEG00000006701 0.8652118 0.4348735
#80  gene:ENSFHEG00000020075 0.8603864 0.4111710
#612 gene:ENSFHEG00000008645 0.8505435 0.4186400
#143 gene:ENSFHEG00000002699 0.8378052 0.3084076
#336 gene:ENSFHEG00000002411 0.8365725 0.4427553

write.csv(hubGenes_filtered, "HubGenes_filtered_darkgreen_salinity.csv", row.names = FALSE)

               
               
#########MODULE TAN (almost significantly associated with population trait)

# Choose module and trait
module = "tan"
trait = "population" #this is the trait that corresponds to origin of population

#Find the ME column for this module
column = match(paste0("ME", module), colnames(MEs))

#Get genes that belong to this module
moduleGenes = moduleColors == module

# Build a table for genes in this module
intramodular = data.frame(
  gene = colnames(datExpr_t)[moduleGenes],     # gene IDs
  MM   = geneModuleMembership[moduleGenes, column],   # module membership
  GS   = geneSignificance[moduleGenes, trait],         # gene significance
  MM_pvalue = MMPvalue[moduleGenes, column],
  GS_pvalue = GS_pvalue[moduleGenes, trait]
)

#Add adjusted p-values
intramodular$MM_padjust = p.adjust(intramodular$MM_pvalue, method = "BH") #Benjamini-Hochberg correction applied
intramodular$GS_padjust = p.adjust(intramodular$GS_pvalue, method = "BH") #Benjamini-Hochberg correction applied


# Plot MM vs GS
sizeGrWindow(7,7)
par(mfrow = c(1,1))
verboseScatterplot(
  abs(geneModuleMembership[moduleGenes, column]),
  abs(geneControlSignificance[moduleGenes, 1]),
  xlab = paste("Module membership in", module, "module"),
  ylab = paste("Gene significance for", trait),
  main = paste("Module membership vs. gene significance\n"),
  cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2,
  col = module
)
####extract top hub genes for the MEtan module


# Rank by high MM and GS
topHubGenes = intramodular[order(-abs(intramodular$MM), -abs(intramodular$GS)), ]

# Top 10 hub genes
head(topHubGenes, 10)
                       #gene        MM        GS    MM_pvalue   GS_pvalue   MM_padjust GS_padjust
#151 gene:ENSFHEG00000006981 0.9121473 0.2378776 4.043342e-15 0.156268304 1.665857e-12  0.4075083
#18  gene:ENSFHEG00000023306 0.8884549 0.1189022 2.174402e-13 0.483351389 4.479268e-11  0.6943905
#121 gene:ENSFHEG00000021993 0.8707211 0.2854422 2.483985e-12 0.086801395 3.411339e-10  0.3311592
#125 gene:ENSFHEG00000022652 0.8642861 0.3142653 5.509988e-12 0.058185370 5.675288e-10  0.2830428
#299 gene:ENSFHEG00000014685 0.8568540 0.2641715 1.316711e-11 0.114117576 1.084970e-09  0.3702082
#155 gene:ENSFHEG00000021672 0.8498176 0.2599087 2.875139e-11 0.120295290 1.941050e-09  0.3765483
#270 gene:ENSFHEG00000018640 0.8475635 0.4422910 3.661739e-11 0.006125875 1.941050e-09  0.1402145
#260 gene:ENSFHEG00000017160 0.8472918 0.3779972 3.769029e-11 0.021073015 1.941050e-09  0.1977471
#220 gene:ENSFHEG00000007329 0.8329899 0.1055921 1.601120e-10 0.533948898 7.329570e-09  0.7308536
#30  gene:ENSFHEG00000009617 0.8260525 0.3716563 3.078125e-10 0.023517782 1.171025e-08  0.1992112


# Save all hub genes for pink module related to population (origin of population variable)
write.csv(topHubGenes, "HubGenes_tan_population.csv", row.names = FALSE)

#Gives you a list of genes that are highly correlated for reporting
hubGenes_filtered = subset(topHubGenes, abs(MM) > 0.8 & abs(GS) > 0.2)
head(hubGenes_filtered)
                       #gene        MM        GS    MM_pvalue   GS_pvalue   MM_padjust GS_padjust
#151 gene:ENSFHEG00000006981 0.9121473 0.2378776 4.043342e-15 0.156268304 1.665857e-12  0.4075083
#121 gene:ENSFHEG00000021993 0.8707211 0.2854422 2.483985e-12 0.086801395 3.411339e-10  0.3311592
#125 gene:ENSFHEG00000022652 0.8642861 0.3142653 5.509988e-12 0.058185370 5.675288e-10  0.2830428
#299 gene:ENSFHEG00000014685 0.8568540 0.2641715 1.316711e-11 0.114117576 1.084970e-09  0.3702082
#155 gene:ENSFHEG00000021672 0.8498176 0.2599087 2.875139e-11 0.120295290 1.941050e-09  0.3765483
#270 gene:ENSFHEG00000018640 0.8475635 0.4422910 3.661739e-11 0.006125875 1.941050e-09  0.1402145

write.csv(hubGenes_filtered, "HubGenes_filtered_tan_population.csv", row.names = FALSE)

    ####NOTE:p-values and adjusted p-values were included but are not the key to identifying
    ####hub genes. Hub genes are based on high MM and high GS, not FDR or padj <0.05. They are
    ####there to include in supplemental but note that these are not parameters that tells you
    ####what the hub genes are. Only MM and GS does. Gene signficance can have a + or - value,
    ####as long as it is not 0. Gene signficance shows how biologically relevant it is.

###Annotate genes for .csv with genes and their corresponding modules by merging filtered hubGenes .csv with gene annotation file (mine is Mummichog.gtf)
###Gene Annotation - Merging highly correlated hub genes (HubGenesFiltered) with gene annotation file from the Atlantic killifish - Mummichog

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("rtracklayer")

library(rtracklayer)

# Import GTF
gtf = import("mummichog.gtf")
# Convert to data frame
gtf_df = as.data.frame(gtf)
# Keep only gene entries
gtf_genes = gtf_df[gtf_df$type == "gene", c("gene_id", "gene_name")]


#######Merge hub genes from module pink


#Read your hub gene table
pinkhubGenes = read.csv("HubGenes_filtered_pink_treatment.csv")

            #This code did not work because I had to remove the prefix "gene:" from the annotation file so that
            #the tables would merge...see code right below on filtering out prefix and downstream analyses
            #Merge annotation into hubGenes
            #pinkhubGenes_annotated = merge(
            #pinkhubGenes,
            #gtf_genes,
            #by.x = "gene",     # column name in your CSV
            #by.y = "gene_id",  # column name in GTF
            #all.x = TRUE       # keep all hub genes
            #)


pinkhubGenes$gene_clean = sub("^gene:", "", pinkhubGenes$gene)
gtf_genes$gene_clean = sub("^gene:", "", gtf_genes$gene_id)

pinkhubGenes_annotated = merge(
  pinkhubGenes,
  gtf_genes,
  by = "gene_clean",
  all.x = TRUE
)

write.csv(pinkhubGenes_annotated, "HubGenes_pink_treatment_annotated.csv", row.names = FALSE)


#######Merge hub genes from module darkgreen


#Read your hub gene table
darkgreenhubGenes = read.csv("HubGenes_filtered_darkgreen_salinity.csv")


darkgreenhubGenes$gene_clean = sub("^gene:", "", darkgreenhubGenes$gene)
gtf_genes$gene_clean = sub("^gene:", "", gtf_genes$gene_id)

darkgreenhubGenes_annotated = merge(
  darkgreenhubGenes,
  gtf_genes,
  by = "gene_clean",
  all.x = TRUE
)

write.csv(darkgreenhubGenes_annotated, "HubGenes_darkgreen_salinity_annotated.csv", row.names = FALSE)


#######Merge hub genes from module tan


#Read your hub gene table
tanhubGenes = read.csv("HubGenes_filtered_tan_population.csv")


tanhubGenes$gene_clean = sub("^gene:", "", tanhubGenes$gene)
gtf_genes$gene_clean = sub("^gene:", "", gtf_genes$gene_id)

tanhubGenes_annotated = merge(
  tanhubGenes,
  gtf_genes,
  by = "gene_clean",
  all.x = TRUE
)

write.csv(tanhubGenes_annotated, "HubGenes_tan_population_annotated.csv", row.names = FALSE)



##################Creating boxplots for significant modules in the pop x treatment interaction (treatment trait), but grouping them based on control (0) or stress (1)
library(ggplot2)
# Example: traitData$treatment is 0 = Control, 1 = Stress



#####MODULE SALMON

#Salmon module 
module = "salmon"
MEname = paste0("ME", module)

# Make a dataframe with ME values + treatment
plotData <- data.frame(
  Treatment = as.factor(traitData$treatment),
  Eigengene = MEs[, MEname]
)

#Plot the boxplot

ggplot(plotData, aes(x = Treatment, y = Eigengene, fill = Treatment)) +
  geom_boxplot(alpha = 1, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  labs(
    title = paste("Module eigengene expression:", module),
    x = "Treatment",
    y = "Module Eigengene Value"
  ) +
  scale_fill_manual(values = c("0" = "#0000FF", "1" = "#FF0000")) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),  # removes major horizontal gridlines
    panel.grid.minor.y = element_blank(),  # removes minor horizontal gridlines
    panel.grid.major.x = element_blank(),  # optional: remove vertical lines
    panel.grid.minor.x = element_blank()   # optional: remove vertical lines
  )


#####MODULE RED

#Red module 
module = "red"
MEname = paste0("ME", module)

# Make a dataframe with ME values + treatment
plotData <- data.frame(
  Treatment = as.factor(traitData$treatment),
  Eigengene = MEs[, MEname]
)

#Plot the boxplot

ggplot(plotData, aes(x = Treatment, y = Eigengene, fill = Treatment)) +
  geom_boxplot(alpha = 1, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  labs(
    title = paste("Module eigengene expression:", module),
    x = "Treatment",
    y = "Module Eigengene Value"
  ) +
  scale_fill_manual(values = c("0" = "#0000FF", "1" = "#FF0000")) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),  # removes major horizontal gridlines
    panel.grid.minor.y = element_blank(),  # removes minor horizontal gridlines
    panel.grid.major.x = element_blank(),  # optional: remove vertical lines
    panel.grid.minor.x = element_blank()   # optional: remove vertical lines
  )


#####MODULE BLUE

#Blue module 
module = "blue"
MEname = paste0("ME", module)

# Make a dataframe with ME values + treatment
plotData <- data.frame(
  Treatment = as.factor(traitData$treatment),
  Eigengene = MEs[, MEname]
)

#Plot the boxplot

ggplot(plotData, aes(x = Treatment, y = Eigengene, fill = Treatment)) +
  geom_boxplot(alpha = 1, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  labs(
    title = paste("Module eigengene expression:", module),
    x = "Treatment",
    y = "Module Eigengene Value"
  ) +
  scale_fill_manual(values = c("0" = "#0000FF", "1" = "#FF0000")) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),  # removes major horizontal gridlines
    panel.grid.minor.y = element_blank(),  # removes minor horizontal gridlines
    panel.grid.major.x = element_blank(),  # optional: remove vertical lines
    panel.grid.minor.x = element_blank()   # optional: remove vertical lines
  )


#####MODULE PINK

#Pink module 
module = "pink"
MEname = paste0("ME", module)

# Make a dataframe with ME values + treatment
plotData <- data.frame(
  Treatment = as.factor(traitData$treatment),
  Eigengene = MEs[, MEname]
)

#Plot the boxplot

ggplot(plotData, aes(x = Treatment, y = Eigengene, fill = Treatment)) +
  geom_boxplot(alpha = 1, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  labs(
    title = paste("Module eigengene expression:", module),
    x = "Treatment",
    y = "Module Eigengene Value"
  ) +
  scale_fill_manual(values = c("0" = "#0000FF", "1" = "#FF0000")) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),  # removes major horizontal gridlines
    panel.grid.minor.y = element_blank(),  # removes minor horizontal gridlines
    panel.grid.major.x = element_blank(),  # optional: remove vertical lines
    panel.grid.minor.x = element_blank()   # optional: remove vertical lines
  )


#####MODULE WHITE

#Salmon module 
module = "white"
MEname = paste0("ME", module)

# Make a dataframe with ME values + treatment
plotData <- data.frame(
  Treatment = as.factor(traitData$treatment),
  Eigengene = MEs[, MEname]
)

#Plot the boxplot

ggplot(plotData, aes(x = Treatment, y = Eigengene, fill = Treatment)) +
  geom_boxplot(alpha = 1, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  labs(
    title = paste("Module eigengene expression:", module),
    x = "Treatment",
    y = "Module Eigengene Value"
  ) +
  scale_fill_manual(values = c("0" = "#0000FF", "1" = "#FF0000")) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),  # removes major horizontal gridlines
    panel.grid.minor.y = element_blank(),  # removes minor horizontal gridlines
    panel.grid.major.x = element_blank(),  # optional: remove vertical lines
    panel.grid.minor.x = element_blank()   # optional: remove vertical lines
  )


#####MODULE TURQUOISE

#Salmon module 
module = "turquoise"
MEname = paste0("ME", module)

# Make a dataframe with ME values + treatment
plotData <- data.frame(
  Treatment = as.factor(traitData$treatment),
  Eigengene = MEs[, MEname]
)

#Plot the boxplot

ggplot(plotData, aes(x = Treatment, y = Eigengene, fill = Treatment)) +
  geom_boxplot(alpha = 1, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  labs(
    title = paste("Module eigengene expression:", module),
    x = "Treatment",
    y = "Module Eigengene Value"
  ) +
  scale_fill_manual(values = c("0" = "#0000FF", "1" = "#FF0000")) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),  # removes major horizontal gridlines
    panel.grid.minor.y = element_blank(),  # removes minor horizontal gridlines
    panel.grid.major.x = element_blank(),  # optional: remove vertical lines
    panel.grid.minor.x = element_blank()   # optional: remove vertical lines
  )
              
