################# WGCNA  - RNA-seq #####################
######### Gulf Killifish - Fundulus grandis - Gill ########

# Antrelle D. Clark - Auburn University Department of Biological Sciences - Fall 2025
# Gratitude to Ally Swank for the WGCNA script
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


## Will display correlations and their p-values
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
