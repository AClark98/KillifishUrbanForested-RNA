##############Killifish RNA-seq dataset################
##### Chapter 1 of dissertation - Gene expression analysis continued...-June 2025;R-studio v4.4.1 #####
##Genetic front loading test for a Reciprocal Transplant Experiment##


#####THIS TAKES PLACE AFTER RUNNING DESEq to get differentially expressed genes based on origin (origin + ~ treatment), treatment (~ origin + treatment), and origin+treatment (~1)
setwd("~/Library/CloudStorage/Box-Box/Bernal_lab/Antrelle/GulfKillifish")

library(VennDiagram)
library(gplots)


#Define "control" and "experimental" conditions:
	#long_bayou_F control conditions: long_bayou_F:tide
	#long_bayoou_F experimental conditions: long_bayou_F:tide_runoff
	#weakley_bayou_U control conditions: weakley_bayou_U:tide_runoff
	#weakley_bayou_U experimental conditions: weakley_bayou_U:tide


#extract normalized expression values using vst (for entire dataset; not pop specific)
# Run variance stabilizing transformation
vsd <- vst(GFKdds, blind=FALSE)
vst_counts <- assay(vsd)


#Subset data based on origin to determine effect of treatment

##### Long Bayou population
LB_dds <- GFKdds[, GFKdds$origin == "long_bayou_F"]
design(LB_dds) <- ~ treatment
LB_dds <- DESeq(LB_dds)
LB_treat_res <- results(LB_dds, contrast = c("treatment", "tide_runoff", "tide"))  # foreign vs home


# Get significantly upregulated genes (treatment-induced)
treatment_induced_LB <- subset(LB_treat_res, padj < 0.05 & log2FoldChange > 1)
treatment_induced_genes_LB <- rownames(treatment_induced_LB)


# Long Bayou = long_bayou_F, home treatment = tide
LB_control_samples <- colnames(GFKdds)[
  GFKdds$origin == "long_bayou_F" & GFKdds$treatment == "tide"
]

# Mean expression of experimental-induced genes from Weakley Bayou in Long Bayou home condition
mean_expr_in_LB <- rowMeans(vst_counts[treatment_induced_genes_LB, LB_control_samples, drop=FALSE])

# View top front-loaded candidates
frontloaded_WB_in_LB <- sort(mean_expr_in_LB, decreasing = TRUE) #this means genes that are induced in WB that are already expressed in LB when in their native conditions
head(frontloaded_WB_in_LB, 10)




##### Weakley Bayou population
WB_dds <- GFKdds[, GFKdds$origin == "weakley_bayou_U"]
design(WB_dds) <- ~ treatment
WB_dds <- DESeq(WB_dds)
WB_treat_res <- results(WB_dds, contrast = c("treatment", "tide_runoff", "tide"))  # foreign vs home


# Get significantly upregulated genes (treatment-induced)
treatment_induced_WB <- subset(WB_treat_res, padj < 0.05 & log2FoldChange > 1)
treatment_induced_genes_WB <- rownames(treatment_induced_WB)


# Weakley Bayou = weakley_bayou_U, home treatment = tide_runoff
WB_control_samples <- colnames(GFKdds)[
  GFKdds$origin == "weakley_bayou_U" & GFKdds$treatment == "tide_runoff"
]

# Mean expression of experimental-induced genes from Long Bayou in Weakley Bayou home condition
mean_expr_in_WB <- rowMeans(vst_counts[treatment_induced_genes_WB, WB_control_samples, drop=FALSE])

# View top front-loaded candidates
frontloaded_LB_in_WB <- sort(mean_expr_in_WB, decreasing = TRUE) #this means genes that are induced in LB that are already expressed in WB when in their native conditions
head(frontloaded_LB_in_WB, 10)




