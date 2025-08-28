##############Killifish RNA-seq dataset################
##### Chapter 1 of dissertation - Gene expression analysis continued...-June 2025;R-studio v4.4.1 #####
##Genetic front loading test for a Reciprocal Transplant Experiment##


##############Killifish RNA-seq dataset################
##### Chapter 1 of dissertation - Gene expression analysis continued...-June 2025;R-studio v4.4.1 #####
##Genetic front loading test for a Reciprocal Transplant Experiment##


#####THIS TAKES PLACE IMMEDIATELY AFTER RUNNING DESEQ2 to get differentially expressed genes based on origin (origin + ~ treatment), treatment (~ origin + treatment), and origin+treatment (~1)
setwd("~/Library/CloudStorage/Box-Box/Bernal_lab/Antrelle/GulfKillifish")

library(VennDiagram)
library(gplots)


#Define "control" and "experimental" conditions (this is not a code; just so that you know what is the control and experimental treatment for each population):
#long_bayou_F control conditions: long_bayou_F:tide
#long_bayoou_F experimental conditions: long_bayou_F:tide_runoff
#weakley_bayou_U control conditions: weakley_bayou_U:tide_runoff
#weakley_bayou_U experimental conditions: weakley_bayou_U:tide


#extract normalized expression values using vst (for entire dataset; not pop specific)
# Run variance stabilizing transformation
vsd <- vst(GFKdds, blind=FALSE)
vst_counts <- assay(vsd)

###########BACKLOADING############ #what genes are significantly turned on in response to new stressor or environment

#Subset data based on origin to determine effect of treatment

##### Long Bayou population
LB_dds <- GFKdds[, GFKdds$origin == "long_bayou_F"]
design(LB_dds) <- ~ treatment
LB_dds <- DESeq(LB_dds)
LB_treat_res <- results(LB_dds, contrast = c("treatment", "tide_runoff", "tide"))  # foreign vs home


# Get significantly upregulated genes (turned-on)
turned_on_LB <- subset(LB_treat_res, padj < 0.05 & log2FoldChange > 1) #the >1 means that the genes change expression (backloaded)
turned_on_genes_LB <- rownames(turned_on_LB)


# Long Bayou = long_bayou_F, home treatment = tide
LB_control_samples <- colnames(GFKdds)[
  GFKdds$origin == "long_bayou_F" & GFKdds$treatment == "tide"
]

# Mean expression of experimental-induced genes from Weakley Bayou in Long Bayou home condition
mean_expr_in_LB <- rowMeans(vst_counts[turned_on_genes_LB, LB_control_samples, drop=FALSE])

# View top front-loaded candidates
backloaded_WB_in_LB <- sort(mean_expr_in_LB, decreasing = TRUE)
head(backloaded_WB_in_LB, 10) #shows top 10 genes

    #gene:ENSFHEG00000006314 gene:ENSFHEG00000008149 gene:ENSFHEG00000016601 
                   #9.839421                7.160313                7.020118 

backloaded_WB_in_LB #to see all genes
    #3 genes significantly induced when WB pop is exposed to LB normal conditions (tide treatment)

#extract list of backloaded genes in LB pop
write.csv(backloaded_WB_in_LB, file = "backloaded_WB_in_LB.csv")


##################################################


##### Weakley Bayou population
WB_dds <- GFKdds[, GFKdds$origin == "weakley_bayou_U"]
design(WB_dds) <- ~ treatment
WB_dds <- DESeq(WB_dds)
WB_treat_res <- results(WB_dds, contrast = c("treatment", "tide_runoff", "tide"))  # foreign vs home


# Get significantly upregulated genes (treatment-induced)
turned_on_WB <- subset(WB_treat_res, padj < 0.05 & log2FoldChange > 1)
turned_on_genes_WB <- rownames(turned_on_WB)


# Weakley Bayou = weakley_bayou_U, home treatment = tide_runoff
WB_control_samples <- colnames(GFKdds)[
  GFKdds$origin == "weakley_bayou_U" & GFKdds$treatment == "tide_runoff"
]

# Mean expression of experimental-induced genes from Long Bayou in Weakley Bayou home condition
mean_expr_in_WB <- rowMeans(vst_counts[turned_on_genes_WB, WB_control_samples, drop=FALSE])

# View top back-loaded candidates
backloaded_LB_in_WB <- sort(mean_expr_in_WB, decreasing = TRUE)
        #head(backloaded_LB_in_WB, 10) #to see top 10


        #gene:ENSFHEG00000021576 gene:ENSFHEG00000010562 gene:ENSFHEG00000009304 
                      #11.446843               11.170573               10.020650 
        #gene:ENSFHEG00000003810 gene:ENSFHEG00000020910 gene:ENSFHEG00000006476 
                       #9.179378                8.772974                8.692392 
        #gene:ENSFHEG00000016212 gene:ENSFHEG00000023005 gene:ENSFHEG00000018965 
                       #8.535428                8.383383                8.190001 
        #gene:ENSFHEG00000001392 
                       #7.909710 

backloaded_LB_in_WB

        #gene:ENSFHEG00000021576 gene:ENSFHEG00000010562 gene:ENSFHEG00000009304 
                      #11.446843               11.170573               10.020650 
        #gene:ENSFHEG00000003810 gene:ENSFHEG00000020910 gene:ENSFHEG00000006476 
                       #9.179378                8.772974                8.692392 
        #gene:ENSFHEG00000016212 gene:ENSFHEG00000023005 gene:ENSFHEG00000018965 
                       #8.535428                8.383383                8.190001 
        #gene:ENSFHEG00000001392 gene:ENSFHEG00000015240 gene:ENSFHEG00000014562 
                       #7.909710                7.812889                7.758287 
        #gene:ENSFHEG00000011984 
                       #7.635522

#13 genes significantly induced when LB pop is exposed to WB normal conditions (tide_runoff treatment) 

#extract list of backloaded genes in LB pop
write.csv(backloaded_LB_in_WB, file = "backloaded_LB_in_WB.csv")



##################################################

###Gene Annotation

#Transfer the .gtf and/or .gff file that corresponds to the reference genome (i.e., from EnsemblID or NCBI) (this step was already done so transferred from my linux terminal working directory)
  #scp aubadc002@asax.asc.edu:/home/aubadc002/killifish/fastqc/paired/Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gff3 /Users/adc0064/Library/CloudStorage/Box-Box/Bernal_lab/Antrelle/GulfKillifish/
      #after transfer, renamed Mummichog.gff3
  #scp aubadc002@asax.asc.edu:/home/aubadc002/killifish/fastqc/paired/Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.112.gtf /Users/adc0064/Library/CloudStorage/Box-Box/Bernal_lab/Antrelle/GulfKillifish/
      #after transfer, renamed Mummichog.gtf


### Extract gene_id and gene_name to create mapping table 

# Read the GTF file as a plain text table
gtf_lines <- readLines("Mummichog.gtf")

# Filter lines that describe genes
gene_lines <- gtf_lines[grepl("\tgene\t", gtf_lines)]

# Extract gene_id and gene_name using regular expressions
gene_info <- do.call(rbind, lapply(gene_lines, function(line) {
  gene_id <- sub('.*gene_id "([^"]+)".*', '\\1', line)
  gene_name <- if (grepl('gene_name "', line)) {
    sub('.*gene_name "([^"]+)".*', '\\1', line)
  } else {
    NA  # If gene_name doesn't exist
  }
  return(data.frame(gene_id = gene_id, gene_name = gene_name, stringsAsFactors = FALSE))
}))

# Remove duplicates
gene_info <- gene_info[!duplicated(gene_info$gene_id), ]

# Save to CSV if you want
write.csv(gene_info, "gene_id_to_gene_name_mapping.csv", row.names = FALSE)

#Merge the mapping table with your backloaded gene data

### Long Bayou (induced in WB when exposed to LB control "tide")
gene_ids_WB_in_LB <- names(backloaded_WB_in_LB)
backloaded_df_WB_in_LB <- data.frame(
  gene_id = gene_ids_WB_in_LB,
  expression = backloaded_WB_in_LB
)
annotated_WB_in_LB <- merge(backloaded_df_WB_in_LB, gene_info, by = "gene_id", all.x = TRUE)
write.csv(annotated_WB_in_LB, "annotated_backloaded_WB_in_LB.csv", row.names = FALSE)


### Weakley Bayou (induced in LB when exposed to WB control "tide_runoff")
gene_ids_LB_in_WB <- names(backloaded_LB_in_WB)
backloaded_df_LB_in_WB <- data.frame(
  gene_id = gene_ids_LB_in_WB,
  expression = backloaded_LB_in_WB
)
annotated_LB_in_WB <- merge(backloaded_df_LB_in_WB, gene_info, by = "gene_id", all.x = TRUE)
write.csv(annotated_LB_in_WB, "annotated_backloaded_LB_in_WB.csv", row.names = FALSE)

#^this was done but the genes IDS have "gene:" in front which prevents accurate mapping to get "gene name" so this was redone below to exclude the prefix "gene:"

### Long Bayou (induced in WB when exposed to LB control "tide")
gene_ids_WB_in_LB <- sub("^gene:", "", names(backloaded_WB_in_LB))
backloaded_df_WB_in_LB <- data.frame(
  gene_id = gene_ids_WB_in_LB,
  expression = as.numeric(backloaded_WB_in_LB)
)
annotated_WB_in_LB <- merge(backloaded_df_WB_in_LB, gene_info, by = "gene_id", all.x = TRUE)
write.csv(annotated_WB_in_LB, "annotated_backloaded_WB_in_LB.csv", row.names = FALSE)


### Weakley Bayou (induced in LB when exposed to WB control "tide_runoff")
gene_ids_LB_in_WB <- sub("^gene:", "", names(backloaded_LB_in_WB))
backloaded_df_LB_in_WB <- data.frame(
  gene_id = gene_ids_LB_in_WB,
  expression = as.numeric(backloaded_LB_in_WB)
)
annotated_LB_in_WB <- merge(backloaded_df_LB_in_WB, gene_info, by = "gene_id", all.x = TRUE)
write.csv(annotated_LB_in_WB, "annotated_backloaded_LB_in_WB.csv", row.names = FALSE)

#^this worked. some still have NA in "gene name" which simply means that it is unannotated in the reference(sister species) genome

# at this point, I could annotate the genes or not
    #skipping for now


###############Visualize data using heatmap
  
###HEATMAP###

#VST-normalized expression matrix (normalizing data again by submitting to different output name)
vsd <- vst(GFKdds, blind = FALSE)
vsd_mat <- assay(vsd)


library(pheatmap)

#read in backloading data
genes_WB_in_LB <- names(backloaded_WB_in_LB)
genes_LB_in_WB <- names(backloaded_LB_in_WB)

#Subset data (VST matrix) for each back-loaded gene list
expr_WB_in_LB <- vsd_mat[genes_WB_in_LB, ]
expr_LB_in_WB <- vsd_mat[genes_LB_in_WB, ]


# Add annotations (e.g., origin, treatment)
anno <- as.data.frame(colData(vsd)[, c("origin", "treatment")])

### Heatmap for WB genes back-loaded in LB

# Define your custom colors
my_colors <- list(
  origin = c("long_bayou_F" = "blue", "weakley_bayou_U" = "red"),     # blue for LB and red for WB for origins
  treatment = c("tide" = "blue", "tide_runoff" = "red")  # blue corresponds to LB control and red corresponds to WB control for treatments
)

pheatmap(expr_WB_in_LB,
         annotation_col = anno,
         annotation_colors = my_colors,
         show_rownames = FALSE,
         cluster_cols = FALSE, #keeps samples ordered; I ran both so I have clustering and non-clustered heatmaps
         main = "WB Genes Back-Loaded in LB")

  #heatmap saved as "Backloading_WB_in_LBhome.pdf" and "Backloading_WB_in_LBhome.png"


### Heatmap for LB genes back-loaded in WB
pheatmap(expr_LB_in_WB,
         annotation_col = anno,
         annotation_colors = my_colors,
         show_rownames = FALSE,
         cluster_cols = FALSE, #keeps samples ordered; I ran both so I have clustering and non-clustered heatmaps
         main = "LB Genes Back-Loaded in WB")

  #heatmap saved as "Backloading_LB_in_WBhome.pdf" and "Backloading_LB_in_WBhome.png"







###########FRONTLOADING############ #what genes are turned on prior to exposure to stressor or new environment

####Long Bayou population


# Get significantly genes that has unchanged expression (already-on)
already_on_LB <- subset(LB_treat_res, padj < 0.05 & log2FoldChange < -1)
already_on_genes_LB <- rownames(already_on_LB)


# Mean expression of already on genes from Weakley Bayou in Long Bayou home condition
mean_expr_in_LB <- rowMeans(vst_counts[already_on_genes_LB, LB_control_samples, drop=FALSE])

# View top front-loaded candidates
frontloaded_WB_in_LB <- sort(mean_expr_in_LB, decreasing = TRUE)
head(frontloaded_WB_in_LB, 10) #shows top 10 genes
frontloaded_WB_in_LB

        #gene:ENSFHEG00000017456 gene:ENSFHEG00000009528 gene:ENSFHEG00000014019 
                       #9.679064                8.013500                7.313976 
        #gene:ENSFHEG00000012714 
                       #6.877167 

#4 genes already on in WB and stays at a consistent expression when exposed to LB normal conditions (tide treatment)

#extract list of frontloaded genes in LB pop
write.csv(frontloaded_WB_in_LB, file = "frontloaded_WB_in_LB.csv")


##################################################


##### Weakley Bayou population

# Get significantly genes that has unchanged expression (already-on)
already_on_WB <- subset(WB_treat_res, padj < 0.05 & log2FoldChange < -1)
already_on_genes_WB <- rownames(already_on_WB)

# Mean expression of already on genes from Weakley Bayou in Long Bayou home condition
mean_expr_in_WB <- rowMeans(vst_counts[already_on_genes_WB, WB_control_samples, drop=FALSE])

# View top front-loaded candidates
frontloaded_LB_in_WB <- sort(mean_expr_in_WB, decreasing = TRUE)
head(frontloaded_LB_in_WB, 10) #shows top 10 genes
frontloaded_LB_in_WB

        #gene:ENSFHEG00000012294 gene:ENSFHEG00000013510 gene:ENSFHEG00000009231 
                       #9.639754                9.469048                9.280234 
        #gene:ENSFHEG00000017456 gene:ENSFHEG00000003293 gene:ENSFHEG00000009823 
                       #8.736212                8.557035                8.519578 
        #gene:ENSFHEG00000012277 gene:ENSFHEG00000003544 gene:ENSFHEG00000023008 
                       #8.101133                7.865141                7.798165 
        #gene:ENSFHEG00000007932 gene:ENSFHEG00000000626 gene:ENSFHEG00000015691 
                       #7.736734                7.582190                7.555056 
        #gene:ENSFHEG00000009196 gene:ENSFHEG00000009514 gene:ENSFHEG00000019924 
                       #7.550657                7.368980                7.288220 
        #gene:ENSFHEG00000009502 
                      #7.161205


#16 genes already on in LB pop and stays at a consistent expression when exposed to WB normal conditions (tide + runoff treatment)


#extract list of frontloaded genes in WB pop
write.csv(frontloaded_LB_in_WB, file = "frontloaded_LB_in_WB.csv")



###Gene Annotation for frontloaded genes


### Long Bayou (induced in WB when exposed to LB control "tide")
gene_ids_WB_in_LB <- sub("^gene:", "", names(frontloaded_WB_in_LB))
frontloaded_df_WB_in_LB <- data.frame(
  gene_id = gene_ids_WB_in_LB,
  expression = as.numeric(frontloaded_WB_in_LB)
)
annotated_WB_in_LB <- merge(frontloaded_df_WB_in_LB, gene_info, by = "gene_id", all.x = TRUE)
write.csv(annotated_WB_in_LB, "annotated_frontloaded_WB_in_LB.csv", row.names = FALSE)


### Weakley Bayou (induced in LB when exposed to WB control "tide_runoff")
gene_ids_LB_in_WB <- sub("^gene:", "", names(frontloaded_LB_in_WB))
frontloaded_df_LB_in_WB <- data.frame(
  gene_id = gene_ids_LB_in_WB,
  expression = as.numeric(frontloaded_LB_in_WB)
)
annotated_LB_in_WB <- merge(frontloaded_df_LB_in_WB, gene_info, by = "gene_id", all.x = TRUE)
write.csv(annotated_LB_in_WB, "annotated_frontloaded_LB_in_WB.csv", row.names = FALSE)


###############Visualize data using heatmap

###HEATMAP###

#VST-normalized expression matrix (normalizing data again by submitting to different output name)
vsd <- vst(GFKdds, blind = FALSE)
vsd_mat <- assay(vsd)


library(pheatmap)

#read in backloading data
genes_WB_in_LBf <- names(frontloaded_WB_in_LB)
genes_LB_in_WBf <- names(frontloaded_LB_in_WB)

#Subset data (VST matrix) for each front-loaded gene list
expr_WB_in_LBf <- vsd_mat[genes_WB_in_LBf, ]
expr_LB_in_WBf <- vsd_mat[genes_LB_in_WBf, ]


# Add annotations (e.g., origin, treatment)
anno <- as.data.frame(colData(vsd)[, c("origin", "treatment")])

### Heatmap for WB genes back-loaded in LB

# Define your custom colors
my_colors <- list(
  origin = c("long_bayou_F" = "blue", "weakley_bayou_U" = "red"),     # blue for LB and red for WB for origins
  treatment = c("tide" = "blue", "tide_runoff" = "red")  # blue corresponds to LB control and red corresponds to WB control for treatments
)

pheatmap(expr_WB_in_LBf,
         annotation_col = anno,
         annotation_colors = my_colors,
         show_rownames = FALSE,
         cluster_cols = FALSE, #keeps samples ordered; I ran both so I have clustering and non-clustered heatmaps
         main = "WB Genes Front-Loaded in LB")

#heatmap saved as "Frontloading_WB_in_LBhome.pdf" and "Frontloading_WB_in_LBhome.png"


### Heatmap for LB genes front-loaded in WB
pheatmap(expr_LB_in_WBf,
         annotation_col = anno,
         annotation_colors = my_colors,
         show_rownames = FALSE,
         cluster_cols = FALSE, #keeps samples ordered; I ran both so I have clustering and non-clustered heatmaps
         main = "LB Genes Front-Loaded in WB")

#heatmap saved as "Frontloading_LB_in_WBhome.pdf" and "Frontloading_LB_in_WBhome.png"


######NOTE: It is normal to get genes in the front-back loading analysis that DOES NOT show up in the pairwise comparison.


#########End of front-back loading analysis#########



