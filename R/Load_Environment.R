###############################################################################
#
# Set up global R environment for single cell analysis of human corticogenesis
#
# Alexandro E. Trevino, Stanford University 2020
#
###############################################################################

###############################################################################
# Import packages
###############################################################################

# Color packages

library(viridis)
library(wesanderson)
library(RColorBrewer)

# Plotting packages

library(circlize)
library(ComplexHeatmap)
library(patchwork)
library(gplots)
library(grid)
library(ggplot2)
library(ggrepel)
library(ggrastr)
library(ggdendro)
library(ggthemes)
library(gghighlight)

# Tidy data packages

library(dplyr)
library(magrittr)
library(reshape2)
library(splitstackshape)
library(broom)
library(stringr)
library(binr)
library(tibbletime)

# Statistics packages

library(e1071)
library(fclust)

# Bioinformatics packages

library(org.Hs.eg.db)
library(clusterProfiler)
library(SingleCellExperiment)
library(GenomicInteractions)
library(GenomicRanges)
library(motifmatchr)
library(JASPAR2018)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(URD)
library(Seurat)

###############################################################################
# Set working directory
###############################################################################

# Make it the cloned git directory
setwd("brainchromatin")

###############################################################################
# Import functions
###############################################################################

source("R/Single_Cell_Functions.R")
source("R/FCM_Functions.R")

###############################################################################
# Define any visualization themes and color palettes
###############################################################################

# My custom color palettes

pals <- read.table(
  file = "misc/colorpals.ai.csv", 
  header = T , 
  sep = ",", 
  fill = T, 
  comment.char = ""
)

# Define and name

plant_palettes <- lapply(
  levels(pals$palette.name), 
  function(x) {
    
    pal = pals %>% 
      filter(palette.name == x) %>% 
      pull(hex) %>%
      as.character()
    
    return(pal)
  }
)
names(plant_palettes) <- levels(pals$palette.name)

# Cluster colors for scRNA and scATAC datasets

color.rna <- read.table(
  file = "data_files/tsv/scRNA_Colors.txt", 
  header = T, 
  sep = "\t", 
  fill = NA, 
  comment.char = "", 
  stringsAsFactors = F
)

color.atac <- read.table(
  file = "data_files/tsv/scATAC_Colors.txt", 
  header = T, 
  sep = "\t", 
  fill = NA, 
  comment.char = "", 
  stringsAsFactors = F
)

# Vector format

color.rna.vector <- color.rna$color
names(color.rna.vector) <- color.rna$cluster

rownames(color.atac) <- color.atac$cluster

color.atac.vector <- color.atac$color
names(color.atac.vector) <- rownames(color.atac)

# Sample age colors 

time.colors <- c('#D7496E', '#FAAD6C', '#D84C27', '#7E2B18')
names(time.colors) <- c('pcw16', 'pcw20', 'pcw21', 'pcw24')

###############################################################################
# Import scRNA-seq data
###############################################################################

# Processed SingleCellExperiments (expression) 
sce <- readRDS('data_files/rds/scRNA_SCE.RDS') 

# scRNA-seq metadata
colData <- colData(sce) %>% as.data.frame()

# scRNA-seq UMAP
rna.umap <- data.frame(sce@int_colData$reducedDims$UMAP)

# Marker genes (Seurat)
markers <- readRDS('data_files/rds/scRNA_Cluster_Markers.RDS')

# Basic unfiltered gene annotations
gene.gr <- readRDS('data_files/rds/AllGenes_GenomicRanges.RDS')

###############################################################################
# Some useful genes
###############################################################################

gs.cellcycle <- read.table(
  file = "gene_sets/gs_KEGG_CELL_CYCLE.txt", 
  skip = 2, 
  stringsAsFactors = F
)[[1]]

sfari.db <- read.table(
  file = "misc/sfari_gene.csv",
  header = T, 
  sep = ",",
  stringsAsFactors = F
)

# Table of TFs
tfs <- read.table(
  file = "gene_sets/geneset_TF_all.txt", 
  header = F, 
  skip = 2, 
  sep = "\t", 
  stringsAsFactors = F
) %>% pull(V1)

# Quickly find any TFs
transcriptionFactors(c('EOMES', 'GAPDH', 'VIM', 'NEUROG2'))

# Flexible or quick plotting functions, e.g.:
geneSetAveragePlot(gs.cellcycle, plot.type = 'averages')
plotrna(c('hopx', 'cryab'))

###############################################################################
# Gene ID mapping
###############################################################################

# Mapping gene IDs and gene symbols (for easier use of gene symbols)
# This annotation file contains metadata about each gene, including the Ensembl 
# ID and the gene symbol

name2ensembl <- names(gene.gr)

# Something funny about this one, unsure why yet
name2ensembl['CDR1'] <- 'ENSG00000281508'

names(name2ensembl) <- elementMetadata(gene.gr)[ ,'gene_name']

ensembl2name <- elementMetadata(gene.gr)[ ,'gene_name']
names(ensembl2name) <- names(gene.gr)

# A quick function to interconvert, e.g.:
convertGeneIDs('NEUROD2', name2ensembl)

###############################################################################
# Import scATAC-seq data
###############################################################################

# SCE object
atac.sce <- readRDS('data_files/rds/scATAC_SCE.RDS')

# Cell metadata
atac.colData <- colData(atac.sce) %>% as.data.frame()

# Integration mapping to RNA
cca <- readRDS('data_files/rds/CCA_Matching.RDS')

# ChromVAR results
cvar <- readRDS('data_files/rds/scATAC_ChromVAR_SCE.RDS')

# Gene activity
ga <- readRDS('data_files/rds/scATAC_GeneActivity.RDS')

# Motif matrix
motif.mat <- readRDS('data_files/rds/scATAC_MotifMatchMatrix.RDS')

# Peak gene linkages

sig.links <- readRDS('data_files/rds/PeakGeneLinks_Significant.RDS')
all.links <- readRDS('data_files/rds/PeakGeneLinks_All.RDS')

###############################################################################
# Read the CCA data (currently ChrAccR gene scores)
###############################################################################

# get a vector that maps accessibility cells to RNA cells
acc2expr <- cca %>% filter(mapping == "nearest_expr_for_acc")
map.acc2expr <- acc2expr$cell_expr
names(map.acc2expr) <- acc2expr$cell_acc

# and vice versa
expr2acc <- cca %>% filter(mapping == "nearest_acc_for_expr")
map.expr2acc <- expr2acc$cell_acc
names(map.expr2acc) <- expr2acc$cell_expr

###############################################################################
# Glia / Fuzzy Clustering 
###############################################################################

# URD object to flood / pseudotime stage + var genes
# Nice functions in URD to get variable genes and pseudotime. 
glial.urd <- readRDS('data_files/rds/scRNA_URD_GlialCells.RDS')

# 
var.glial.gene.ids <- glial.urd@var.genes
var.glial.genes <- convertGeneIDs(var.glial.gene.ids, ensembl2name)
pseudotime <- glial.urd@pseudotime$pseudotime # associated PT

# Sampled pseudobulks (From KNN graph) for RNA:
#   a. read
spb.rna.glia.nn <- readRDS('data_files/rds/RNA_GlialPseudobulks.RDS')
glial.pb.mat <- spb.rna.glia.nn$matrix

#   b. Function to scale, remove NA, etc with matrices
clean.glia.pb <- cleanMatrix(
  glial.pb.mat, 
  scaling = "row", 
  cuts = NULL
)

#   c. Subtract cell cycle genes (for fuzzy clusters)
no.cc <- (!var.glial.genes %in% gs.cellcycle)

#   d. define metadata for pseudobulk aggregates
glia.colData <- glial.urd@meta

glial.groupings <- lapply(
    X = spb.rna.glia.nn$sampleIdxL,
    FUN = dominantCluster, 
    groupings = glial.urd@meta$seurat_clusters
) %>% Reduce('c', .)

glial.ages <- lapply(
    X = spb.rna.glia.nn$sampleIdxL,
    FUN = dominantCluster, 
      groupings = glial.urd@meta$Age
) %>% Reduce('c', .)

glial.ages <- factor(glial.ages, levels = names(time.colors))
glial.colors <- color.rna %>% filter(cluster %in% glial.groupings)

# This matrix just includes all genes - not just variable-in-glia genes
# For more plotting flexibility
full.glial.pb.mat <- readRDS('data_files/rds/Matrix_Glial_KNNpseudoBulk_RNA_Full.RDS')

###############################################################################
# And ATAC pseudobulks 
###############################################################################

spb.cca.atac <- readRDS('data_files/rds/SamplePseudoBulk_Glial_KNNpseudobulk_ATAC_GA.RDS')
glial.pb.ga.mat <- readRDS('data_files/rds/Matrix_Glial_KNNpseudoBulk_ATAC_GeneActivity.RDS')
glial.pb.atac.mat <- readRDS('data_files/rds/Matrix_Glial_KNNpseudoBulk_ATAC_Accessibility.RDS')
glial.pb.cvar.mat <- readRDS('data_files/rds/Matrix_Glial_KNNpseudoBulk_ATAC_ChromVAR.RDS')

atac.nn.glia <- findNearestNeighbor(rownames(glial.urd@meta), nn.mapping = map.expr2acc)
pre.grouping <- atac.colData[atac.nn.glia, ][['Iterative.LSI.Clusters']]
names(pre.grouping) <- atac.nn.glia

glial.groupings.atac <- lapply(
    X = spb.cca.atac,
    FUN = dominantCluster, 
    groupings = pre.grouping
) %>% Reduce('c', .)

glial.colors.atac <- color.atac %>% filter(cluster %in% glial.groupings.atac)


###############################################################################
# TF match possibilities and correlations table
###############################################################################

tf.match.possibilities <- readRDS('data_files/rds/TF_MotifExpressionCorrelation.RDS')

###############################################################################
###############################################################################
###############################################################################




