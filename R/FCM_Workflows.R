###############################################################################
# 
# Fuzzy C-means clustering implementation on scRNA-seq data (or pseudo-bulk)
# Figure 4, 5, 6 Trevino et al 2020 (bioRxiv)
# 
# Alexandro E. Trevino 2020
# Stanford University
#
###############################################################################

###############################################################################
# Import Functions
###############################################################################

source('R/Single_Cell_Functions.R')
source('R/FCM_Functions.R')
source('R/Load_Environment.R')

###############################################################################
# 1. Set up the analysis
###############################################################################

# Read our main fuzzy clustering object 
fuzz <- readRDS('data_files/rds/FCM_Object.RDS')

# ~OR~ #

# Create a fresh one from the data:
# you can also compute your own pseudobulks from your own clusters etc.

trial.seed <- 2
set.seed(trial.seed)

pseudo.bulks <- readRDS("data_files/rds/RNA_GlialPseudobulks.RDS")
pbcounts <- pseudo.bulks$matrix
mat <- cleanMatrix(pbcounts, scaling = "row", trim = NULL)


# Set up the project "object"

fuzz <- list()

fuzz[['Trial.Name']]    <- 'descriptive_trial_name'
fuzz[['Random.Seed']]   <- trial.seed
fuzz[['Data.Matrix']]   <- mat
fuzz[['Counts.Matrix']] <- pbcounts
# fuzz[['colData']]       <- glia.colData
fuzz[['Pseudobulks']]   <- pseudo.bulks[c("sampleIdxL", "grouping")]

# You can put the main UMAP or anything else here:
#fuzz[['Projections']][['SCE']][['DF']] <- rna.umap

# Create a dir for this project
dir.create(sprintf("fuzzy_clusters/%s", fuzz[['Trial.Name']]), showWarnings = F)

###############################################################################
# 2. Do fuzzy clustering
###############################################################################

# Define parameters

fuzz <- setFuzzyParameters(
  c = 14,
  m = 1.25,
  object = fuzz
)

# Fuzzy clustering

fuzz <- doFuzzyClustering(fuzz, iter.max = 50)

###############################################################################
# 3. Set threshold, get module overlaps, and get module gene memberships
###############################################################################

# Set overlap threshold

fuzz <- setThreshold(0.2, fuzz) 
print(min(maxMemberhip(fuzz)))

# Plot membership heatmap

plotMembershipHeatmap(fuzz)
plotGenesAtThreshold(fuzz)
plotGeneMembership(fuzz)

# Get thresholded membership in modules - and respective weights from FCM 
# membership matrix

fuzz <- addMemberGeneList(fuzz)

fuzz <- addMemberGeneWeights(fuzz)

# Compute overlaps

fuzz <- computeOverlaps(fuzz, method = 'ratio')
fuzz <- computeOverlaps(fuzz, method = 'number')
fuzz <- computeOverlaps(fuzz, method = 'jaccard')
fuzz <- computeOverlaps(fuzz, method = 'scaled')

plotOverlapPCA(fuzz)

###############################################################################
# 4. Compute GO enrichments 
###############################################################################

fuzz <- addGOEnrichments(fuzz, universe = 'glia')
fuzz <- addGOEnrichments(fuzz)

# can also use BP if you want

fuzz <- addGOEnrichments(fuzz, universe = 'glia', ontology = 'BP')
fuzz <- addGOEnrichments(fuzz, ontology = 'BP')

plotGOfcm(object = fuzz, max.terms = 6, ontology = 'MF')
plotGOfcm(object = fuzz, max.terms = 6, ontology = 'BP')

###############################################################################
# 6. Add pseudotime and plot heatmaps
#
# Calculate PT for each module 
# Get pseudotime from URD (c8, glia only)
#
###############################################################################

fuzz <- computeModulePseudotime(fuzz)

# Sort samples and modules by pseudotime, and plot a heatmap of the centers mat

mod.sample.mat <- cleanMatrix(
  mat = getCenterMatrix(fuzz), 
  scaling = 'none', 
  cuts = c(0.05,0.95)
)

# Metadata & colors for the heatmap

ca <- HeatmapAnnotation(
  df = fuzz$colData[ , 2:4],
  col = list(
    group = glia.color.rna.vector,
    age = time.colors,
    pseudotime = circlize::colorRamp2(
      breaks = seq(0,1, 0.1), 
      colors = viridis_pal(option = 'C')(11)
    ) 
  )
)

# could order by other metadata here

cl.pt.ord <- fuzz$colData %>%
  group_by(group) %>%
  summarize(mean.pseudotime = mean(pseudotime)) %>% 
  arrange(mean.pseudotime) %>%
  mutate(ord.group = factor(group, levels = group)) %>%
  pull(ord.group)

ord.ind <- fuzz$colData %>%
  mutate(group = factor(group, levels = cl.pt.ord),
         ord = seq_len(n())) %>%
  arrange(group, pseudotime) %>%
  pull(ord)

# just automate
mod.order <- c(10,3,8,12,6,9,11,5,14,2,13,4,1,7)

pdf(
  file = 'fuzzy_clusters/Module_by_Sample_Heatmap.pdf', 
  height = 4.5, 
  width = 6.5
)

Heatmap(
  matrix = mod.sample.mat, 
  row_order = mod.order, 
  column_order = ord.ind, 
  width = unit(7, 'cm'), height = unit(7, 'cm'),
  top_annotation = ca,
  use_raster = T
)
dev.off()

# Plot select genes in the heatmap

select.genes <- c('top2a', 'nfic', 'nr2f1', 'lmo2', 'foxj1', 'aqp4', 'MBP') %>% 
  toupper()

sel.gene.ids <- select.genes %>% convertGeneIDs(., name2ensembl)

hm.name <- sprintf('fuzzy_clusters/%s/Select_Genes_Sample_Heatmap.pdf', fuzz$Trial.Name)
pdf(
  file = hm.name,
  height = 3,
  width = 6.5
)
Heatmap(
  matrix = trimQuantiles(fuzz$Data.Matrix[sel.gene.ids, ]), 
  cluster_rows = F, 
  row_labels = select.genes,
  column_order = ord.ind,
  width = unit(7, 'cm'), 
  height = unit(3, 'cm'),
  use_raster = T,
  col = brewer.pal(9,'Reds')
)
dev.off()

# Gene version of heatmap

gene.order <- do.call(c, getMemberGeneList(fuzz)[order(fuzz$Module.Pseudotime)])
names(gene.order) <- NULL

id.order <- convertGeneIDs(genes.order, name2ensembl)
module.df <- getMemberGeneDF(fuzz)

genemod.sample.mat <- cleanMatrix(
  mat = getDataMatrix(fuzz)[id.order, ],
  scaling = 'none',
  cuts = c(0.05, 0.95)
)

pdf(
  file = sprintf('fuzzy_clusters/%s/Genes_and_Modules_by_Sample_Heatmap.pdf', fuzz$Trial.Name),
  height = 6.5, 
  width = 6.5
)
Heatmap(
  matrix = genemod.sample.mat, 
  row_order = id.order,
  show_row_names = F,
  split = module.df[,1],
  column_order = order(fuzz$colData$pseudotime), 
  width = unit(7, 'cm'), height = unit(12, 'cm'),
  top_annotation = ca,
  use_raster = T
)
dev.off()

# saveFCM(fuzz)

###############################################################################
# 7. Re-project pseudobulks from FCM space
###############################################################################

# UMAP (optimization elsewhere)
my.umap.config <- umap::umap.defaults
my.umap.config[["min_dist"]] <- 0.4

# Adds 'UMAP' and 'Projections'[['UMAP']]
fuzz <- fcmUMAP(fuzz, proj.name = 'UMAP', my.umap.config)

###############################################################################
# 8. Re-cluster projected data (fuzz$centers). 
#     With Seurat code assist from Fabian
###############################################################################

resolutions <- c(0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.2)

fuzz <- fcmProjectionReclustering(
  object = fuzz, 
  proj.name = 'UMAP',
  resolutions = resolutions, 
  find.centroids = F
)

fuzz <- setReclusterResolution(
  object = fuzz, 
  proj.name = 'UMAP',
  resolution = 'cr_0.8'
)

fuzz <- fcmProjectionReclustering(fuzz, find.centroids = T)

fuzz <- addFCMProjectionDataFrame(
  object = fuzz, 
  proj.name = 'UMAP', 
  metadata = fuzz$colData
)

# Plot all the covariates till now

p1 <- plotProjection(
  object = fuzz, 
  proj.name = 'UMAP',
  covariate = 'group', 
  plot.title = 'Original clusters', 
  color.palette = glial.colors$color
)
p2 <- plotProjection(
  object = fuzz, 
  covariate = 'age', 
  plot.title = 'Sample age', 
  color.palette = time.colors
)
p3 <- plotProjection(
  object = fuzz, 
  covariate = 'pseudotime', 
  plot.title = 'Pseudotime', 
  color.palette = viridis_pal(option = 'C')(100)
)
p4 <- plotProjection(
  object = fuzz, 
  covariate = 'clusters', 
  plot.title = 'Louvain reclustering', 
  color.palette = stata_pal()(fuzz$Projections$UMAP$Reclustering$nclust),
  alpha = 0.3,
  add.layers = fuzz$Projections$UMAP$Reclustering$Centroids,
  add.geoms = c(
    'geom_point(data = add.layers, shape = 8, size = 3, color = \'black\')',
    'geom_label_repel(data = add.layers, aes(label = clusters), size = 2, box.padding = .05, show.legend = F)',
    'guides(color = guide_legend(override.aes = list(alpha = 1)))'
    )
) 

patch <- (p1 + p2) / (p3 + p4)
pdf(
  file = sprintf("fuzzy_clusters/%s/UMAP.pdf", fuzz[['Trial.Name']]), 
  width = 9,
  height = 9,
  useDingbats = F
)
print(patch)
dev.off()

###############################################################################
# 9. Project ***module*** centroids:
###############################################################################

fuzz <- computeModuleExpressionCentroids(fuzz, threshold = .1)

fuzz <- drawModuleConnectivity(fuzz, overlap.method = 'ratio')
fuzz <- drawModuleConnectivity(fuzz, overlap.method = 'number')
fuzz <- drawModuleConnectivity(fuzz, overlap.method = 'jaccard')
fuzz <- drawModuleConnectivity(fuzz, overlap.method = 'filter')

fuzz <- weightModuleConnectivity(fuzz, 'ratio', threshold = .08)
fuzz <- weightModuleConnectivity(fuzz, 'number', threshold = 80)
fuzz <- weightModuleConnectivity(fuzz, 'jaccard', threshold = .2)
fuzz <- weightModuleConnectivity(fuzz, 'filter', threshold = .1)

plotModuleConnectivity(fuzz, 'ratio')
plotModuleConnectivity(fuzz, 'number')
plotModuleConnectivity(fuzz, 'jaccard')
plotModuleConnectivity(fuzz, 'filter')

# Elbow test to find a cutoff / threshold. Jaccard probably the most sensible
# metric here

connections <- lapply(seq(0,1,.025), function(x) { 
  sum(fuzz$Overlaps$jaccard > x) 
  }) %>% Reduce('c', .)

pdf(
  file = sprintf("fuzzy_clusters/%s/Jaccard_Threshold_Elbow_Test.pdf", trial.name),
  width = 4.5, 
  height = 4.5, 
  useDingbats = F
)
plot(seq(0,1,.025), connections)
abline(v = 0.2, col = 'red')
dev.off()

###############################################################################
# 10. Plot module expression, in FCM space (see how modules are expressed)
###############################################################################

plotModuleExpression(
  object = fuzz, 
  plot.object = fuzz[['Data.Matrix']], 
  plot.object.class = 'matrix', 
  proj.name = 'UMAP', 
  dir.name = NULL, 
  weight.expression = T, 
  scale.expression = T,
  convert = T,
  point.size = 2,
  aspectratio = 1
)

# in the Main UMAP now

plotModuleExpression(
  object = fuzz, 
  plot.object = sce, 
  plot.object.class = 'sce', 
  proj.name = 'SCE', 
  dir.name = NULL, 
  weight.expression = T, 
  scale.expression = T
)

###############################################################################
# 11. Project other data into the FCM space
###############################################################################

# Wrangle data matrices and metadata

# ATAC PB metadata

atac.pb.colData <- data.frame(
  pseudotime = pb.pseudotimes, 
  group = factor(glial.groupings.atac, levels = glial.colors.atac$cluster)
)
fuzz[['pb.ATAC.colData']] <- atac.pb.colData

# Fix matrices of ATAC-seq data, which don't quite match features. 
# --> missing some lincRNA from gene activity calculation. 
atac.matrix <- imputeMissingFeatures(
  object = fuzz, 
  newdata.matrix = glial.pb.ga.mat, 
  newdata.name = 'Gene.Activity.Matrix', 
  convert.ids = T,
  impute.fun = 'median'
)

# Project the data

fuzz <- projectCells(
  object = fuzz,
  proj.name = 'ATAC.UMAP',
  cell.matrix = atac.matrix,
  metadata = atac.pb.colData,
  add.matrix = T
)

# Add projection dataframes to plot these new projections
a.p1 <- plotProjection(
  object = fuzz, 
  proj.name = 'ATAC.UMAP',
  covariate = 'group',
  plot.title = 'ATAC-seq clusters -- pseudobulk gene score NNs', 
  color.palette = glial.colors.atac$color
)
 
pdf(
  file = sprintf('fuzzy_clusters/%s/ATAC_UMAP.pdf', fuzz[['Trial.Name']]),
  width = 9, height = 5 ,
  useDingbats = F
)

print(a.p1)

dev.off()

# Plot these new projections

# Also, add linked peaks to each module - Regulatory information

fuzz <- addLinkedATACModules(
  object = fuzz, 
  link.object = sig.links, 
  bg.link.object = all.links,
  make.unique = F
)

# Plot peak accessibility in the main UMAP

fcm.atac.dir <- sprintf(
  "fuzzy_clusters/%s/ModuleAccessibility_thres%s/", 
  fuzz$Trial.Name, 
  as.character(fuzz$Threshold)
)

dir.create(fcm.atac.dir, showWarnings = F)
fcm.atac <- fuzz$Module.Linked.Peaks
fcm.background.atac <- fuzz$Module.Background.Linked.Peaks

lapply(
  seq_along(fcm.atac),
  function(X) {
    pdf(sprintf("%s/Module_%s.pdf", fcm.atac.dir, names(fcm.atac)[X]))
    plotAccessibility(
      peaks = fcm.atac[[X]],
      background.peaks = fcm.background.atac[[X]]
    )
    dev.off()
  }
)

# see plotExpression; plotActivity; plotCvar to plot specific genes / loci / 
# motifs with a "fuzz" object (or any object)

###############################################################################
# 12. Plot motif enrichments (in linked peaks) across modules
###############################################################################

# Construct a base plot as a backdrop - to plot this other data on the centroids

segments <- fuzz$Overlap.Segments.Filtered$jaccard
base.df <- fuzz$Projections$UMAP$DF

base.gg <- ggplot(base.df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(size = 2, shape = 16, color = 'grey') +
  guides(color = guide_legend(override.aes = list(alpha = 1)),
         alpha = guide_legend(title = str_to_title('Jaccard index'))) +
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2),
               data = segments, 
               inherit.aes = F, 
               color = 'black') +
  theme_classic() +
  theme(axis.line = element_blank(),
        panel.border = element_rect(fill = NA)) +
  theme(aspect.ratio = 1)

#

fcm <- getMemberGeneList(fuzz)

# Find linked peaks within the modules

mod.peaks <- do.call(
  rbind,
  lapply(seq_along(fcm), function(x) {
    mname <- paste0('m', as.character(x))
    genes <- fcm[[x]]
    pks <- findLinkedPeaks(genes, object = sig.links, make.unique = T)
    out <- data.frame(
      cbind(
        module = mname,
        peaks = pks
      )
    )
  })
)

# Find linked peaks *to top ranked genes* within the modules

num.genes.to.rank <- 50

mod.peaks <- do.call(
  rbind,
  lapply(seq_along(fuzz$Member.Gene.Weights), function(x) {
    
    mname <- paste0('m', as.character(x))
    weights <- sort(fuzz$Member.Gene.Weights[[x]], decreasing = T)[1:num.genes.to.rank]
    genes <- convertGeneIDs(names(weights), ensembl2name)
    
    pks <- findLinkedPeaks(genes, object = sig.links, make.unique = T)
    
    out <- data.frame(
      cbind(
        module = mname,
        peaks = pks
      )
    )
    
  })
)


# Compute enrichments across these clusters

mod.motif.enr <- hypergeometricTestList(
  groupings = mod.peaks$module, 
  match_matrix = motif.mat[mod.peaks$peaks, ], 
  lower_tail = 'both'
) %>% mutate(p.value = pmin(enrichment.p, depletion.p)/2,
             p.adj = p.adjust(p.value, method = 'bonf'),
             neg.log10.p = -log10(p.adj),
             log2FoldEnrichment = log2((matches / draws) / (bg.matches / bg.non.matches)),
             motif.name = add_split_name_vector(names, "_", 2))

match.mod.motif.enr <- mod.motif.enr %>%
  group_by(group) %>%
  slice_min(enrichment.p, n = 3) %>%
  inner_join(tf.match.possibilities, by = c("motif.name", "names" = "motif.id")) %>%
  group_by(group, motif.name) %>%
  slice_max(glia.pb.rho) %>%
  dplyr::select(group, motif.name, gene.symbol, glia.pb.rho, log2FoldEnrichment, enrichment.p, p.adj, neg.log10.p) 

memmotif.dir <- sprintf('fuzzy_clusters/%s/Motif_Module_Plots/', getTrialName(fuzz))
dir.create(memmotif.dir, showWarnings = F)

# Selected motifs of interest here

sel.motif.ids <- c('ID4', 'NFIA', 'ASCL1', 'EGR1', 'FOS', 'HES5', 'OLIG2', 
                   'MEIS2', 'NHLH1', 'SOX21', 'SOX10', 'PBX1', 'NEUROD1', 
                   'NEUROG2', 'EOMES')

motif.dfs <- lapply(sel.motif.ids, function(x) {
  df <- motifModuleDF(x)
  out <- cbind(fuzz$Module.Centroids, df)
  }
)

pdf(
  file = sprintf('fuzzy_clusters/%s/Motif_Module_Plots/Structured.pdf', fuzz$Trial.Name),
  width = 5.5,
  height = 4,
  useDingbats = F
)

lapply(
  seq_along(motif.dfs),
  function(x) {
    gg <- base.gg +
      geom_point(size = 7, data = motif.dfs[[x]], shape = 21, aes(fill = neg.log10.p)) +
      # geom_blank(data = motif.dfs[[x]], aes(fill = -neg.log10.)) +
      geom_text(size = 3, data = motif.dfs[[x]], aes(label = group)) +
      scale_fill_gradientn(colors = rev(brewer.pal(9, 'BrBG'))) +
      labs(title = sel.motif.ids[x]) 
      # xlim(c(-8,-1)) +
      # ylim(c(4,11))
  }
)

dev.off()


# Plot enrichments 

plotMotifsInModules('ASCL1')
plotMotifsInModules('OLIG2')
plotMotifsInModules('HES5')
plotMotifsInModules('EGR1')
plotMotifsInModules('FOS')
plotMotifsInModules('NFIA')
plotMotifsInModules('ID4')





