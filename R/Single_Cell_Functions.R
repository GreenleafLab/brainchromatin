###############################################################################
#
#   Single Cell ATAC and RNA 
#   Global environment functions
#
#     -  Plotting functions
#     -  Bioinformatics functions
#     -  Utility functions
#
#   Written OR adapted by Alexandro E. Trevino 2017-2020 
#   Stanford University
#   Thanks to Fabian Mueller (FM), the source of several key functions
# 
###############################################################################

###############################################################################
#
# Plotting functions
#
# #
# # # # # #
# # # # # # # # # #
# # # # # # # # # # # # # # # # # # # #
###############################################################################

#' Produces a two-color highlight of selected cells
#' 
#' Currently this just highlights a group of pre-defined cells. In the future 
#' it should also flexibly plot data or metadata in a defined group of 
#' highlighted cells. 
#' 
#' @param cells A vector of cell IDs corresponding to rownames in dim.reduction
#' @param plot.title A character string
#' @param dim.reduction E.g. a UMAP dataframe. Rownames should be cell IDs
#' @param highlight.color A character string specifying a color 
#' @param aspectratio Plot aspect ratio
#' @param point.size The plot point point size, passed to ggplot2
highlightCells <- function(cells, 
                           plot.title = "",
                           dim.reduction = rna.umap,
                           highlight.color = 'black',
                           aspectratio = 1.25,
                           point.size = 0.3) {
  
  pl.df <- dim.reduction %>%
    mutate(Highlight = ifelse(rownames(dim.reduction) %in% cells, T, F)) 
  
  dim1 <- grep('UMAP', colnames(pl.df), value = T)[1]
  dim2 <- grep('UMAP', colnames(pl.df), value = T)[2]
  
  gg <- ggplot(pl.df, aes_string(x = dim1, y = dim2, color = 'Highlight')) +
    geom_point(size = point.size, shape = 16) +
    theme_classic() +
    labs(title = plot.title) +
    scale_color_manual(values = c("grey67", highlight.color)) +
    theme(aspect.ratio = aspectratio,
          plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          axis.line.x = element_blank(),
          axis.line.y = element_blank(),
          legend.position = 'none')
  
  print(gg)
  # todo: separate UMAP plotting into another more flexible function
}

#' Function to plot gene expression from a SingleCellExperiment or Seurat 
#' object with a dimensionality reduction. Easier for me and more flexible than
#' Seurat. Uses ggplot2.
#' 
#' Takes a vector of genes and a SCE object, as well as other parameters, and 
#' produces an expression plot. Option to make a PDF directly.
#' 
#' @param genes a vector of genes
#' @param object a SingleCellExperiment object (maybe Seurat too)
#' @param object.class the type of object. "sce" or "seurat", default "sce"
#' @param projection if your object is just a matrix, this is the coordinate 
#' projection
#' @param assay.name In Seurat object, what is the Assay being used.
#' @param embedding In Seurat object, what is the Embedding name
#' @param convert Whether or not to convert gene symbols to ENSEMBL IDs first
#' @param plot.title Plot title
#' @param guide.name Legend title (ggplot2)
#' @param plot.type Whether you are plotting 1 gene, the average of the genes 
#' listed, or individual panels for each of the genes listed (currently uses 
#' ggplot2 facets, cowplot grids coming)
#' @param w A vector of weights for gene expression averages. For example, 
#' cluster membership scores. 
#' Automatically computes a weighted mean if not NULL.
#' @param scaled Whether or not to scale gene expression across cells
#' @param color.palette A vector of colors, passed to scale_color_gradientn()
#' @param aspectratio The desired aspect ratio of each plot/panel
#' @param rastr Whether or not to use ggrastr to plot (plays poorly with 
#' aspectratio and with server RStudio due to Cairo issues. Meh.)
#' @param point.size custom point size
#' @param trim For plotting large datasets, trims the uppermost and lowermost
#' percentiles, eg c(0.01, 0.99) trims the highest and lowest percentile, so
#' the color scale is not driven by outlier points
#' @param lims X / Y limits of coordinates. subset your projection space
#' @param print Whether the plot is passed through "print"
#' @param num.panel.rows For facet plots, how many rows of facets?
geneSetAveragePlot <- function(genes, 
                               object = sce, 
                               object.class = "sce",
                               projection = NULL,
                               assay.name = "RNA",
                               embedding = 'umap',
                               convert = T,
                               plot.title = "Gene Expression", 
                               guide.name = "Gene expression",
                               plot.type = c("single", "averages", "panels"), 
                               w = NULL,
                               scaled = F,
                               color.palette = plant_palettes[['Iris']], 
                               aspectratio = 1.25,
                               rastr = F,
                               point.size = 0.3,
                               trim = NULL,
                               lims = NULL,
                               print = T,
                               num.panel.rows = NULL) {
  
  # exception in gene names...
  
  genes <- validGenes(genes)
  
  # Gene ID input handling
  if(convert == T) {
    gene.ids <- convertGeneIDs(genes, id.mapping.table = name2ensembl)
  } else {
    gene.ids <- genes
  }
  
  filt.ind <- filterGenesNotInData(
    genes = gene.ids, 
    object = object, 
    to.return = 'indices'
  )
  gene.ids <- gene.ids[filt.ind]
  genes <- genes[filt.ind]
  
  # Object class handling
  
  object.class <- tolower(object.class)
  
  if (object.class == "seurat") {
    rna.mat <- logcounts(as.SingleCellExperiment(object))
    umap <- Seurat::Embeddings(object, embedding)
  }
  if (object.class == "sce") {
    rna.mat <- logcounts(object)
    umap <- SingleCellExperiment::reducedDim(sce, "UMAP")
  }
  if (object.class == 'matrix') {
    rna.mat <- object
    if (is.null(projection)) {
      stop("Must provide a projection with matrix class")
    }
    umap <- projection
  }

  # Check again
  if (!all(gene.ids %in% rownames(rna.mat))) {
    print(gene.ids[!gene.ids %in% rownames(rna.mat)])
  }
  
  # Make submatrix
  gene.exp <- rna.mat[gene.ids, ]
  dims <- grep('UMAP', colnames(umap), value = T)
  dim1 <- dims[1]
  dim2 <- dims[2]
  
  umap <- umap[ , c(dim1, dim2)]

  # Scale
  if (scaled == T) {
    gene.exp <- rowScale(gene.exp)
  }
  
  # Trim
  if (! is.null(trim)) {
    gene.exp <- trimQuantiles(gene.exp, cuts = trim)
  }

  # Handle plot types. Make plot data frames
  plot.type = match.arg(plot.type)
  
  if(plot.type == "single" || (length(gene.ids) == 1)) {
    pl.df <- data.frame(umap, geneexpr = gene.exp)
    plot.title <- genes
  } 
  
  if (plot.type == "averages" & (length(gene.ids) > 1)) {
    
    if (is.null(w)) {
      gene.exp <- colMeans(gene.exp)
    } else {
      w <- w[filt.ind]
      gene.exp <- apply(
        X = gene.exp, 
        MARGIN = 2,
        FUN = weighted.mean,
        w = w
      )
    }
    
    pl.df = data.frame(umap, geneexpr = gene.exp)

  }
  
  if(plot.type == "panels" & length(gene.ids) >1 ) {
    
    gene.exp = t(as.matrix(gene.exp))
    colnames(gene.exp) = genes
    
    cat(sprintf("%s genes; %s cells\n", dim(gene.exp)[2], dim(gene.exp)[1]))
    
    pl.df2 = as.data.frame(cbind(umap, gene.exp))
    pl.df = melt(pl.df2, id.vars = c(dim1, dim2), variable.name = "Gene", value.name = "geneexpr")
    
    if (is.null(num.panel.rows) ) {
      num.panel.rows <- pmax(1, floor(length(genes)/3))
    }

  }

  # Plot
  gg = ggplot(pl.df, aes_string(x = dim1, y = dim2, color = 'geneexpr')) +
    theme_classic() +
    scale_color_gradientn(colours = color.palette, name = guide.name) +
    labs(title = plot.title) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(colour = "black", fill=NA, size=1))
  
  if (rastr == F) {
    gg = gg + 
      geom_point(size = point.size, shape = 16) +
      theme(aspect.ratio = aspectratio)
  } else {
    gg = gg + ggrastr::geom_point_rast()
  }
  
  if (plot.type == "panels") {
    gg = gg + facet_wrap(~Gene, nrow = num.panel.rows)
  }
  
  if (!is.null(lims)) {
    x.lims = lims[[1]]
    y.lims = lims[[2]]
    gg = gg + xlim(x.lims) + ylim(y.lims)
  }
  
  if (print == T) { 
    print(gg)
  } else {
    return(gg)
  }
  
}

#' Convenient plotting for a small set of genes
#' 
#' Wrapper around geneSetAveragePlot
#' @param genes Genes to plot
plotrna <- function(genes) {
  
  if (length(genes) > 1 ) { 
    ptype = 'panels' 
  } else {
    ptype = 'single'
  }
  return(geneSetAveragePlot(genes = genes, plot.type = ptype, trim = c(0.01,.99)))

}

#' Plot raw or normalized *aggregate accessibility in a projection space. 
#' 
#' Takes a peakset (must be more than just one) and optionally, a background
#' peakset (there is some batch effect present in single cell counts), and 
#' plots average accessibility across those peaks for each single cell
#' 
#' You could also easily use geneSetAveragePlot (above) for this...
#' 
#' @param peaks A vector of peak names
#' @param background.peaks A vector of peak names. Accessibility will be 
#' normalized against this background accessibility
#' @param object The ATAC-seq object or matrix
#' @param plot.type Whether to average or sample or what (not implemented)
#' @param projection The projection. Should line up with the object above
#' @param color.palette A vector of colors
#' @param plot.title A character string
#' @param aspectratio The aspect ratio of the final plot
#' @param point.size Custom point size
plotAccessibility <- function(peaks, 
                              background.peaks = NULL,
                              object = atac.sce,
                              plot.type = 'averages',
                              projection = atac.umap,
                              color.palette = plant_palettes$Unnamed,
                              plot.title = "Accessibility",
                              aspectratio = 1.25,
                              point.size = .3) {
  
  
  if(length(peaks) < 2) {
    stop("Must provide more than one peak")
  }
  
  ind <- object@rowRanges$name %in% peaks
  
  accessibility <- colSums(assay(object)[ind, ])
  
  if (! is.null(background.peaks)) {
    
    bg.ind <- object@rowRanges$name %in% background.peaks
    bg.accessibility <- colSums(assay(object)[bg.ind, ])
    accessibility <- accessibility / bg.accessibility
  }
  
  pl.df <- data.frame(projection, Accessibility = accessibility)
  
  gg <- ggplot(pl.df, aes(x = UMAP1, y = UMAP2, color = Accessibility)) +
    theme_classic() +
    scale_color_gradientn(colours = color.palette, name = "Accessibility") +
    labs(title = plot.title) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          aspect.ratio = aspectratio) +
    geom_point(size = point.size, shape = 16)
  
  print(gg)
}


#' Plot ChromVAR deviations in a umap
#' 
#' Takes a motif name and plots the deviation Z-scores
#' 
#' @param motif.name Name of the motif. It's passed to grep so know your motifs
#' and be unique
#' @param plot.title Plot title
#' @param file.name Deprecated
#' @param cvar.sce A SingleCellExperiment containing the chromVAR dev Z scores
#' @param umap The embedding to use
#' @param aspectratio plot aspect ratio
#' @param width deprecated
#' @param height deprecated
#' @param color.palette Color scheme
plotChromVAR <- function(motif.name, 
                         plot.title = paste0(motif.name, " ChromVAR"), 
                         file.name = paste0(paste(motif.name, sep = "_"), ".pdf"), 
                         cvar.sce = cvar, 
                         umap = atac.umap, 
                         aspectratio = 1.25, 
                         width = 6, 
                         height = 9,
                         color.palette = plant_palettes$`Photinia fraseri B`) {
  
  # plot.type = match.arg(plot.type)
  
  name = grep(motif.name, rownames(cvar.sce@assays$data[["z"]]), value = T)[1]
  devs = cbind(cvar.sce@assays$data[["z"]][name, ])
  colnames(devs) = name
  
  pl.df2 = as.data.frame(cbind(umap, devs))
  pl.df = melt(pl.df2, id.vars = c("UMAP1", "UMAP2"), variable.name = "Motif", value.name = "Deviation.Z.scores")
  
  gg <- ggplot(pl.df, aes(x = UMAP1, y = UMAP2, color = Deviation.Z.scores)) +
    geom_point(size = .3, shape = 16) +
    geom_blank(aes(color = -Deviation.Z.scores)) +
    theme_classic() +
    scale_color_gradientn(colours = color.palette, name = "Deviation\nZ-score") +
    labs(title = plot.title) +
    theme(aspect.ratio = aspectratio,
          plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(colour = "black", fill=NA, size=1), 
          axis.line.x = element_blank(), 
          axis.line.y = element_blank())
  
  pdf(file.name, width = width, height = height, useDingbats = F);  print(gg)  ; dev.off()
  
}

#' Volcano plot of motif enrichments
#' 
#' Format takes a data frame with columns for enrichment and p value
#' @param dat.subset An identifier that maps to the 'grouping.var'. Subsets DF.
#' @param grouping.var DF column of data subgroups. Can be used to pass through
#' e.g. lapply, and print separate plot images; or could be used to facet the
#' data
#' @param test.df A dataframe with the plotting variables. 
#' @param enrichment.var DF column that stores fold enrichments (x axis)
#' @param p.value.var DF column that stores p-values (y axis)
#' @param do.neg.log10 T or F, whether to log-transform P values prior to 
#' plotting.
#' @param dep.n Integer indicating how many top depletion results to highlight
#' @param enr.n As dep.n, but for enrichment.
#' @param facet Whether or not to employ ggplot faceting
#' @param title.string A string for the plot title.
volcanoPlot <- function(dat.subset = NULL,
                        grouping.var = "group",
                        test.df, 
                        enrichment.var = "log2.enrichment",
                        p.value.var = "log10.p.adj",
                        do.neg.log10 = F,
                        dep.n=6, 
                        enr.n=6, 
                        facet = F,
                        title.string="") {
  
  if (!is.null(dat.subset)) {
    test.df <- test.df %>% filter(grouping.var == dat.subset)
  }
  
  print(head(test.df[[p.value.var]]))
  if (do.neg.log10 == T) {
    test.df[[p.value.var]] <- -log10(test.df[[p.value.var]])
  }
  
  dff1 <- test.df %>%
    filter(!! sym(enrichment.var) <= 0) %>%
    slice_max(!! sym(p.value.var), n = dep.n)
  # slice_min((!! sym(enrichment.var)), n = dep.n)
  dff2 <- test.df %>% 
    filter(!! sym(enrichment.var) > 0 ) %>%
    slice_max(!! sym(p.value.var), n = enr.n) 
  # slice_max((!! sym(enrichment.var)), n = enr.n)
  dff <- rbind(dff1, dff2)
  
  max.match <- max(test.df$matches)
  min.match <- min(test.df$matches)
  
  gg <- ggplot(test.df, aes_string(x = enrichment.var, y = p.value.var)) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    geom_point(size = 2, shape = 21, color = "black", stroke = 0.2) +
    #, aes(fill = matches)) +
    geom_label_repel(data = dff, aes(label = short.name), size = 3) +
    theme_classic() +
    scale_fill_viridis(option="inferno", limits = c(min.match, max.match)) +
    labs(title = title.string)
  
  if (facet == T) { gg <- gg + facet_wrap(~cluster.num, scales = "free_y")}
  
  return(gg)
}


###############################################################################
#
# Bioinformatics functions
#
# #
# # # # # #
# # # # # # # # # #
# # # # # # # # # # # # # # # # # # # #
###############################################################################

plotEnhancerTiming <- function(gene, 
                               expression.mat = clean.glia.pb, 
                               atac.mat = glial.pb.atac.mat,
                               peak.set = peaks.gr,
                               links.object = sig.links, 
                               atac.ranges = rowRanges(atac.sce),
                               col.ord = order(atac.pb.colData$pseudotime)) {
  
  print("Setting up gene expression")
  gene <- toupper(gene)
  geneid <- convertGeneIDs(gene, name2ensembl)
  genexp <- data.frame(Gene.Expression = expression.mat[geneid, ])
  
  print("Setting up peaks")
  peaks <- findLinkedPeaks(gene, links.object, make.unique = T)
  peaks.gr <- atac.ranges[atac.ranges$name %in% peaks]
  peak.class <- distalProximalGenic(trim(peak.set), gene.gr)
  
  print("Smoothing ATAC matrix")
  raw.mat <- as.matrix(atac.mat[peaks, col.ord])
  
  smooth.mat <- do.call(rbind, lapply(
    apply(raw.mat, 1, lowess),
    function(x) { x$y } 
  ))
  colnames(smooth.mat) <- colnames(raw.mat)
  
  print("Getting peak order")
  maxes <- apply(smooth.mat, 1, function(x) which(x == max(x))[1]) %>%
    Reduce(c, .)
  
  ord.mat <- smooth.mat[order(maxes) , ]
  
  print(med(ord.mat))
  # return(ord.mat)
  print("Cleaning matrix")
  mat <- cleanMatrix(
    mat = ord.mat, 
    trim = T, 
    cuts = c(2,98), 
    scaling = 'row'
  )
  
  print("Classifying peaks")
  disp.mat <- mat[rowSums(mat) != 0, ]
  peak.class <- data.frame(class = peak.class[rowSums(mat) != 0])
  
  color.vector <- plant_palettes$`Strelitzia reginae`[1:3]
  names(color.vector) <- c('distal','proximal','genic')
  
  print("Creating heatmap")
  LA <- HeatmapAnnotation(
    df = peak.class,
    which = 'row', 
    col = list(class = color.vector)
  )
  HA <- HeatmapAnnotation(
    df = genexp, 
    col = list(Gene.Expression = colorRamp2(
      breaks = seq(min(disp.mat), max(disp.mat), length.out = 6), 
      colors = plant_palettes$Iris
    )),annotation_height = unit(1, 'cm')
  )
  
  H <- Heatmap(
    matrix = disp.mat, 
    cluster_rows = F, 
    cluster_columns = F, 
    top_annotation = HA,
    left_annotation = LA
  )
  print(H)
  
}

#' handling of orf gene names
validGenes <- function(genes) {
  
  out <- ifelse(grepl('orf', genes), genes, toupper(genes))
  return(out)
}

#

#' not implemented ... need something like this though
isEnsemble <- function(genes) {
  
  if(all(grepl('ENS[[:upper:]]+[[:digit:]]{11}', rownames(rna.mat)))) {
    return(TRUE)
  } else
    return(FALSE)
  
}

#' Quickly ID cells that have enough expression of some genes
cellsAboveThreshold <- function(genes, 
                                threshold = 0.05, 
                                thresh.as.quantile = T, 
                                object = logcounts(sce), 
                                convert = T) {
  
  genes <- validGenes(genes)
  
  # Gene ID input handling
  if(convert == T) {
    gene.ids <- convertGeneIDs(genes, id.mapping.table = name2ensembl)
  } else {
    gene.ids <- genes
  }
  
  list.cells <- lapply(gene.ids, function(x) {
    
    # gene threshold
    if (thresh.as.quantile == T) {
      th.x <- quantile(object[x, ], probs = threshold)
    } else {
      th.x <- threshold
    }
    
    return(which(object[x, ] > th.x))
    
  })
  
  return(Reduce(intersect, list.cells))
  
}


#' Prioritize motifs' regulatory effect at a given gene
#' 
#' Not really effective yet, unless there are lots of peaks linked 
keyMotifs <- function(genes, 
                      num.obj = sig.links, 
                      denom.obj = all.links, 
                      denomenator = NULL) {
  
  # Get linked peaks 
  
  linked.peaks <- findLinkedPeaks(genes, object = num.obj)
  
  # Compute the number of motif matches in linked peaks
  
  numerator <- colSums(motif.mat[linked.peaks, ])
  
  # If there's no explicit background, compute it with variable glial genes
  
  if (is.null(denomenator)) {
    denomenator <- colMeans(motif.mat[findLinkedPeaks(var.glial.genes, denom.obj), ])
  }
  
  # The ratio of matches at the gene of interest, to the average number of 
  # matches per gene in the background 
  
  ratio <- numerator / denomenator
  
  out <- data.frame(
    motif.id = colnames(motif.mat),
    motif.name = add_split_name_vector(colnames(motif.mat), "_", 2),
    n.motifs = numerator,
    bg.motifs = denomenator,
    ratio = ratio,
    Z.score = scale_this(ratio),
    stringsAsFactors = F
  ) %>% arrange(-ratio)
  
  return(out)
  
}


#' Need to wrangle Gencode in order to get RNA-seq TPM values 
#' (comparison across samples is critical)
#' 
#' May be useful again some time 
getGencodeGeneLengths <- function(path = 'gbm/tcga/rna/counts/gencode.v22.annotation.gtf') {
  
  gencode <- rtracklayer::import.gff(path)
  exons <- gencode[gencode$type == 'exon' & seqnames(gencode) != 'chrM']
  
  # A non-overlapping exon set
  
  reduced.exons <- reduce(split(exons, ~gene_id))
  
  gene.lengths <- lapply(seq_along(reduced.exons), function(x) {
    sum(width(ranges(reduced.exons[[x]])))
  }) %>% Reduce(c, .)
  
  names(gene.lengths) <- names(reduced.exons)
  
  return(gene.lengths)
}

#' from FM
getUwotModelFromSeurat <- function(object, 
                                   assayUsed = DefaultAssay(object), 
                                   reductionUsed = "pca", ...){
  
  require(uwot)
  
  paramL <- object@commands[[paste0("RunUMAP.",assayUsed,".",reductionUsed)]]@params
  
  if (is.null(paramL[["dims"]])) logger.error("Invalid dims")
  reqParamNames <- c(
    "n.neighbors", "n.components", "metric", "n.epochs", "learning.rate", "min.dist",
    "spread", "set.op.mix.ratio", "local.connectivity", "repulsion.strength",
    "negative.sample.rate", "a", "b", "uwot.sgd", "seed.use"
  )
  
  # assign NULL to missing parameters
  for (pn in reqParamNames){
    if (!is.element(pn, names(paramL))) paramL[pn] <- list(NULL)
  }
  
  X <- Embeddings(object[[reductionUsed]])[, paramL[["dims"]]]
  
  # make sure to use the same 'random' numbers
  
  if (!is.null(paramL$seed.use)) {
    set.seed(seed = paramL$seed.use)
  }
  
  umapRes <- umap(
    X = X,
    n_neighbors = as.integer(paramL$n.neighbors),
    n_components = as.integer(paramL$n.components),
    metric = paramL$metric,
    n_epochs = paramL$n.epochs,
    learning_rate = paramL$learning.rate,
    min_dist = paramL$min.dist,
    spread = paramL$spread,
    set_op_mix_ratio = paramL$set.op.mix.ratio,
    local_connectivity = paramL$local.connectivity,
    repulsion_strength = paramL$repulsion.strength,
    negative_sample_rate = paramL$negative.sample.rate,
    a = paramL$a,
    b = paramL$b,
    fast_sgd = paramL$uwot.sgd,
    ret_model=TRUE, # this is the important part
    ...
  )
  
  return(umapRes)
}

#' from FM
projectMatrix_SeuratUMAP <- function(X, 
                                     object, 
                                     umap_model = NULL,
                                     assayUsed = DefaultAssay(object), 
                                     dataType = "counts",
                                     return.object = F) {
  
  # X <- GetAssayData(object=object[["RNA"]], slot="counts") # for testing
  # colnames(X) <- paste0("cell", 1:ncol(X))
  
  if (!is.element(dataType, c("counts", "normalized", "scaled", "PCA"))){
    stop(sprintf("Invalid dataType: %s", dataType))
  }
  
  paramL_norm <- object@commands[[paste0("NormalizeData.",assayUsed)]]@params
  paramL_scale <- object@commands[[paste0("ScaleData.",assayUsed)]]@params
  
  if (!paramL_scale[["do.scale"]]) stop("Don't know how to deal with unscaled data yet")
  
  paramL_pca <- object@commands[[paste0("RunPCA.",assayUsed)]]@params
  paramL_umap <- object@commands[[paste0("RunUMAP.",assayUsed,".",paramL_pca[["reduction.name"]])]]@params
  seuQ <- X # query Seurat dataset
  
  if (class(X) != "Seurat"){
    seuQ <- CreateSeuratObject(counts = X, project = "seuQ", assay = assayUsed)
  }
  
  print(dim(seuQ))
  if (is.element(dataType, c("counts"))){
    cat("Normalizing ...\n")
    seuQ <- NormalizeData(
      seuQ,
      assay = paramL_norm$assay,
      normalization.method = paramL_norm$normalization.method,
      scale.factor = paramL_norm$scale.factor,
      verbose = FALSE
    )
  }
  if (!is.element(dataType, c("scaled", "PCA"))){
    cat("Scaling ...\n")
    seuQ <- ScaleData(
      seuQ,
      assay = paramL_scale$assay,
      features = paramL_scale$features,
      do.scale = paramL_scale$do.scale,
      do.center = paramL_scale$do.center,
      model.use = paramL_scale$model.use,
      verbose = FALSE
    )
  }
  if (!is.element(dataType, c("PCA"))){
    
    cat("PC projection ...\n")
    X <- GetAssayData(seuQ[[assayUsed]], slot="scale.data")
    features <- rownames(object@reductions$pca@feature.loadings)
    
    if (any(!features %in% rownames(X))){
      stop("Could not find all features in X for the PC projection")
    }
    
    X <- t(X[features,])
    
  }
  
  projM <- Loadings(object, reduction = paramL_pca[["reduction.name"]])
  pcaCoord_proj <- X %*% projM
  
  # pcaCoord_orig <- Embeddings(object[[paramL_pca[["reduction.name"]]]]) # compare the original
  
  cat("Retrieving UMAP model ...\n")
  if (is.null(umap_model)) {
    umapRes <- getUwotModelFromSeurat(object, assayUsed=assayUsed, reductionUsed=paramL_pca[["reduction.name"]])
  } else {
    umapRes <- umap_model
  }

  cat("UMAP projection ...\n")
  umapCoord_orig <- Embeddings(object[[paramL_umap[["reduction.name"]]]])
  
  return(list(orig.umap = umapCoord_orig, proj.pca = pcaCoord_proj))
  
  umapCoord_proj <- uwot::umap_transform(pcaCoord_proj[ , paramL_umap[["dims"]]], umapRes)
  
  rownames(umapCoord_proj) <- rownames(pcaCoord_proj)
  colnames(umapCoord_proj) <- colnames(umapCoord_orig)
  
  res <- list(
    pcaCoord=pcaCoord_proj,
    pcsUsed=paramL_umap[["dims"]],
    umapCoord=umapCoord_proj
  )
  
  if (return.object == T) {
    s.out <- list(objectrat = seuQ)
    res <- append(res, s.out)
  }
  
  return(res)
}

#' Grab genes from a particular gene ontology
#' Wrapped around 'clusterProfiler'
#' 
#' @param x GO path ID
#' @param OrgDb e.g. 'org.Hs.eg.db' database
#' @param ont ontology to use
#' @param keytype symbol? Entrez?
getGOgeneSet <- function(x, 
                         OrgDb = "org.Hs.eg.db", 
                         ont = "BP", 
                         keytype = "SYMBOL") {
  
  goObj <- clusterProfiler:::get_GO_data(OrgDb=OrgDb, ont=ont, keytype=keytype)
  mapObj <- goObj$PATHID2EXTID
  
  if (!is.element(x, names(mapObj))){
    tt <- names(goObj$PATHID2NAME)
    names(tt) <- goObj$PATHID2NAME
    x <- tt[x]
  }
  
  return(sort(unique(mapObj[[x]])))
}

#' Louvain clustering using Seurat
#' 
#' From Fabian Mueller. Takes a matrix and a clustering resolution
#' @param X a matrix
#' @param resolution A clustering resolution
clusterSeurat <- function(X, resolution = 0.8) {
  
  dummyMat <- matrix(11, ncol=ncol(X), nrow = 11)
  colnames(dummyMat) <- colnames(X)
  rownames(dummyMat) <- paste0("df", 1:nrow(dummyMat))
  
  sObj <- Seurat::CreateSeuratObject(
    counts = dummyMat, 
    project = 'clustering', 
    min.cells = 0, 
    min.features = 0, 
    assay = "clust"
  )
  sObj[["clustpca"]] <- Seurat::CreateDimReducObject(
    embeddings = t(X), 
    key="dim_", 
    assay="clust"
  )
  sObj <- Seurat::FindNeighbors(
    object = sObj, 
    reduction = "clustpca", 
    assay = "clust", 
    verbose = FALSE
  )
  sObj <- Seurat::FindClusters(
    object = sObj, 
    resolution = resolution, 
    verbose = FALSE
  )
  cluster.assignment <- getClusterFactor(sObj@meta.data[,"seurat_clusters"])
  return(cluster.assignment)
}

#' Naming clusters and factorizing 
#' 
#' From FM, supports clusterSeurat
#' @param x a vector
#' @param prefix a prefix to add to factors 
getClusterFactor <- function(x, prefix="c") {
  if (is.factor(x)){
    if (!all(grepl(paste0("^", prefix), levels(x)))) levels(x) <- paste0(prefix, levels(x))
  } else if (is.numeric(x) || is.integer(x)){
    x <- factor(paste0(prefix, x), levels=paste0(prefix, sort(unique(x))))
  } else {
    x <- factor(x)
    if (!all(grepl(paste0("^", prefix), levels(x)))) levels(x) <- paste0(prefix, levels(x))
  }
  return(x)
}

#' Make NN-based pseudobulk samples
#' 
#' Takes the data and makes pseudo-bulk samples to reduce the complexity of 
#' many cells for plotting / interpretation. Adapted from FM.
#' @param grouping A vector of groupings, which will inform the proportions
#' of cells randomly selected to seed pseudobulks
#' @param mat A matrix, e.g. an expression matrix, across which cell data will
#' be aggregated using 'aggr'
#' @param percSub The downsampling percentage. A value of 0.1 with 10000 cells
#' will produce 1000 pseudo-bulk samples
#' @param nSample The number of cells to sample per PB
#' @param aggr a function that aggregates cells.
samplePseudoBulkNN <- function(grouping, 
                               nn.object,
                               mat = NULL, 
                               percSub = 0.1, 
                               nSample = 50, 
                               aggr = "sum"){
  
  # Start by getting the cluster groupings
  cat('Get groupings \n')
  nPerGroup <- table(grouping)
  nPerGroup <- nPerGroup[nPerGroup > 0] # exclude unoccupied factor levels if present
  
  # if not specified, take half of the minimum group size as sample size
  if (is.null(nSample)) {
    nSample <- ceiling(min(nPerGroup) / 2)
    logger.info(c("Using a sample size of", nSample))
  }
  
  nSamplePerGroup <- ceiling(nPerGroup * percSub)
  
  # Apply over all groups
  
  cat('Apply over all groups \n')
  
  sampleIdxL <- lapply(
    X   = names(nSamplePerGroup), 
    FUN = function(gn) {
      
      # Indices of all the cells in the group
      idx <- which(grouping == gn)

      # N: num cells in group, or the sample size
      N <- min(c(length(idx), nSample))
      
      rr <- lapply(
        1:nSamplePerGroup[gn], 
        FUN = function(i) {
          
          random.cell <- sample(idx, 1)
          
          nn <- sampleNearestNeighbors(
            cell = random.cell, 
            knn.object = nn.object, 
            N = N,
            output = 'indices'
          )
          return(nn)
        }
      )
      
      return(rr)
      
    }
  )
  
  sampleIdxL <- unlist(sampleIdxL, recursive = FALSE)
  
  aggrFun <- NULL
  if (aggr=="sum"){
    aggrFun <- function(X){rowSums(X, na.rm=TRUE)}
  } else if (aggr=="mean"){
    aggrFun <- function(X){rowMeans(X, na.rm=TRUE)}
  } else {
    logger.error(c("Unknown aggregation function:", aggr))
  }
  resM <- NULL
  if (!is.null(mat)){
    # id a matrix is supplied, aggregate across the matrix
    resM <- do.call(
      "cbind", 
      lapply(
        X = sampleIdxL, 
        FUN=function(x) {
          aggrFun(mat[,x,drop=FALSE])
        }
      )
    )
  }
  
  res <- list(
    matrix = resM,
    sampleIdxL = sampleIdxL,
    grouping = grouping
  )
  #class(res) <- "SampledPseudobulk"
  return(res)
}

#' A function to sample nearest neighbors of a particular cell
#' 
#' Takes a cell ID or index, and a KNN object (output from URD in this case), 
#' and returns a specified number of cells that are nearest neighbors of that 
#' cell. Cell selection can be biased according to a function of the cell-cell
#' distances (FUN)
#' @param cell A cell ID or idx
#' @param knn.object An S4 object with names nn.idx, nn.cells, and nn.dists,
#' providing NN info
#' @param N The number of cells to sample
#' @param FUN A quoted string, which evaluates to an expression for weighting
#' the probability of randomly selecting cells. Allows you to weight by 
#' distance for example.
#' @param output String, either 'indices' or not, which specifies whether 
#' numeric cell indices or IDs will be used.
sampleNearestNeighbors <- function(cell, 
                                   knn.object, 
                                   N = 50, 
                                   FUN = "(1/exp(dists)) * (1/sum(1/exp(dists)))",
                                   output = 'indices') {
  
  
  if (output == 'indices') {
    cells <- knn.object[['nn.idx']][cell, ]
  }
  else {
    cells <- knn.object[['nn.cells']][cell, ]
  }
  
  dists <- knn.object[['nn.dists']][cell, ]
  
  sample(
    x = cells, 
    size = N, 
    replace = F, 
    prob = eval(parse(text = FUN))
  )
  
}

#' Retrieve the dominant grouping present in a group of single cells
#' 
#' This is useful for e.g. getting the most prominent cluster in a pseudobulk
#' sample
#' @param cells a list of cell names or indices into 'groupings'
#' @param groupings a vector of labels, for example cluster labels
dominantCluster <- function(cells, groupings) {
  cl.tbl <- table(groupings[cells])
  out <- names(cl.tbl)[cl.tbl == max(cl.tbl)]
  if(length(out) == 1) 
    return(out)
  else {
    return(out[1])
  }
}


#' For a given TF motif, return all peaks with that motif
#' 
#' Takes a character string and a binary motif matrix (rows are peaks, columns
#' are the motif names) and returns a vector of characters - the peak names - 
#' that have that motif
#' @param tf
#' @param motif.matrix
peaksContainingMotif <- function(tf, motif.matrix) {
  
  matches <- grep(tf, colnames(motif.matrix), value = T)
  
  if(length(matches) != 1) {
    
    sprintf("%s string matches for that motif", as.character(length(matches)))
    
    return(NA)
  } else {
    
    tf <- matches
    out <- rownames(motif.matrix)[motif.matrix[,tf]]
    
    return(out)
  }
}



#' For a gene list, get the expression-weighted pseudotime mean 
#'  
#' Purpose is to see "when" the genes are expressed
#' 
#' @param genes A vector of genes
#' @param object A data matrix
#' @param scale.expression Whether or not to standardize expression prior to 
#' computing weighted means.
#' @param pseudotimes A vector of pseudotimes for cells in SCE object - subject
#' to and matching indices with "cell.subset" if applicable
#' @param cell.subset
#' @param id.type Handles gene input into the expression matrix, which indexed 
#' as Ensembl IDs by default
#' @param p Exponent of weights vector. You can exaggerate the expression here
genePseudotimes <- function(genes, 
                            object, 
                            scale.expression = F,
                            pseudotimes,
                            cell.subset,
                            id.type = 'symbol',
                            p = 1) {
  
  if (id.type == 'symbol') {
    genes <- convertGeneIDs(genes, name2ensembl)
  }
  
  filt.ind <- filterGenesNotInData(
    genes = genes,
    object = object,
    to.return = 'indices'
  )
  gene.ids <- genes[filt.ind]
  
  rna.mat <- object[gene.ids , cell.subset]
  
  if (scale.expression == T) {
    rna.mat <- rowScale(rna.mat)
  }
  
  weighted.pt <- do.call(
    c,
    lapply(
      seq_along(gene.ids),
      function(x) {
        if (x %% 100 == 0) { cat(sprintf("%s...", as.character(x))) }
          weighted.mean(x = pseudotimes, w = (rna.mat[gene.ids[x], ])^p)
      }
    )
  )
  cat("\n")
  return(weighted.pt)
  
}

#' Calculate a density score for Linked Peaks near a gene
#' 
#' Density score is simply n.peaks per megabase. 
#' @param gene A gene
#' @param object An object to pass to findLinkedPeaks.
#' @param peak.ranges A GRanges object that gives peak locations (by peak name)
getPeakDensities <- function(gene, object, peak.ranges) {
  # todo
}

#' Make NN- or factor-based pseudobulk samples
#' 
#' Takes the data and makes pseudo-bulk samples to reduce the complexity
#' of many many cells for plotting / interp
#' From FM
samplePseudoBulk <- function(grouping, 
                             mat = NULL, 
                             percSub = 0.1, 
                             nSample = 50, 
                             aggr = "sum"){
  
  # if (class(grouping)=="SampledPseudobulk"){
  # sampling already done. Use preexisting samples
  # sampleIdxL <- grouping$sampleIdxL
  # grouping <- grouping$grouping
  # } else {
  # Sample from the grouping if not done already
  nPerGroup <- table(grouping)
  nPerGroup <- nPerGroup[nPerGroup > 0] # exclude unoccupied factor levels if present
  # if not specified, take half of the minimum group size as sample size
  if (is.null(nSample)) {
    nSample <- ceiling(min(nPerGroup) / 2)
    logger.info(c("Using a sample size of", nSample))
  }
  nSamplePerGroup <- ceiling(nPerGroup * percSub)
  
  sampleIdxL <- lapply(
    X   = names(nSamplePerGroup), 
    FUN = function(gn) {
      
      idx <- which(grouping==gn)
      N <- min(c(length(idx), nSample))
      
      rr <- lapply(
        1:nSamplePerGroup[gn], 
        FUN = function(i) {
          sort(sample(idx, N))
        }
      )
      return(rr)
    }
  )
  
  sampleIdxL <- unlist(sampleIdxL, recursive = FALSE)
  
  aggrFun <- NULL
  if (aggr=="sum"){
    aggrFun <- function(X){rowSums(X, na.rm=TRUE)}
  } else if (aggr=="mean"){
    aggrFun <- function(X){rowMeans(X, na.rm=TRUE)}
  } else {
    logger.error(c("Unknown aggregation function:", aggr))
  }
  resM <- NULL
  if (!is.null(mat)){
    # id a matrix is supplied, aggregate across the matrix
    resM <- do.call("cbind", lapply(sampleIdxL, FUN=function(x){
      aggrFun(mat[,x,drop=FALSE])
    }))
  }
  res <- list(
    matrix=resM,
    sampleIdxL=sampleIdxL,
    grouping=grouping
  )
  #class(res) <- "SampledPseudobulk"
  return(res)
}

#' Compute differential counts or norm counts between two groups of cells
#' 
#' 
differentialWilcox <- function(group1, 
                               group2, 
                               dat = sce@assays$data$logcounts,
                               min.counts = 0,
                               atac = F,
                               diffs = 'fold change') {
  
  cat("\nFiltering zero-count features across both groups...\n\n")
  group.union <- base::union(group1, group2)
  condition <- rowSums(dat[ , group.union]) > min.counts
  dat <- dat[condition , ]
  print(dim(dat))
  
  cat(
    sprintf(
      "Computing rank-sums on %s features...\n", 
      as.character(sum(condition))
    )
  )
  
  # this is weird, todo: fix this hardcoded thing
  if (atac == T) {
    rns <- rowRanges(atac.sce)$name
    rns <- rns[condition]
  } else{
    rns <- rownames(dat)
  }
  
  out <- do.call(
    rbind, 
    lapply(
      1:nrow(dat), 
      function(x) {
        if(x %% 10 == 0) {
          cat(paste0(c("\t", as.character(x), " features\n")))
        }

        gp1 <- dat[x, group1]
        gp2 <- dat[x, group2]
        avg1 <- mean(gp1)
        avg2 <- mean(gp2)
        
        wt <- wilcox.test(gp1, gp2)
        
        if (diffs == 'fold change') {
          diff <- log2(avg1 / avg2)
          
          out <- data.frame(
            gene.id = rns[x], 
            W.statistic = wt$statistic, 
            group.1.avg = avg1,
            group.2.avg = avg2,
            avg.log2.FC = diff,
            p.value = wt$p.value, 
            row.names = F,
            stringsAsFactors = F
          )
          
        } else {
          diff <- avg1 - avg2
          
          out <- data.frame(
            gene.id = rns[x], 
            W.statistic = wt$statistic, 
            group.1.avg = avg1,
            group.2.avg = avg2,
            differential = diff,
            p.value = wt$p.value, 
            row.names = F,
            stringsAsFactors = F
          )
          
        }
      }
    )
  )
  out[['p.adj']] <- p.adjust(out[['p.value']], method = "BH")
  return(out)
}

#' Retrieve accessibility peaks linked to target genes
#' 
#' For use with a data frame containing peak to gene linkages, already 
#' filtered. Or, for a background set, all possible linkages.
#' @param genes a vector of gene symbols. 
#' @param object A data frame where each row is a peak to gene link
#' @param make.unique Whether to return only 1 instance of each peak, 
#' regardless of how many linkages there are.
findLinkedPeaks <- function(genes, 
                            object, 
                            make.unique = T) {
  pks <- object %>% filter(gene.symbol %in% genes) %>% pull(peak.name)
  if(make.unique == T) {
    pks <- unique(pks)
  }
  return(pks)
}

#' Retrieve accessibility peaks linked to target genes
#' 
#' For use with a data frame containing peak to gene linkages, already 
#' filtered. Or, for a background set, all possible linkages.
#' @param peaks a vector of peak ids. 
#' @param object A data frame where each row is a peak to gene link
#' @param make.unique Whether to return only 1 instance of each gene, 
#' regardless of how many linkages there are.
findLinkedGenes <- function(peaks, 
                            object, 
                            make.unique = T) {
  genes <- object %>% filter(peak.name %in% peaks) %>% pull(geneSymbol)
  if(make.unique ==T) {
    genes <- unique(genes)
  }
  return(genes)
}

#' Computes differences between two peak sets.
#' @param peak.set1
#' @param peak.set2

#' Find all cells expressing any of the target gene. Can set the threshold
#' 
#' This function converts a given gene symbol to ENSEMBL. easy to adjust that 
#' if needed. 
#' Also currently uses log counts for threshold units
#' 
#' @param x a gene name. 
findTargetCellsHelper <- function(x, 
                                  threshold = 0){
  
  id <- convertGeneIDs(x, name2ensembl)
  cells <- rownames(colData(sce))[sce@assays$data$logcounts[id,] > threshold]
  return(cells)
  
}

#' Add an arbitrary feature/vector to a data frame. Expression or vector
#' must be (or evaluate to) the number of rows in the projection. Projection
#' will be coerced to data.frame. 
#' @param projection A projection of cell coordinates like UMAP coordinates
#' @param expr an expression, which will be parsed and evaluated to produce
#' a vector of values; or a vector, with the same number of elements as there 
#' are rows in the projection.
#' @param col.name the column name of the new variable.
addFeatureToProjection <- function(projection, 
                                   expr, 
                                   col.name = "variable") {
  
  if(class(expr) == "character") {
    added <- eval(parse(text = expr))
    print(head(added))
  } 
  if(class(expr) == "vector") {
    added <- expr
  }
  
  return(data.frame(projection, col.name = added))
  
}

#' Deprecated
#' Given a list of genes, remove names that aren't valid or not present
#' 
#' Alternatively, return indices that correspond to the valid genes. 
#' @param genes a vector of gene names
#' @param object a data matrix.
#' @param key.type The input type
#' @param convert.back Whether to return genes converted back to symbols
#' @param to.return What to return
filterGenesNotInDataOld <- function(genes, 
                                 object = logcounts(sce), 
                                 key.type = "ensembl", 
                                 convert.back = F,
                                 to.return = 'sym') {
  

  if(grepl("sym|name", key.type)) {
    genes <- convertGeneIDs(genes, name2ensembl)
  }
  indices <- geneInData(genes, object = object)
  genes <- genes[indices]
  if(convert.back == T) {
    genes <- convertGeneIDs(genes, ensembl2name)
  }
  if(grepl("ind|indices", to.return)) {
    return(indices)
  } else {
    return(genes)
  }
}

#' Given a list of genes, remove names that aren't valid or not present
#' 
#' Alternatively, return indices that correspond to the valid genes. 
#' @param genes a vector of gene names
#' @param object a data matrix.
#' @param to.return What to return
filterGenesNotInData <- function(genes, 
                                 object = logcounts(sce), 
                                 to.return = 'sym') {
  
  indices <- geneInData(genes, object = object)
  genes <- genes[indices]
  
  if(grepl("ind|indices", to.return)) {
    return(indices)
  } else {
    return(genes)
  }
}

#' Given a gene and a SCE dataset, test if the gene is present
#' @param gene A gene ID
#' @param object A matrix or SCE object
geneInData <- function(gene, 
                       object = logcounts(sce)) {
  return(gene %in% rownames(object))
}

#' Retrieve average expression of a gene set in single cells
#' 
#' @param genes a vector of gene names
#' @param object a SingleCellExperiment object
#' @param convert whether or not to convert names to ENSEMBL IDs
averageExpression <- function(genes, 
                              object = sce, 
                              convert = T) {
  
  if (convert == T) { 
    gene.ids = convertGeneIDs(genes, name2ensembl) 
  } else {
    gene.ids = genes
  }
  
  if (length(gene.ids) == 0) {
    gene.exp <- NA
  } 
  if (length(gene.ids) == 1) {
    gene.exp <- logcounts(sce)[gene.ids, ]
  }
  if (length(gene.ids) > 1) {
    gene.exp <- colMeans(logcounts(sce)[gene.ids, ])
  }
  
  return(gene.exp)
  
}

#' Function to retrieve significantly correlated peaks (by name) related to a 
#' gene list
correlatedPeaks <- function(genes, 
                            corr.object) {
  
}

#' Helper function to perform row-wise correlations between two matched matrices
#' 
#' Computes correlations. Typically this will be comparing two cross-mapped datasets
#' e.g. accessibility-derived gene activity and RNA expression
#' @param mat1 A numeric matrix
#' @param mat2 A (dim-matched) numeric matrix
#' @param method Pearson, spearman, or other correlation (input to "cor")
rowwiseCorrelationsHelper <- function(mat1, 
                                      mat2, 
                                      method) {
  out <- do.call(c, lapply(1:nrow(mat1), function(x) {
    cor(mat1[ x , ], mat2[ x , ], method = method)
  }))
  return(out)
}

#' Performd row-wise correlations between two matched matrices
#' 
#' Typically this will be comparing two cross-mapped datasets
#' e.g. accessibility-derived gene activity and RNA expression
#' @param mat1 A numeric matrix
#' @param mat2 A (dim-matched) numeric matrix
#' @param method Pearson, spearman, or other correlation (input to "cor")
#' @param max.load The computation slows considerably at scale. This chunks the
#' original matrices into bins of a maximum size, which are computed separately
#' and concatenated.
rowwiseCorrelations <- function(mat1, 
                                mat2, 
                                method = "pearson", 
                                max.load = NULL) {
  
  if(dim(mat1)!=dim(mat2)) {
    cat("Matrices do not match\n")
    return(NULL)
  } 
  
  if(!is.null(max.load)) {
    
    cat("splitting sub-matrices...\n")
    n.bins <- round(nrow(mat1)/max.load)
    b <- binr::bins(1:nrow(mat1), target.bins = n.bins, max.breaks = n.bins+5)
    lows <- b[['binlo']]
    highs <- b[['binhi']]
    
    out <- do.call(c, lapply(1:length(lows), function(x){
      
      cat(sprintf("bin %s:", x))
      bin.range <- lows[x]:highs[x]
      cat(sprintf("\t%s-%s\n", as.character(lows[x]), as.character(highs[x]) ))
      rowwiseCorrelationsHelper(mat1[bin.range, ], mat2[bin.range, ], method = method)
      
    }))
  } else {
    out <- rowwiseCorrelationsHelper(mat1, mat2, method = method)
  }
  
  names(out) <- rownames(mat2)
  
  return(out) 
}

#' Hypergeometric test helper function
#' 
#' Performs the actual computation. 
#' @param ii An identifier for the test group. It's used to return a list of 
#' indices from "test.groups". Indices will subset a draw/match matrix 
#' ("test.match_matrix") into foreground (matching ii) and background 
#' (!matching ii)
#' @param test.groups Atomic vector of groups corresponding to the rows in 
#' "test.match_matrix"
#' @param test.match_matrix A binary matrix or sparse matrix of matches/
#' nonmatches
#' @param test.lower_tail T or F, whether to test for the lower tail 
#' (depletion) or not 
computeHypergeometricTest <- function(ii, 
                                      test.groups, 
                                      test.match_matrix, 
                                      test.lower_tail) {
  
  # variable names: enrichment or depletion?
  suffix <- ifelse(test.lower_tail==TRUE, "depletion", "enrichment")
  pvalue.text <- paste0(suffix, ".p")
  group.text <- paste0(suffix, ".group")
  
  # define foreground & bg
  draw <- test.groups == ii
  bg <- test.groups != ii
  
  cat(paste0("Group ", ii, ",\t", "size ", sum(draw), "\n"))
  
  # phyper params
  k <- colSums(test.match_matrix[draw, ]) # Number of matches in test group
  m <- colSums(test.match_matrix[bg, ]) # Number of matches in background
  n <- length(bg) - m # Number of non-matches in background
  N <- sum(draw) # Number of draws
  
  # output
  df <- data.frame(tidy(phyper(k-1, m, n, N, 
                               log.p = F, 
                               lower.tail = test.lower_tail)), 
                   matches = k, 
                   draws = N, 
                   bg.matches = m, 
                   bg.non.matches = n, 
                   group = ii) %>%
    dplyr::rename((!!pvalue.text) := "x")
  
  return(df)
  
}

#' Hypergeometric test for enrichment in a list of groups. 
#' 
#' @param groupings A vector of possible groupings, e.g. cluster designations.
#' Defines the group assignment for each row of a match matrix. For example: 
#' c("foreground", "foreground", "background", "background", "foreground")
#' @param match_matrix A binary matrix or sparse matrix of matches/
#' nonmatches
#' @param lower_tail T, F, or "both". Whether to test for lower tail enrichment
#' (depletion) or not.
hypergeometricTestList <- function(groupings, 
                                   match_matrix, 
                                   lower_tail = F) {
  
  tbl <- table(groupings)
  the.groups <- do.call(c, lapply(names(tbl), function(x) {
    if(tbl[x] < 10) {
      msg <- sprintf("Only %s draws in group %s, skipping group\n", as.character(tbl[x]), x)
      warning(msg)
      return(NULL)
    } else {
      return(x)
    }
  }))
  
  the.groups <- the.groups[!is.null(the.groups)]
  
  if (is.logical(lower_tail)) {
    cat("single\n")
    output <- do.call(rbind, lapply(X = the.groups, 
                                    FUN = computeHypergeometricTest, 
                                    test.groups = groupings, 
                                    test.match_matrix = match_matrix, 
                                    test.lower_tail = lower_tail))
  } else {
    cat("\nEnrichment: \n")
    out <- do.call(rbind, lapply(X = the.groups, 
                                 FUN = computeHypergeometricTest, 
                                 test.groups = groupings, 
                                 test.match_matrix = match_matrix, 
                                 test.lower_tail = F))
    cat("\nDepletion: \n")
    out.lt <- do.call(rbind, lapply(X = the.groups, 
                                    FUN = computeHypergeometricTest, 
                                    test.groups = groupings, 
                                    test.match_matrix = match_matrix, 
                                    test.lower_tail = T)) 
    
    out.lt <- out.lt %>% dplyr::select(depletion.p, names, group)
    
    output <- inner_join(out, out.lt)
    
  }
  
  return(output)
}

# Chi squared test for matches in a foreground vs background motif matrix
bg_chisq_test <- function(fg_motif_matrix, 
                          bg_motif_matrix) {
  k <- colSums(fg_motif_matrix)
  m <- colSums(bg_motif_matrix)
  n <- nrow(bg_motif_matrix)
  N <- nrow(fg_motif_matrix)
  df <- do.call(rbind, lapply(1:length(k), function(x) {
    cbind(names=colnames(fg_motif_matrix)[x], 
          tidy(chisq.test(matrix(c(k[x], N, m[x], n), nrow = 2))),
          log2.fold.enrich = log2((k[x]/N)/(m[x]/n)),
          fg.matches = k[x],
          fg.size = N,
          bg.matches = m[x],
          bg.size = n)
  }))
}

#' Subset SCE cells by their metadata
#' 
#' Given a match (like "cluster1") and the right columns etc, this will
#' subset the data appropriately. Any metadata column.
#' @param metadata.match
#' @param mat
#' @param metadata
#' @param metadata.column
#' @param grep
subsetCellsByMetadata <- function(metadata.match, 
                                  mat, 
                                  metadata, 
                                  metadata.column) {
  
  if (!length(metadata.match) == length(metadata.column)) {
    stop("Option metadata.match should the same length as metadata.column")
  }
  
  index.list <- lapply(
    1:length(metadata.match),
    function(x) {
      colnames(mat)[ metadata[[ metadata.column[[x]] ]] %in% metadata.match[[x]] ]
    }
  )
  
  indices <- Reduce(intersect, index.list)
  sub.mat <- mat[ , indices]
  
  return(sub.mat)
}



motif_mat_from_df <- function(df, 
                              genome = "hg38", 
                              motifs = motifs, 
                              bkgd = "subject", 
                              pval.cutoff = 5e-5, 
                              row.names = T, 
                              row.name.col = "peak.name") {
  
  df.r <- makeGRangesFromDataFrame(
    df, 
    starts.in.df.are.0based = T, 
    ignore.strand = T, 
    keep.extra.columns = T
  )
  
  df.motifs <- matchMotifs(
    motifs, 
    df.r, 
    genome = genome, 
    out = "scores", 
    p.cutoff = pval.cutoff, 
    bg = bkgd
  )
  
  df.mat <- as.matrix(motifMatches(df.motifs))
  
  if (row.names == T) { 
    rownames(df.mat) <- df[[row.name.col]] 
  }
  
  return(df.mat)
}

#' Modified from Alicia Schep's Motifmatchr function. Now enables new versions of JASPAR
#' to be called. let's see if this holds up over time...
#' 
#' @param species Two word species designation
#' @param collection CORE collection or Unvalidated
#' @version "JASPAR2016", "JASPAR2018", "JASPAR2020" etc
getNewJasparMotifs <- function(species = "Homo sapiens", collection = "CORE", version = "JASPAR2018") {
  
  opts <- list()
  opts['species'] <- species
  opts['collection'] <- collection
  
  motifs <- TFBSTools::getMatrixSet(eval(parse(text=paste0(version, "::", version))), opts = opts)
  
  if (!isTRUE(all.equal(TFBSTools::name(motifs), names(motifs)))) { 
    names(motifs) <- paste(names(motifs), TFBSTools::name(motifs), sep = "_")
  }
  
}

#' CCA nearest neighbor. 
findNearestNeighbor <- function(cell.id, nn.mapping = map.acc2expr) {
  return(mapValues(x = cell.id, mapping.table = nn.mapping))
}

#' Given a list of gene IDs, map the values to another ID type using a table
#' 
#' @param genelist A list of genes (vector)
#' @param id.mapping.table A named vector, where names are of the input ID type, and
#' values are of the desired output ID type
convertGeneIDs <- function(genelist, id.mapping.table) {
  return(mapValues(x = genelist, mapping.table = id.mapping.table))
}

#' Test whether peaks are distal or proximal to gene
#' 
#' Function takes two granges files, peaks and gene model, and applies a TSS cutoff
#' @param peak.granges
#' @param gene.granges
#' @param tss.cutoff 
distalProximalGenic <- function(peak.granges, gene.granges, tss.cutoff = c(1000,200)) {
  
  # tss.gr <- GenomicRanges::resize(gene.granges, width = 1, fix = "start")
  promoter.gr <- GenomicRanges::promoters(x = gene.granges, upstream = tss.cutoff[1], downstream = tss.cutoff[2])
  dist.to.promoter <- elementMetadata(GenomicRanges::distanceToNearest(peak.granges, promoter.gr))[['distance']]
  dist.to.gene <- elementMetadata(GenomicRanges::distanceToNearest(peak.granges, gene.granges))[['distance']]
  
  iter <- 1:length(peak.granges)
  peak.class <- do.call(
    c,
    lapply(
      iter,
      function(x) {
        if(x %% 50000 == 0) {
          cat(sprintf("%s peaks\n", as.character(x)))
        }
        
        if(dist.to.promoter[x] == 0) {
          out = "proximal"
        } else if (dist.to.gene[x] == 0) {
          out = "genic"
        } else {
          out = "distal"
        }
        return(out)
        
      }
    )
  )
  return(peak.class)
}

getTFFamilies <- function(name1, 
                          name2, 
                          lookup.table) {
  
  fam1.1 = lookup.table$motif.family[lookup.table$gene.symbol == name1] 
  fam1.2 = lookup.table$motif.family[lookup.table$gene.symbol == name2]
  fam2.1 = lookup.table$cluster.id[lookup.table$gene.symbol == name1]
  fam2.2 = lookup.table$cluster.id[lookup.table$gene.symbol == name2]
  fam3.1 = lookup.table$binding.domain[lookup.table$gene.symbol == name1]
  fam3.2 = lookup.table$binding.domain[lookup.table$gene.symbol == name2]
  
  l.cl <- length(intersect(fam1.1,fam1.2))
  l.jasp <- length(intersect(fam2.1,fam2.2))
  l.htfs <- length(intersect(fam3.1,fam3.2))
  
  return(ifelse(sum(l.cl, l.jasp, l.htfs)==0, F, T))
}


###############################################################################
#
# General Utitlity Functions
#
# #
# # # # # #
# # # # # # # # # #
# # # # # # # # # # # # # # # # # # # #
###############################################################################

rolling_mean <- rollify(mean, window = 7, na_value = 0)

binarize <- function(x, thres) {
  x <- ifelse(x < thres, 0, 1)
  return(x)
}

plen <- function(x) {
  return(pull(x) %>% length())
}

puniq <- function(x) {
  return(pull(x) %>% unique() %>% length())
}

zeroOne <- function(x){
  (x - min(x)) / (max(x) - min(x))
}

#' Generic function to map values using a mapping table
#' 
#' @param x a vector
#' @param mapping.table a named vector, where names correspond to the possible
#' values of x, and vector values correspond to their mapping (the output) 
mapValues <- function(
  x, 
  mapping.table
) {
  
  out <- mapping.table[x]
  
  if(any(is.na(out))) { 
    na <- sum(is.na(out)) 
    cat(sprintf("Warning: mapping returned %s NA values\n", as.character(na)))
  }
  
  return(out)
}

#' Clean a matrix for plotting. 
#' 
#' Removes NA, Inf, and NaN
#' Trims quantiles
#' Scales rows or columns
#' @param mat A numeric matrix
#' @param remove Character or char vector describing classes to remove
#' @param trim Whether to trim quantiles or not (default = TRUE)
#' @param cuts Quantile values to cut
#' @param scaling Character describing the style of scaling to employ
cleanMatrix <- function(
  mat, 
  remove = c("na", "inf", "nan"), 
  cuts = c(.1,.99), 
  scaling = c("row", "col", "both", "none")
) {
  
  remove = tolower(remove)
  scaling = tolower(scaling)
  
  if("na" %in% remove) {
    mat[which(is.na(mat))] <- 0
    cat("Removed NA\n")
  } 
  if("inf" %in% remove) {
    mat[mat==Inf] <- max(mat[!mat==Inf], na.rm = T)
    cat("Removed Inf\n")
  } 
  if("nan" %in% remove) {
    mat[is.nan(mat)] <- 0
    cat("Removed NaN\n")
  } 
  if(is.null(remove)) {
    mat <- mat
  }
  
  if(!is.null(cuts)) {
    mat <- trimQuantiles(mat, cuts = cuts)
  }
  
  if(scaling == "row"){
    mat <- rowScale(mat)
  }
  if(scaling == "col"){
    mat <- colScale(mat)
  }
  if(scaling == "both"){
    mat <- scale(mat)
  } 
  if(is.null(scaling) || scaling == "none"){
    mat <- mat
  }
  
  return(mat)
  
}

#' Trim the extreme quantiles of a numeric matrix or vector
#' Makes large heatmap plots more interpretable, since outlier values rescale color values
#' for the whole heatmap
#' 
#' @param mat a numeric matrix or vector to trim
#' @cuts a length-2 vector of quantiles, low and high, to trim from the object 
#' (they will be replaced by min and max respectively)
trimQuantiles <- function(x, 
                          cuts = c(.1,.99)) {
  
  if((
    any(is.null(x)) || 
    any(is.na(x)) ||
    any(is.infinite(x))
  )) {
    warning("Object contains non-numeric values.")
  }
  
  if(sum(cuts >= 1) == 2 & sum(cuts <= 100) ==2) {
    cuts <- cuts / 100
  }
  
  quantile.cuts <- quantile(x, cuts, na.rm = T)
  x[x < quantile.cuts[1]] <- quantile.cuts[1]
  x[x > quantile.cuts[2]] <- quantile.cuts[2]
  
  return(x)
}

#' View the first rows *and columns of of an object,
#' 
#' Handles matrices smartly.
#' Handles lists smartly, e.g., displays the first n elements
#' of the first list item, instead of printing the entirity of the first 
#' list elements.
#' @param X an atomic object
#' @param n the number of elements ^2 to display
med <- function(X, 
                n = 6) {
  
  if( !class(X) == "list" & (is.vector(X) || is.null(dim(X)) )) {
    len <- length(X)
    cl <- class(X)
    cat(sprintf("Vector \nlength %s\n, %s", len, cl))
    
    n <-  min(n, len)
    out <-  X[1:n]
    return(out)
  } 
  
  if(is.matrix(X) || is.data.frame(X) || length(dim(X)) == 2) {
    rows <- nrow(X)
    cols <- ncol(X)
    type <- ifelse(is.data.frame(X), "Data Frame", ifelse(is.matrix(X), "Matrix", "2D object"))
    cat(sprintf("%s dim %s x %s\n", type, rows, cols))
    
    n1 = min(n, rows)
    n2 = min(n, cols)
    
    out = X[1:n1, 1:n2]
    return(out)
  } 
  
  if(class(X) == "list") {
    
    len = length(X)
    tag = ifelse(len > n, "...", "")
    ln = min(n, len)
    
    n.names = names(X)[1:ln]
    cat(sprintf("List length %s\n", len))
    cat(sprintf("Names: %s%s\n", str_c(n.names, sep = ", "), tag))
    
    out <- lapply(X[1:ln], function(x) {
      return(med(x, n = n))
    })
    return(out)
    
  } else{
    cat("Invalid input ...\n")
  }
}

#' Scale/zscore a vector of values. Compatible with dplyr::group_by() 
#' 
#' @param x a vector of numeric values
scale_this <- function(x) {
  stdev <- sd(x, na.rm=TRUE)
  if(stdev == 0) {
    return(0)
  } else {
    return((x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE))
  }
}

#' Split a character vector and take the desired piece
#' 
#' @param x a character vector
#' @param split.char A character / pattern to split each string on
#' @param idx the index of the desired piece
#' @example 
#' vec <- c("M10284_Fox", "M10285_Ram")
#' add_split_name_vector(vec, "_", 2)
#' [1] "Fox", [2] "Ram"
add_split_name_vector <- function(x, 
                                  split.char, 
                                  idx) {
  split.names <- strsplit(as.character(x), split.char)
  y <- as.character(sapply(split.names, "[", idx))
  return(y)
}

#' Scale rows of a matrix (row z scores)
#'
#' @param mat a numeric matrix
rowScale <- function(mat) {
  
  if (is.null(nrow(mat))) {
    return(scale_this(mat))
  }
  
  rn = rownames(mat)
  cn = colnames(mat)
  
  out = do.call(rbind, lapply(1:nrow(mat), function(x) {
    s = scale_this(mat[x, ])
  }))
  
  rownames(out) = rn
  colnames(out) = cn
  
  return(out)
}

#' Scale columns of a matrix (column z scores)
#'
#' @param mat a numeric matrix
colScale <- function(mat) {
  
  rn = rownames(mat)
  cn = colnames(mat)
  
  out = do.call(cbind, lapply(1:ncol(mat), function(x) {
    s = scale_this(mat[ , x])
  }))
  
  rownames(out) = rn
  colnames(out) = cn
  
  return(out)
}

#' Randomly sample, preserving original order
orderedRandomSample <- function(values, N){
  
  size <- length(values)
  
  values[sapply(1:size, function(i){
    select <- as.logical(rbinom(1,1,N/(size+1-i)))
    if(select) N <<- N - 1
    select
  })]
  
}

#' Randomly matrix columns or rows
#' 
#' @param mat A matrix of numeric values
#' @param MARGIN 1 indicating rows, 2 indicating columns
#' @param subset.size The size of the sample. If less than one, returns a fraction
subsetMatrix <- function(mat, MARGIN = 2, subset.size = 10000) {
  if(MARGIN==2) {
    samp <- orderedRandomSample(1:ncol(mat), subset.size)
    new.mat <- mat[ ,samp]
  } else{
    samp <- orderedRandomSample(1:nrow(mat), subset.size)
    new.mat <- mat[samp , ]
  }
  return(new.mat)
}


#' This doesn't work yet.
#' Diagonal reording of matrix columns. 
#' 
#' @param mat A matrix of values
#' @param element_order indicates the desired order of the matrix rows. If NULL, then it defaults to the existing row order.
#' @param subset_matrix Integer value indicating whether to subset columns.
#' @param top_left_first indicates the diagonal directionality. If TRUE (default), high values will run from top left to bottom right.
#' The default ordering function is a weighted mean, with column values as the vector and the row index as the weights
#' 
reorderMatrixDiagonally <- function(mat, 
                                    element_order, 
                                    subset_matrix = NULL, 
                                    top_left_first = T) {
  
  if(is.null(subset_matrix)) {
    mat <- mat
  } else{
    mat <- mat[, sample(x = 1:ncol(mat), size = subset_matrix, replace = F)]
  }
  
  
  if(is.null(element_order)) {
    element_order <- 1:nrow(mat)
  }
  
  
  
  ord.mat <- mat[element_order, ]
  wm = apply(ord.mat,
             MARGIN = 2, 
             FUN = weighted.mean, 
             w = 1:nrow(ord.mat))
  peak.order = order(wm)
  
  if(top_left_first == F) {
    peak.order = rev(peak.order)
  }
  
  out.mat <-  ord.mat[, peak.order]
  return(out.mat)
  
}

#' Really great memory management tool. Quickly list the largest memory
#' objects in your environment. 
#' 
#' Remove them with rm() or remove() then garbage collect "gc()"
#' 
.ls.objects <- function (pos = 1, 
                         pattern, 
                         order.by,
                         decreasing=FALSE, 
                         head=FALSE, 
                         n=5) {
  napply <- function(names, fn) sapply(names, function(x)
    fn(get(x, pos = pos)))
  names <- ls(pos = pos, pattern = pattern)
  obj.class <- napply(names, function(x) as.character(class(x))[1])
  obj.mode <- napply(names, mode)
  obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
  obj.size <- napply(names, object.size)
  obj.dim <- t(napply(names, function(x)
    as.numeric(dim(x))[1:2]))
  vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
  obj.dim[vec, 1] <- napply(names, length)[vec]
  out <- data.frame(obj.type, obj.size, obj.dim)
  names(out) <- c("Type", "Size", "Rows", "Columns")
  if (!missing(order.by))
    out <- out[order(out[[order.by]], decreasing=decreasing), ]
  if (head)
    out <- head(out, n)
  out
}
lsos <- function(..., n=10) {
  .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}



###############################################################################
#' Sample pseudobulks from matrix and compute differential counts (DESeq2)
#' 
#' From Fabian Mueller. FM
#' 
#' @param X a matrix 
#' @param idx1 a vector of indices for cells
#' @param idx2 a vector of indices for cells (Group 2)
#' 
scDifferential_pseudobulk <- function(X,
                                      idx1, 
                                      idx2, 
                                      cm.g1 = NULL,
                                      cm.g2 = NULL,
                                      nCellsPerBulk = 100, 
                                      nBulkPerGroup = 20, 
                                      rowranges = NULL){
  
  if (is.null(cm.g1)) {
    nCells <- c(length(idx1), length(idx2))
    names(nCells) <- paste0("group", 1:2)
    nCellsPerBulk_cur <- nCellsPerBulk
    
    doBoostrap <- FALSE
    
    if (any(nCells < nCellsPerBulk)){
      nCellsPerBulk_cur <- min(c(nCells, nCellsPerBulk))
      logger.warning(c("Few cells detected per group", "--> selecting only", nCellsPerBulk_cur, "cells for sampling"))
      doBoostrap <- TRUE
    }
  
    idx_sampled_g1 <- lapply(1:nBulkPerGroup, FUN = function(i){
      sample(idx1, nCellsPerBulk_cur, replace = doBoostrap)
    })
    
    idx_sampled_g2 <- lapply(1:nBulkPerGroup, FUN = function(i){
      sample(idx2, nCellsPerBulk_cur, replace = doBoostrap)
    })
    
    cm.g1 <- do.call("cbind", lapply(idx_sampled_g1, FUN = function(iis){
      rowSums(X[,iis])
    }))
    
    cm.g2 <- do.call("cbind", lapply(idx_sampled_g2, FUN=function(iis){
      rowSums(X[,iis])
    }))
    
    colnames(cm.g1) <- paste("sample", "g1", 1:nBulkPerGroup, sep = "_")
    colnames(cm.g2) <- paste("sample", "g2", 1:nBulkPerGroup, sep = "_")
    
  }
  
  
  cat('Generating DESeq2 Dataset\n')
  
  print(dim(cm.g1))
  print(dim(cm.g2))
  
  col.df <- data.frame(
    sampleId = c(colnames(cm.g1), colnames(cm.g2)), 
    group    = rep(c("group1", "group2"), times = c(ncol(cm.g1), ncol(cm.g2)))
  )
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = cbind(cm.g1, cm.g2),
    colData = col.df,
    design = as.formula(paste0("~", paste("group", collapse = "+")))
  )
  
  if (!is.null(rowranges)) {
    SummarizedExperiment::rowRanges(dds) <- rowranges
  }
  
  dds <- DESeq2::DESeq(dds)
  dm <- data.frame(DESeq2::results(dds, contrast = c("group", "group1", "group2")))
  
  cat('Preparing results \n')
  rankMat <- cbind(
    # rank(-dm[,"baseMean"]), na.last="keep", ties.method="min"),
    rank(-abs(dm[,"log2FoldChange"]), na.last = "keep", ties.method = "min"),
    rank(dm[,"pvalue"], na.last = "keep", ties.method = "min")
  )
  
  dm[,"cRank"] <- matrixStats::rowMaxs(rankMat, na.rm = FALSE)
  # dm[,"cRank"] <- rowMaxs(rankMat, na.rm=TRUE)
  dm[!is.finite(dm[,"cRank"]),"cRank"] <- NA
  dm[,"cRank_rerank"] <- rank(dm[,"cRank"], na.last = "keep", ties.method = "min")
  
  l2fpkm <- log2(DESeq2::fpkm(dds, robust = TRUE)+1)
  grp1.m.l2fpkm <- rowMeans(l2fpkm[, colnames(cm.g1), drop = FALSE], na.rm = TRUE)
  grp2.m.l2fpkm <- rowMeans(l2fpkm[, colnames(cm.g2), drop = FALSE], na.rm = TRUE)
  
  vstCounts <- SummarizedExperiment::assay(DESeq2::vst(dds, blind = FALSE))
  grp1.m.vst <- rowMeans(vstCounts[, colnames(cm.g1), drop = FALSE], na.rm = TRUE)
  grp2.m.vst <- rowMeans(vstCounts[, colnames(cm.g2), drop = FALSE], na.rm = TRUE)
  
  diffTab <- data.frame(
    geneId          = rownames(dm),
    geneName        = elementMetadata(rowranges[rownames(dm)])[,"gene_name"],
    log2BaseMean    = log2(dm[,"baseMean"] + 1e-4),
    meanLog2Fpkm_g1 = grp1.m.l2fpkm,
    meanLog2Fpkm_g2 = grp2.m.l2fpkm,
    meanVstCount_g1 = grp1.m.vst,
    meanVstCount_g2 = grp2.m.vst,
    dm,
    stringsAsFactors = FALSE
  )
  
  res <- list(
    # nCells=nCells,
    dataset=dds,
    tab=diffTab
  )
  return(res)
  
}

###############################################################################
#
# External Functions
#
# #
# # # # # #
# # # # # # # # # #
# # # # # # # # # # # # # # # # # # # #
###############################################################################

# Some packages are pretty close to what you need, but I needed to monkey 
# around with them more than they would allow. Here is the source code for
# mfuzz functions copied in. 

###############################################################################
# Package: Mfuzz 
###############################################################################

#' A function to compute overlaps between fuzzy clusters
#' 
#' Output is a N.clusters x N.clusters matrix. Each element i, j is the sum of 
#' the vector product U[ i] * U[ ,j], where U is the membership matrix. This 
#' data is then column scaled. This is adapted from the Mfuzz package
computeFuzzyOverlaps <- function(U) {
  
  O <- matrix(
    data = NA, 
    nrow=dim(U)[2],
    ncol=dim(U)[2]
  )
  
  for (i in 1:dim(U)[2]) {
    for (j in 1:dim(U)[2]){
      O[i,j] <-  sum(U[,i] * U[,j])
    }
  }
  for (i in 1:dim(O)[2]){
    O[,i] <- O[,i] / sum(O[,i])
  }
  return(O)
}


overlap <- function(cl) {
  
  O <- matrix(
    data = NA, 
    nrow=dim(cl[[4]])[2],
    ncol=dim(cl[[4]])[2]
  )
  
  for (i in 1:dim(cl[[4]])[2]) {
    for (j in 1:dim(cl[[4]])[2]){
      O[i,j] <-  sum(cl[[4]][,i] * cl[[4]][,j])
    }
  }
  for (i in 1:dim(O)[[2]]){
    O[,i] <- O[,i] / sum(O[,i])
  }
  O
}

overlap.plot <- function(cl,
                         overlap,
                         thres = 0.1,
                         scale = TRUE, 
                         magni = 30,
                         P = NULL) {
  
  x <- prcomp(cl[[1]], scale=TRUE)
  if (!missing(P)){
    x[[2]] <- P
  }
  
  x[[5]] <- t(t(x[[2]]) %*% t(cl[[1]]))
  plot(
    x[[5]][,1],
    x[[5]][,2],
    xlab = "PC1",
    ylab = "PC2",
    type = "n"
  )
  
  for (i in 1:dim(x[[5]])[[1]]){
    for (j in 1:dim(x[[5]])[[1]]){
      if (thres < overlap[i,j]){
        lines(
          c(x[[5]][i,1],x[[5]][j,1]),
          c(x[[5]][i,2],x[[5]][j,2]),
          col = "blue",
          lwd = magni * overlap[i,j]
        )
      }
    }
  }
  
  for (i in 1:dim(x[[5]])[[1]]){
    points(x[[5]][i,1],x[[5]][i,2], pch=20, cex=4, col="red", lwd=2)
    points(x[[5]][i,1],x[[5]][i,2], pch=21, cex=4, col="red", lwd=2)
    text(x[[5]][i,1],x[[5]][i,2], i, font=2)
  }
  
  x[[2]]
}


plotFuzzyOverlap <- function(prototype,
                             overlap,
                             thres = 0.1,
                             scale = TRUE, 
                             magni = 30, 
                             P = NULL) {
  
  x <- prcomp(cl[[1]], scale=TRUE)
  if (!missing(P)){
    x[[2]] <- P
  }
  
  x[[5]] <- t(t(x[[2]]) %*% t(cl[[1]]))
  plot(
    x[[5]][,1],
    x[[5]][,2],
    xlab = "PC1",
    ylab = "PC2",
    type = "n"
  )
  
  for (i in 1:dim(x[[5]])[[1]]){
    for (j in 1:dim(x[[5]])[[1]]){
      if (thres < overlap[i,j]){
        lines(
          c(x[[5]][i,1],x[[5]][j,1]),
          c(x[[5]][i,2],x[[5]][j,2]),
          col = "blue",
          lwd = magni * overlap[i,j]
        )
      }
    }
  }
  
  for (i in 1:dim(x[[5]])[[1]]){
    points(x[[5]][i,1],x[[5]][i,2], pch=20, cex=4, col="red", lwd=2)
    points(x[[5]][i,1],x[[5]][i,2], pch=21, cex=4, col="red", lwd=2)
    text(x[[5]][i,1],x[[5]][i,2], i, font=2)
  }
  
  x[[2]]
}

mfuzz.plot <- function(mat,
                       cl,
                       mfrow=c(1,1),
                       colo,
                       min.mem = 0,
                       time.labels,
                       new.window = TRUE) {
  
  # function for plotting the clusters 
  clusterindex <- cl[[3]]
  memship <- cl[[4]]
  memship[memship < min.mem] <- -1 
  colorindex <- integer(dim(mat)[[1]])
  if (missing(colo)){
    colo <- c("#FF8F00",
              "#FFA700", "#FFBF00", "#FFD700", "#FFEF00", "#F7FF00", "#DFFF00", "#C7FF00",
              "#AFFF00", "#97FF00", "#80FF00", "#68FF00", "#50FF00", "#38FF00", "#20FF00",
              "#08FF00", "#00FF10", "#00FF28", "#00FF40", "#00FF58", "#00FF70", "#00FF87",
              "#00FF9F", "#00FFB7", "#00FFCF", "#00FFE7", "#00FFFF", "#00E7FF", "#00CFFF",
              "#00B7FF", "#009FFF", "#0087FF", "#0070FF", "#0058FF", "#0040FF", "#0028FF",
              "#0010FF", "#0800FF", "#2000FF", "#3800FF", "#5000FF", "#6800FF", "#8000FF",
              "#9700FF", "#AF00FF", "#C700FF", "#DF00FF", "#F700FF", "#FF00EF", "#FF00D7",
              "#FF00BF", "#FF00A7", "#FF008F", "#FF0078", "#FF0060", "#FF0048", "#FF0030",
              "#FF0018")
    
  }
  
  colorseq <- seq(0, 1, length = length(colo))
  
  for (j in 1:max(clusterindex)){
    tmp <- mat[clusterindex==j, , drop=FALSE]# thanks Ian for the fix
    tmpmem <- memship[clusterindex==j,j]
    
    if (((j-1)%% (mfrow[1] * mfrow[2]))==0){
      
      if (new.window) X11()
      par(mfrow=mfrow)
      
      if (sum(clusterindex==j)==0) {
        ymin <- -1; ymax <- +1;
      } else {
        ymin <- min(tmp);ymax <- max(tmp);    
      }
      
      plot.default(x=NA,xlim=c(1,dim(exprs(eset))[[2]]), ylim= c(ymin,ymax),
                   xlab="Time",ylab="Expression changes",main=paste("Cluster",j),axes=FALSE)
      if (missing(time.labels)){
        axis(1, 1:dim(exprs(eset))[[2]],c(1:dim(exprs(eset))[[2]]))
        axis(2)
      } else {
        axis(1, 1:dim(exprs(eset))[[2]],time.labels)
        axis(2)
      } 
      
      
    } else {
      
      if (sum(clusterindex==j)==0) {
        ymin <- -1; ymax <- +1;
      } else {
        ymin <- min(tmp);ymax <- max(tmp);    
      }
      
      
      plot.default(x=NA,xlim=c(1,dim(exprs(eset))[[2]]), ylim= c(ymin,ymax),
                   xlab="Time",ylab="Expression changes",main=paste("Cluster",j),axes=FALSE)
      
      if (missing(time.labels)){
        axis(1, 1:dim(exprs(eset))[[2]],c(1:dim(exprs(eset))[[2]]))
        axis(2)
      } else {
        axis(1, 1:dim(exprs(eset))[[2]],time.labels)
        axis(2)
      } 
      
      
    }
    
    
    if (!(sum(clusterindex==j)==0)){
      for (jj in 1:(length(colorseq)-1)){
        tmpcol <- (tmpmem >= colorseq[jj] & tmpmem <= colorseq[jj+1])
        if (sum(tmpcol)> 0) {
          tmpind <- which(tmpcol)
          for (k in 1:length(tmpind)){
            lines(tmp[tmpind[k],],col=colo[jj])
          }
        }
      }}
  }
}