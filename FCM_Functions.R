###############################################################################
#
#   Fuzzy C-means (FCM) Clustering Functions
#
#   Written or adapted by Alexandro E. Trevino 2017-2020 
#   Stanford University
#   special mention to Fabian Mueller as the source of some key functions
# 
###############################################################################



###############################################################################
#' Set paramters for fuzzy clustering
#' 
#' @param c The number of centers
#' @param m The fuzziness parameter
#' @param object an FCM project object
setFuzzyParameters <- function(c, m, object) {
  object[['Parameters']] <- list(centers = c, m = m)
  return(object)
}


###############################################################################
#' Get centers
getNumCenters <- function(object) {
  return(object[['Parameters']][['centers']])
}


###############################################################################
#' Get fuzziness parameter
getM <- function(object) {
  return(object[['Parameters']][['m']])
}


###############################################################################
#' Get data matrix from FCM object
#' 
#' @param object
getDataMatrix <- function(object) {
  if (is.null(object[['Data.Matrix']])) {
    stop("No data matrix!")
  }
  return(object[['Data.Matrix']])
}

###############################################################################
#' Function to perform fuzzy C-means clustering. This is the e1071 
#' implementation.
#' 
#' @param object An FCM object
#' @param iter.max The max number of iterations
#' @param verbose To print iteration info
#' @param dist The distance metric used
#' @param method Cmeans
doFuzzyClustering <- function(object,
                              iter.max = 30, 
                              verbose = T, 
                              dist = 'euclidean', 
                              method = 'cmeans') {
  
  if (is.null(object[['Parameters']])) {
    stop("Set parameters first")
  }
  
  if (is.null(object[['Data.Matrix']])) {
    stop("No data matrix")
  }
  
  mat <- getDataMatrix(object)
  c <- getNumCenters(object)
  m <- getM(object)
  
  fuzz <- e1071::cmeans(
    x = mat, 
    centers = c, 
    iter.max = iter.max, 
    verbose = verbose, 
    dist = dist, 
    method = method, 
    m = m
  )
  
  object[['FCM.Result']] <- fuzz
  
  return(object)
}


###############################################################################
getFclustObject <- function(object) {
  if (is.null(object[['FCM.Result']])) {
    stop('No FCM result found. Have you run doFuzzyClustering?')
  }
  
  return(object[['FCM.Result']])
}


###############################################################################
#' Get centers from FCM project
#' 
#' @param object An FCM object
getCenterMatrix <- function(object) {
  return(object[['FCM.Result']][['centers']])
}


###############################################################################
#' Get membership matrix from FCM project
#' 
#' @param object An FCM object
getMembership <- function(object) {
  return(object[['FCM.Result']][['membership']])
}


###############################################################################
#' Get trial/project name from FCM project
#' 
#' @param object An FCM object
getTrialName <- function(object) {
  return(object[['Trial.Name']])
}


###############################################################################
#' Plot a heatmap of gene membership by FCM module
#' 
#' Also plots the corresponding heatmap of binarized membership, given a 
#' threshold value
#' 
#' @param object A FCM project object
plotMembershipHeatmap <- function(object, threshold = NULL) {
  
  fuzz.membership <- getMembership(object)
  trial.name      <- getTrialName(object)
  
  if (is.null(threshold)) {
    threshold <- getThreshold(object)
  }
  
  threshold.char  <- as.character(threshold)
  
  ms1 <- Heatmap(
    matrix = t(fuzz.membership), 
    col = plant_palettes$`Jasminium polyanthum A`,
    show_column_names = F,
    name = "Membership"
  )
  
  ms2 <- Heatmap(
    matrix = t(binarize(fuzz.membership, threshold)), 
    col = c("white", "black"),
    show_column_names = F,
    column_order = column_order(ms1),
    cluster_columns = column_dend(ms1),
    name = sprintf("Threshold %s", threshold.char)
  )
  
  pdf.name <- sprintf(
    "fuzzy_clusters/%s/MembershipHeatmap_thres%s.pdf", 
    trial.name, 
    threshold.char
  )
  
  pdf(
    file = pdf.name,
    width = 9, 
    height = 4
  )
  patch <- ms1 + ms2
  print(patch)
  dev.off()
  
}


plotNumGenesInNModules <- function(object) {
  
  m.matrix   <- getMembership(object)
  threshold  <- getThreshold(object)
  trial.name <- getTrialName(object)
  
  df <- apply(
    X      = binarize(m.matrix, threshold), 
    MARGIN = 1,
    FUN    ='sum'
  )
  
  df <- data.frame(n = factor(df))
  gg <- ggplot(df, aes(x = n)) +
    geom_bar(color = 'black', fill = 'grey') +
    theme_classic() +
    theme(aspect.ratio = 1,
          title = element_text(hjust = 0.5)) +
    labs(title = "Num. genes belonging to n modules")
  
  pdf.name <- sprintf(
    'fuzzy_clusters/%s/Genes_Belonging_To_N_Modules_thres%s.pdf',
    trial.name,
    as.character(threshold)
  )
  pdf(
    file = pdf.name,
    width = 4.5,
    height = 4.5
  )
  print(gg)
  dev.off()
}
# apply(binarize(getMembership(fuzz), getThreshold(fuzz)), 1, 'sum')


###############################################################################
#' Get a list of genes belonging to each fuzzy cluster / module
#' 
#' Expects ensemble gene IDs
#' @param object A FCM project object
addMemberGeneList <- function(object) {
  
  fuzz.membership <- getMembership(object)
  threshold       <- getThreshold(object)
  
  fcm.membership <- apply(
    X = fuzz.membership, 
    MARGIN = 2, 
    FUN = function(x) {
      convertGeneIDs(rownames(fuzz.membership)[x > threshold], ensembl2name)
    }
  )
  
  object[['Member.Genes']] <- fcm.membership
  return(object)
}


###############################################################################
getMemberGeneList <- function(object) {
  return(object[['Member.Genes']])
}


###############################################################################
getMemberGeneDF <- function(object) {
  
  fcm <- object[['Member.Genes']]
  
  df <- do.call(
    rbind,
    lapply(
      seq_along(fcm),
      function(x) {
        mname <- paste0('m', as.character(x))
        genes <- fcm[[x]]
        ids <- names(fcm[[x]])
        out <- data.frame(
          cbind(
            module = mname,
            gene.symbol = genes,
            gene.id = ids
          )
        )
        return(out)
      }
    )
  )
  
  df$module <- factor(
    x      = df$module,
    levels = paste0('m', as.character(names(sort(object[['Module.Pseudotime']]))))
  )
  return(df)
}

###############################################################################
addMemberGeneWeights <- function(object) {
  
  fcm.membership <- getMemberGeneList(object)
  fuzz.membership <- getMembership(object)
  
  fcm.weights <- lapply(
    seq_along(fcm.membership),
    function(x) {
      genes <-  names(fcm.membership[[x]])
      memberships <- fuzz.membership[genes, x]
      return(memberships)
    }
  )
  
  object[['Member.Gene.Weights']] <- fcm.weights
  return(object)
}


###############################################################################
getMemberGeneWeights <- function(object) {
  return(object[['Member.Gene.Weights']])
}


###############################################################################
#' Returns the current threshold
#' 
#' @param object
getThreshold <- function(object) {
  if (is.null(object[['Threshold']])) {
    stop("Threshold has not been set")
  } else {
    return(object[['Threshold']])
  }
}


###############################################################################
#' Sets the current threshold
#' 
#' @param object
setThreshold <- function(threshold, object) {
  
  maxm <- maxMemberhip(object, threshold = threshold)
  
  if (any(maxm < threshold)) {
    print('Empty clusters found')
  }
  
  object[['Threshold']] <- threshold
  
  return(object)
}


###############################################################################
#' Computes overlaps from a fuzzy clustering object
#' 
#' If n is the number of centers, the result will be an n x n matrix describing
#' the extent of overlap between the clusters. N[i,j] != N[j, i]. This is based
#' on ratios by default. However, Jaccard and scaled Jaccard indices are also
#' available (in dev)
#' 
#' @param object An FCM object
#' @param method If 'ratio' (default) then the output is the fraction of genes
#' in cluster C_i that overlap with cluster C_j. If 'number', the raw number
#' of genes is returned. If 'jaccard', a Jaccard index is computed. If
#' 'scaled', then a scaled Jaccard is returned
#' @pre.filter An expression which will be evaluated and passed to lapply() as
#' the argument FUN = . 
computeOverlaps <- function(object, method = 'ratio', pre.filter = NULL) {
  
  fcm.obj <- getFclustObject(object)
  
  if (method == 'ratio') {
    overlaps <- overlap(fcm.obj)
  } else if (method == 'number') {
    
    if (is.null(getMemberGeneList(object))) {
      stop('No gene list added. Run addMemberGeneList()')
    }
    
    fcm.membership <- getMemberGeneList(object)
    
    overlaps <- lapply(
      fcm.membership,
      function(i) {
        lapply(
          fcm.membership,
          function(j) {
            length(intersect(i, j))
          }
        ) %>% Reduce('cbind', .)
      }
    ) %>% Reduce('rbind', .)
    
    colnames(overlaps) <- rownames(overlaps) <- names(fcm.membership)
    
  } else if (method == 'jaccard' || method == 'scaled' || method == 'filter') {
    
    fcm.membership <- getMemberGeneList(object)
    
    if (!is.null(pre.filter)) {
      fcm.membership <- lapply(fcm.membership, pre.filter)
    }
    
    overlaps <- lapply(
      fcm.membership,
      function(i) {
        lapply(
          fcm.membership,
          function(j) {
            length(intersect(i, j)) / length(union(i, j))
          }
        ) %>% Reduce('cbind', .)
      }
    ) %>% Reduce('rbind', .)
    
    colnames(overlaps) <- rownames(overlaps) <- names(fcm.membership)
    
    if (method == 'scaled') {
      overlaps <- zeroOne(overlaps)
    }
    
  } else {
    stop("Not a valid overlap method")
  }
  
  object[['Overlaps']][[method]] <- overlaps
  return(object)
}


###############################################################################
#' Get overlaps computed from FCM object
#' 
#' @param object An FCM object
getOverlaps <- function(object, method = 'ratio') {
  return(object[['Overlaps']][[method]])
}


###############################################################################
#' Plot connections between fuzzy modules, using the Mfuzz ratio product 
#' overlaps, and using the Mfuzz overlap PCA plotting method
#' 
#' @param object An FCM object
#' @param threshold Whether to override the current threshold in the FCM object
plotOverlapPCA <- function(object, threshold = NULL, method = 'ratio') {
  
  if (is.null(getOverlaps(object, method = 'ratio'))) {
    stop('No overlap product matrix computed. Run computeOverlaps first.')
  }
  
  if (is.null(threshold)) {
    threshold <- getThreshold(object)
  } else {
    threshold <- threshold
  }
  
  pdf.name <- sprintf(
    "fuzzy_clusters/%s/OverlapPCA_thres%s.pdf", 
    getTrialName(object), 
    as.character(threshold)
  )
  
  pdf(
    file = pdf.name,
    height = 5,
    width = 5,
    useDingbats = F
  )
  overlap.plot(
    cl = getFclustObject(object), 
    overlap = getOverlaps(object, method = method), 
    thres = threshold
  )
  
  dev.off()
}


###############################################################################
#' Plot the number of genes in each module, and the number of genes belonging
#' to any module, at a given threshold. Helps find the membership value 
#' threshold such that every gene is included downstream, or not.
#' 
#' @param object An FCM project
#' @param threshold Whether to override the current threshold in the FCM object
plotGenesAtThreshold <- function(object, threshold = NULL) {
  
  if (is.null(threshold)) {
    threshold <- getThreshold(object)
  } else {
    threshold <- threshold
  }
  
  fuzz.membership <- getMembership(object)
  trial.name      <- getTrialName(object)
  
  max.membership <- sort(apply(fuzz.membership, 1, function(x) { max(x) }))
  genes.per.module <- apply(fuzz.membership, 2, function(x) { sum(x > threshold) })
  
  pdf.name <- sprintf(
    "fuzzy_clusters/%s/GenesAtThreshold_%s.pdf",
    getTrialName(object), 
    as.character(threshold)
  )
  pdf(pdf.name, width = 8, height = 4)
  
  par(mfrow = c(1,2))
  
  plot(max.membership, ylim = c(0,1))
  abline(h = threshold)
  barplot(genes.per.module) 
  dev.off()
  
  return(NULL)
}


###############################################################################
#' Plot GO enrichments 
#' 
#' @param object An FCM object
#' @param max.terms The max number of Go terms to plot per cluster
plotGOfcm <- function(object, 
                      ontology = 'MF', 
                      max.terms = 6, 
                      color.pal = viridis_pal(option = 'E')(100),
                      ncol = 1,
                      width = 7,
                      height = 12,
                      sel.clusters = NULL,
                      facet = TRUE,
                      make.file = TRUE) {
  
  fcm.go <- getGOEnrichments(object, ontology = ontology)
  
  # Subset to top enrichments and compute odds ratio
  
  gg.fcm.go <- fcm.go@compareClusterResult %>%
    rowwise() %>%
    mutate(gene.ratio = eval(parse(text = GeneRatio)), 
           bg.ratio = eval(parse(text = BgRatio)),
           Fold.Enrichment = gene.ratio / bg.ratio) %>%
    group_by(Cluster) %>%
    slice_head(n = max.terms) 
  
  if (! is.null(sel.clusters)) {
    gg.fcm.go %<>% filter(Cluster %in% sel.clusters)
  }
  
  go.gg <- ggplot(gg.fcm.go, aes(x = -log10(p.adjust), y = Description, fill = log2(Fold.Enrichment))) +
    geom_bar(stat = 'identity') +
    theme_classic() +
    scale_fill_gradientn(colors = color.pal, limits = c(0,NA)) +
    scale_y_discrete(position = 'right')
  
  go.gg2 <- ggplot(gg.fcm.go, aes(x = -log10(p.adjust), y = Description, fill = log2(Fold.Enrichment), size = Count)) +
    geom_point(shape = 21) +
    theme_classic() +
    
    scale_fill_viridis(option = 'A')
  
  if (facet == T) {
    go.gg <- go.gg + facet_wrap(~ Cluster, ncol = ncol, scales = "free_y", strip.position = 'left')
    go.gg2 <- go.gg + facet_wrap(~ Cluster, ncol = ncol, scales = 'free_y')
  }
  

  if (make.file == T) {
    pdf.name <- sprintf(
      "fuzzy_clusters/%s/GeneOntologies_%s_thres%s.pdf", 
      getTrialName(object), 
      ontology,
      as.character(getThreshold(object))
    )
    
    pdf(pdf.name, width = width, height = height, useDingbats = F)
    print(go.gg)
    print(go.gg2)
    dev.off()
  } else {
    print(go.gg)
    print(go.gg2)
  }
  
}


#' Add GO enrichments to FCM object.
#' 
#' Uses the gene memberships computed by getMemberGeneList
#' @param object An FCM object
#' @param fun clusterProfiler parameter, which function to use for terms
#' @param org.db clusterProfiler parameter, OrgDb to use for terms
#' @param key.type clusterProfiler parameter, SYMBOL or ENSEMBL etc.
#' @param ontology clusterProfiler parameter, BP/MF/CC
#' @param universe clusterProfiler parameter, the background of genes to use
addGOEnrichments <- function(object,
                             fun = 'enrichGO', 
                             org.db = org.Hs.eg.db, 
                             key.type = 'SYMBOL',
                             ontology = 'MF',
                             universe = NULL) {
  
  fcm.membership <- getMemberGeneList(object)
  
  if (!is.null(universe)) {
    universe <- do.call(c, fcm.membership) %>% unique()
  }
  if (is.null(object[['Ontology.Enrichments']])) {
    object[['Ontology.Enrichments']] <- list()
  }
  
  fcm.go <- compareCluster(
    geneClusters = fcm.membership, 
    fun = fun, 
    OrgDb = org.db, 
    keyType = key.type,
    universe = universe
  )
  
  object[['Ontology.Enrichments']][[ontology]] <- fcm.go
  return(object)
}


###############################################################################
#' Get GO enrichments from FCM object
#' 
#' @param object
#' @param ontology The ontology name
getGOEnrichments <- function(object, ontology = 'MF') {
  
  if (is.null(object[['Ontology.Enrichments']][[ontology]])) {
    stop("No ontology enrichments computed. Run addGOEnrichments")
  }
  
  return(object[['Ontology.Enrichments']][[ontology]])
  
}


###############################################################################
#' This function plots the module gene expression in a projection of choice,
#' from an expression dataset of choice
#' 
#' @param object An FCM object
#' @param plot.object A matrix, SCE, or otherwise expression data set
#' @param plot.object.class A character vector specifying what class dataset it
#' is
#' @param objection A projection, e.g. the coordinates of a UMAP, with rows 
#' corresponding to the columns of the expression object
#' @param objection.name A name for the projection
#' @param dir.name The output directory name. 
#' @param weight.expression Whether or not to weight the gene expression by 
#' cluster membership
#' @param scale.expression Whether or not to scale expression prior to plotting
#' @param ... Additional arguments passed to geneSetAveragePlot()
plotModuleExpression <- function(object, 
                                 plot.object = sce, 
                                 plot.object.class = 'sce',
                                 proj.name = 'SCE',
                                 dir.name = NULL,
                                 weight.expression = F,
                                 scale.expression = T,
                                 ...) {
  
  trial.name <- getTrialName(object)
  threshold.char <- as.character(getThreshold(object))
  
  if (is.null(dir.name)) {
    dir.name <- sprintf(
      "fuzzy_clusters/%s/ModuleExpression_%s_thres%s", 
      trial.name,
      proj.name,
      threshold.char
    )
  }
  
  dir.create(dir.name, showWarnings = F)
  
  fcm.membership <- getMemberGeneList(object)
  if (weight.expression == T) {
    fcm.weights <- getMemberGeneWeights(object)
  } else {
    fcm.weights <- NULL
  }
  
  projection <- object[['Projections']][[proj.name]][['DF']]
  
  lapply(
    seq_along(fcm.membership),
    function(X) {
      
      genes <- fcm.membership[[X]]
      weights <- fcm.weights[[X]]
      module <- names(fcm.membership)[X]
      
      cat(sprintf("Plotting module %s\n", module))
      
      pdf.name <- sprintf(
        "%s/%s_Module_%s.pdf", 
        dir.name, 
        proj.name, 
        module
      )
      
      pdf(file = pdf.name, useDingbats = F)
      geneSetAveragePlot(
        genes = genes, 
        w = weights,
        object = plot.object, 
        object.class = plot.object.class,
        projection = projection,
        scaled = scale.expression,
        plot.type = 'averages',
        plot.title = sprintf("Fuzzy Module %s Expression", module),
        ...
      )
      dev.off()
    }
  )
}


###############################################################################
#' For each FCM gene module, compute the weighted mean of gene-wise pseudotimes
#' in the module. The weights are gene memberships in the module. Pseudotime 
#' is first computed on pseudobulks if not already done
#' 
#' @param object An FCM object
computeModulePseudotime <- function(object) {
  
  if (is.null(object[['sc.pseudotime']])) {
    stop("No pseudotime added")
  }
  
  if (is.null(object[['Pseudobulks']])) {
    stop("No pseudobulk information")
  }
  
  pseudotime  <- object[['sc.pseudotime']]
  pseudobulks <- object[['Pseudobulks']][['sampleIdxL']]
  
  # Calculate bulk pseudotimes (from pseudobulk indices and SC pseudotime)
  if (is.null(object[['pb.pseudotime']])) {
    pb.pseudotimes <- lapply(
      pseudobulks,
      function(x) {
        mean(pseudotime[x])
      }
    ) %>% Reduce(c, .) %>% zeroOne(.)
  } else {
    pb.pseudotimes <- object[['pb.pseudotime']]
  }
  
  # Get relevant objects
  fcm.membership <- getMemberGeneList(object)
  fcm.weights <- getMemberGeneWeights(object)
  mat <- getDataMatrix(object)
  
  module.pseudotime <- lapply(
    seq_along(fcm.membership),
    function(x) {
      
      genes <- names(fcm.membership[[x]])
      weights <- fcm.weights[[x]]
      
      mpt <- genePseudotimes(
        genes = genes, 
        object = logcounts(sce), 
        pseudotimes = pseudotime,
        cell.subset = rownames(object[['sc.colData']]),
        id.type = 'ensembl'
      )
      
      out <- weighted.mean(mpt, w = weights)
      return(out)
    }
  )
  
  module.pseudotime <- Reduce(c, module.pseudotime)
  names(module.pseudotime) <- names(fcm.membership)
  
  object[['Module.Pseudotime']] <- module.pseudotime
  return(object)
  
}


###############################################################################
#' Reproject samples into Fuzzy Clustering space, using the centers matrix
#' 
#' UMAP version
#' @param object An FCM object
#' @param objection.name A name for this projection
#' @param umap.configuration A umap::umap() configuration object
fcmUMAP <- function(object, 
                    proj.name = 'UMAP', 
                    umap.configuration) {
  
  centers <- getCenterMatrix(object)
  
  object[['UMAP']] <- umap::umap(
    t(centers),
    config = umap.configuration,
    method = 'naive'
  )
  
  umap.df <- data.frame(
    UMAP1 = object[['UMAP']][['layout']][ , 1],
    UMAP2 = object[['UMAP']][['layout']][ , 2]
  )
  
  if (is.null(object[['Projections']])) {
    object[['Projections']] <- list()
  }
  
  object[['Projections']][[proj.name]][['DF']] <- umap.df
  object[['Projections']][[proj.name]][['subspace']] <- centers
  
  return(object)
}


###############################################################################
#' Reproject samples into Fuzzy Clustering space, using the centers matrix
#' 
#' PCA version
#' @param object An FCM object
#' @param umap.configuration A umap::umap() configuration object
fcmPCA <- function(object) {
  
  centers <- getCenterMatrix(object)
  
  object[['PCA']] <- prcomp(t(centers))
  
  pca.df <- data.frame(
    PC1 = object[['PCA']][['x']][ , 1],
    PC2 = object[['PCA']][['x']][ , 2]
  )
  
  if (is.null(object[['Projections']])) {
    object[['Projections']] <- list()
  }
  
  object[['Projections']][[proj.name]][['DF']] <- pca.df
  object[['Projections']][[proj.name]][['subspace']] <- centers
  
  return(object)
}


###############################################################################
#' Perform clustering in FCM space on the pseudobulks or cells
#' 
#' Takes the object, resolutions, and automatically performs at several 
#' resolutions. Can also be used to find reclustering centroids, if a 
#' resolution has already been set. 
fcmProjectionReclustering <- function(object, 
                                      proj.name = 'UMAP',
                                      resolutions = seq(0.2,1.2, 0.2), 
                                      find.centroids = T) {
  
  projection <- object[['Projections']][[proj.name]]

  if (is.null(projection)) {
    stop('Projection not found.')
  }
  
  if (find.centroids == F) {
    
    subspace <- projection[['subspace']]
    rownames(subspace) <- paste0("module", 1:nrow(subspace))
    colnames(subspace) <- paste0("cell", 1:ncol(subspace))
    
    clusts <- do.call(
      "data.frame", 
      lapply(
        X = resolutions, 
        FUN = function(cr) {
          clust.assignments <- clusterSeurat(subspace, resolution = cr)
          return(clust.assignments)
        }
      )
    )
    colnames(clusts) <- paste0('cr_', resolutions)
    
    clusterings <- as.list(clusts)
    
    object[['Projections']][[proj.name]][['Reclustering']] <- clusterings
    
    print(apply(clusts, 2, table))
  }
  
  if (find.centroids == T) {
    
    clusterings <- projection[['Reclustering']]
    resolution <- clusterings[['Resolution']]
    
    if (is.null(resolution)) {
      stop("No resolution set. \n")
    }
    
    if (is.null(clusterings)) {
      stop("No reclustering performed. Run first with find.centroids = F")
    }
    
    clust.result <- clusterings[[ resolution ]]
    
    if (is.null(clust.result)) {
      stop("Resolution not found...")
    }
    
    dim1 <- projection[['DF']][ , 1]
    dim2 <- projection[['DF']][ , 2]
    
    cluster.centroids <- do.call(
      rbind,
      lapply(
        X = levels(clust.result), 
        FUN = function(x) {
          
          cells <- which(clust.result == x)
          x0 <- mean(dim1[cells])
          y0 <- mean(dim2[cells])
          
          return(c(x0,y0))
        }
      )
    )
    colnames(cluster.centroids) <- c("UMAP1", "UMAP2")
    cluster.centroids <- as.data.frame(cluster.centroids)
    nclust <- nrow(cluster.centroids)
    
    cluster.centroids$clusters <- paste0("c", seq(0, nclust - 1) )
    
    # Add to object.
    object[['Projections']][[proj.name]][['Reclustering']][['Centroids']] <- cluster.centroids
  }
  
  return(object)
}


###############################################################################
#' Make a tidy dataframe from UMAP and PCA projections. 
#' 
#' Takes the pseudobulk colData, the UMAP and PCA projections, and a Louvain
#' reclustering of the UMAP space, with resolution
addFCMProjectionDataFrame <- function(object, 
                                      metadata,
                                      proj.name = 'UMAP') {
  
  projection <- object[['Projections']][[proj.name]]
  reclustering <- projection[['Reclustering']]
  cluster.resolution <- reclustering[['Resolution']]
  clusters <- reclustering[[cluster.resolution]]
  
  if (is.null(clusters)) {
    stop('No reclustering added yet.')
  }
  
  proj.df <- projection[['DF']]
  out.df <- cbind(
    proj.df,
    metadata,
    clusters = clusters
  )
  
  out.df <- out.df[ , !duplicated(colnames(out.df))] 
  
  object[['Projections']][[proj.name]][['DF']] <- out.df
  return(object)
}


###############################################################################
#' Set Louvain reclustering resolution
#' 
#' @param object AN FCM object
#' @param clustering.name The name of the clustering
#' @param resolution A character vector or index specifying the resolution,
#' should correspond to the names or indices of object[['Reclustering']].
setReclusterResolution <- function(object, proj.name, resolution) {
  object[['Projections']][[proj.name]][['Reclustering']][['Resolution']] <- resolution
  l <- 
    length(
      unique(
        object[['Projections']][[proj.name]][['Reclustering']][[resolution]]
      )
    )
  object[['Projections']][[proj.name]][['Reclustering']][['nclust']] <- l
  return(object)
}


###############################################################################
#' Compute the expression-based centroids of each module
#' 
#' Finds the 'average cell' of a cluster and finds its projection into pre-
#' computed module UMAP space. Also connects each module by the extent of 
#' overlap with other modules, from fuzzy clustering
#' 
#' Think of modules as a PCA
#'
#' Centroids: [Modules X Cells] * T[Genes * Cells] * [GeneMembership X Modules]
#' result: A 14 x 14 matrix
computeModuleExpressionCentroids <- function(object, scale.basis = T, threshold = NULL) {
  
  # Get objects
  
  centers <- getCenterMatrix(object)
  membership <- getMembership(object)
  
  if (scale.bases == T) {
    mat <- getDataMatrix(object)
  } else {
    mat <- object[['Counts.Matrix']]
    if (is.null(mat)) {
      stop('No counts matrix located ( object[[\'Counts.Matrix\']]')
    }
  }
  
  if (!is.null(threshold)) {
    membership[membership < threshold] <- 0
  }
  
  # Compute projection
  
  module.space <- scale(centers %*% t(mat) %*% membership)
  centroid.proj <- predict(object[['UMAP']], module.space)
  colnames(centroid.proj) <- c('UMAP1', 'UMAP2')
  
  # Make data frame
  
  module.pseudotime <- object[['Module.Pseudotime']]
  
  centroid.df <- data.frame(
    centroid.proj,
    pseudotime = module.pseudotime,
    module.name = names(module.pseudotime)
  )
  rownames(centroid.df) <- centroid.df$module.name
  
  object[['Module.Centroids']] <- centroid.df
  return(object)
}


###############################################################################
#' Draw a connectivity graph.
#' 
#' using the specified overlap method, makes a dataframe connecting module 
#' centroids. 
#' @param object An FCM object
#' @param overlap.method Which overlap calculation to use. See computeOverlap.
drawModuleConnectivity <- function(object, overlap.method = 'ratio') {
  
  if (is.null(object[['Module.Centroids']])) {
    stop('No module centroids computed. Run computeModuleExpressionCentroids')
  }
  
  centroid.df <- object[['Module.Centroids']]
  
  fuzz.overlap <- object[['Overlaps']][[overlap.method]]
  
  overlap.df <- fuzz.overlap %>%
    melt(value.name = 'overlap', varnames = c('from', 'to'))
  
  overlap.df$x1 <- centroid.df[overlap.df$from, 'UMAP1']
  overlap.df$x2 <- centroid.df[overlap.df$to, 'UMAP1']
  overlap.df$y1 <- centroid.df[overlap.df$from, 'UMAP2']
  overlap.df$y2 <- centroid.df[overlap.df$to, 'UMAP2']
  
  object[['Overlap.Segments']][[overlap.method]] <- overlap.df
  return(object)
}


###############################################################################
#' Filter module connectivity by a threshold, and weight the connectivity 
#' between the modules. Prepares an object for easy plotting.
#' 
#' @param threshold Whether or not to override the object membership threshold
weightModuleConnectivity <- function(object, overlap.method, threshold = NULL) {
  
  if (is.null(threshold)) {
    threshold <- getThreshold(object)
  }
  
  if (is.null(object[['Overlap.Segments']][[overlap.method]])) {
    stop('No overlap segments found for method. Run drawModuleConnectivity')
  }
  
  overlap.df <- object[['Overlap.Segments']][[overlap.method]]
  
  filt.overlap.df <- overlap.df %>% filter(from != to & overlap > threshold)
  filt.overlap.df$alpha <- filt.overlap.df$overlap
  
  object[['Overlap.Segments.Filtered']][[overlap.method]] <- filt.overlap.df
  return(object)
}


###############################################################################
#' Plot connectivity between module expression centroids
#' 
#' @param object An FCM project
#' @param overlap.method The method used to compute overlaps
#' @param color.varaible The covariate used to color points underneath the 
#' module connectivity graph. 
plotModuleConnectivity <- function(object, 
                                   overlap.method = 'ratio', 
                                   proj.name = 'UMAP',
                                   color.variable = 'pseudotime') {
  
  if (is.null(object[['Projections']][[proj.name]][['DF']])) {
    stop('No tidy projection data found.')
  }
  
  if (is.null(object[['Module.Centroids']])) {
    stop('No module centroids found.')
  }
  
  if (is.null(object[['Overlap.Segments.Filtered']][[overlap.method]])) {
    stop('No overlap segments found.')
  }
  
  projections <- object[['Projections']][[proj.name]][['DF']]
  centroid.df <- object[['Module.Centroids']]
  overlap.df  <- object[['Overlap.Segments.Filtered']][[overlap.method]]
  trial.name  <- object[['Trial.Name']]
  
  cp <- ggplot(projections, aes_string(x = 'UMAP1', y = 'UMAP2', color = color.variable)) +
    geom_point(alpha = 0.2, size = 2, shape = 16) +
    guides(color = guide_legend(override.aes = list(alpha = 1)),
           alpha = guide_legend(title = str_to_title(overlap.method))) +
    geom_point(data = centroid.df, color = 'black', shape = 8, size = 1.5) +
    geom_label_repel(data = centroid.df, color = 'black', aes(label = paste0('m', module.name)), size = 3, label.padding = 0.1) +
    geom_segment(data = overlap.df, inherit.aes = F, aes(x = x1, y = y1, xend = x2, yend = y2, alpha = alpha)) +
    theme_classic() +
    umap.theme +
    theme(aspect.ratio = 1,
          legend.position = "right") +
    scale_color_viridis(option = "C") +
    labs(title = "Module centroids and connectivity")
  
  pdf.name <- sprintf(
    'fuzzy_clusters/%s/Projected_Module_Connectivity_%s.pdf', 
    trial.name, 
    overlap.method
  )
  pdf(
    file = pdf.name,
    width = 5,
    height = 5,
    useDingbats = F
  )
  print(cp)
  dev.off()
  return(cp)
}


###############################################################################
#' Plot covariates in the FCM projection
#' 
#' Covariates should come from the projection data frame. You can pass titles
#' and color palettes and alphas in, and you can specify *which* projection 
#' gets used. For instance, single cell, or pseudobulk. 
#' @param object an FCM object
#' @param covariate A metadata column
#' @param plot.title A character string
#' @param color.palette A vector of colors. If the covariate is identified as
#' a factor, it will use the vector as is. If the covariate is identified as
#' numeric / continuous, it will interpolate the colors.
#' @param alpha Plot point alpha
#' @param objection Name (in the FCM object) of the projection to use
#' @param add.layers A dataframe that can be concisely passed to 
#' additional.geoms. 
#' @param add.geoms Character strings that explicitly define additional
#' ggproto objects, which will be parsed and added to the output ggplot.
plotProjection <- function(object, 
                           proj.name = 'UMAP',
                           covariate = c('group', 'age', 'pseudotime', 'clusters'),
                           plot.title = c('Original cluster', 'Age', 'Pseudotime', 'Reclustering'),
                           color.palette,
                           alpha = 1,
                           point.size = 1,
                           add.layers = NULL,
                           add.geoms = NULL) {
  
  DF <- object[['Projections']][[proj.name]][['DF']]
  
  if (is.null(DF)) {
    stop('No tidy projection data frame found.')
  }
  
  if (is.null(DF[[covariate]])) {
    stop('Covariate not found in projection data frame.')
  }
  
  if (grepl('umap', tolower(proj.name))) {
    dim1 <- 'UMAP1'
    dim2 <- 'UMAP2'
  } else if (grepl('pc', tolower(proj.name))) {
    dim1 <- 'PC1'
    dim2 <- 'PC2'
  } else {
    dim1 <- 'UMAP1'
    dim2 <- 'UMAP2'
  }
  
  if (class(DF[[covariate]]) == 'factor' | class(DF[[covariate]]) == 'character') {
    color <- 'scale_color_manual(values = color.palette)'
    color.expr <- eval(parse(text = color))
  }
  
  if (class(DF[[covariate]]) == 'numeric') {
    color <- 'scale_color_gradientn(colors = color.palette)'
    color.expr <- eval(parse(text = color))
  }
  
  p <- ggplot(DF, aes_string(x = dim1, y = dim2, color = covariate)) +
    geom_point(alpha = alpha, size = point.size, shape = 16) +
    theme_classic() +
    umap.theme +
    theme(aspect.ratio = 1,
          legend.position = "right") +
    labs(title = plot.title) +
    color.expr
  
  if (!is.null(add.geoms)) {
    
    geom.expr <- lapply(add.geoms, function(x) {
      eval(parse(text = x))
    })
    
    p <- p + geom.expr
  }
  
  return(p)
}


###############################################################################
addLinkedATACModules <- function(object, 
                                 link.object, 
                                 bg.link.object = NULL,
                                 make.unique = F) {
  
  fcm.membership <- getMemberGeneList(object)
  
  fcm.atac <- lapply(
    fcm.membership, 
    function(x) {
      findLinkedPeaks(genes = x, object = link.object, make.unique = make.unique)
    }
  )
  
  if (!is.null(bg.link.object)) {
    fcm.background.atac <- lapply(
      fcm.membership, 
      function(x) {
        findLinkedPeaks(genes = x, object = bg.link.object, make.unique = make.unique)
      }
    )
    
    object[['Module.Background.Linked.Peaks']] <- fcm.background.atac
  }
  
  object[['Module.Linked.Peaks']] <- fcm.atac
  
  return(object)
}







###############################################################################
#' Impute missing features from an orthogonal dataset
#' 
#' Allows easy projection into the FCM space. Number of features must match.
#' Function also outputs the list of genes that needed imputing (that were 
#' missing from the new query dataset). Returns the actual matrix.
#' 
#' @param object An FCM object
#' @param newdata.matrix A new data matrix to project (missing values will be
#' imputed)
#' @param convert.ids Whether or not to convert IDs from symbols to ENSEMBL
#' @param impute.fun Not fully implemented -- a function to use for imputation
#' @param verbose Whether to print the genes that are being imputed
imputeMissingFeatures <- function(object = NULL, 
                                  mat = NULL,
                                  newdata.matrix,
                                  convert.ids = T,
                                  impute.fun = 'median',
                                  verbose = T) {
  
  if (is.null(object)) {
    mat <-  mat
  } else {
    mat <- getDataMatrix(object)
  }

  mat2 <- newdata.matrix
  
  if (convert.ids == T) {
    rownames(mat2) <- convertGeneIDs(rownames(mat2), name2ensembl)
  }
  
  missing.data <- rownames(mat)[! rownames(mat) %in% rownames(mat2)]
  
  fix.mat <- do.call(
    rbind,
    lapply(seq_along(rownames(mat)), function(x) {
      
      id <- rownames(mat)[x]
      
      if (! id %in% rownames(mat2)) {
        
        if (impute.fun == 'median') {
          out <- rep(median(mat2), ncol(mat2))          
        } else if (impute.fun == 'zero') {
          out <- rep(0, ncol(mat2))
        }

      } else {
        out <- mat2[id, ]
      }
        
      return(out)
    })
  )
  rownames(fix.mat) <- rownames(mat)
  
  if (convert.ids == T) {
    missing.data <- convertGeneIDs(missing.data, ensembl2name)
  }
  
  if (verbose == T) {
    cat(
      sprintf(
        'Imputed genes not found in query: %s\n\n', 
        as.character(length(missing.data))
      )
    )
    if (length(missing.data) == 0) {
      report <- 'None.'
    } else {
      report <- paste0(missing.data, collapse = ', ')
    }
    cat(report)
    cat('\n')
  }
  
  return(fix.mat)
}


###############################################################################
#' Project data into the FCM space. 
#' 
#' This could be single cell RNA seq or ATAC seq data, or something else that
#' can be coerced to match the features of the FCM object. This function adds
#' an element to the FCM object [['Projections']][[output.name]]
#' 
#' The goal is a module X cell matrix, where each row (module) is scaled across
#' all the cells. For single cells: 
#' [Modules X Genes] * [Genes X Cells] = [Modules X Cells] 
#' 
projectCells <- function(object, 
                         cell.matrix,
                         proj.name,
                         metadata,
                         add.matrix = T,
                         as.separate.object = F) {
  
  # get matrices
  gene.ids <- rownames(getDataMatrix(object))
  membership <- getMembership(object)
  
  # calculate membership matrix product
  filt.cell.matrix <- cell.matrix[gene.ids, ]
  product.mat <- rowScale(t(membership) %*% filt.cell.matrix)
  
  pred <- predict(object[['UMAP']], t(product.mat))
  
  # add matrix
  
  if (is.null(object[['Projections']][[proj.name]])) {
    object[['Projections']][[proj.name]] <- list()
  }
  
  if (add.matrix == T) {
    
    object[['Projections']][[proj.name]][['Data.Matrix']] <- filt.cell.matrix
    object[['Projections']][[proj.name]][['subspace']] <- product.mat
    
  }
  
  pred.df <- data.frame(
    UMAP1 = pred[ ,1], 
    UMAP2 = pred[ ,2]
  )
  
  proj.df <- cbind(
    pred.df,
    metadata
  )
  
  object[['Projections']][[proj.name]][['DF']] <- proj.df
  
  if (as.separate.object == F) {
    return(object)
  } else {
    return(list(
      Data.Matrix = filt.cell.matrix, 
      subspace = product.mat,
      DF = proj.df
    ))
  }
}

###############################################################################
#' Project data from a subset of genes into the FCM space. 
#' 
#' Other features will be zeroed out ~ note that this is a kind of trick
#' 
#' Projection mat could be single cell RNA seq or ATAC seq data, or something 
#' else that can be coerced to match the features of the FCM object. This 
#' function can add elements to the FCM object [['Projections']][[output.name]]
#' 
projectSubset <- function(genes,
                          object = fuzz, 
                          base.matrix,
                          proj.name,
                          metadata,
                          reference.metadata,
                          impute.function = 'zero',
                          plot.density = .5,
                          add.matrix = F,
                          as.separate.object = T) {

  cat('Imputing...\n')
  imputed.matrix <- imputeMissingFeatures(
    object = object, 
    newdata.matrix = base.matrix[genes, ], 
    convert.ids = T,
    impute.fun = impute.function,
    verbose = F
  )
  
  cat('Projecting...\n')
  projection <- projectCells(
    object = object, 
    cell.matrix = imputed.matrix, 
    proj.name = proj.name, 
    metadata = metadata,
    add.matrix = F,
    as.separate.object = F
  )
  
  out <- projection[['Projections']][[proj.name]][['DF']] 
  
  cat('Joining with reference...\n')
  plot.df <- out %>%
    inner_join(reference.metadata, by = c('pseudotime', 'group'))
  
  cat('Plotting...\n')
  
  # From here this code is hard-written for our particular application...
  
  pdf.name <- sprintf('fuzzy_clusters/%s/%s_projection.pdf', getTrialName(object), proj.name)
  plot.name1 <- sprintf('Gene subset projection, all samples: %s', proj.name)
  plot.name2 <- sprintf('Gene subset projection, branches: %s', proj.name)
  
  pdf(pdf.name, useDingbats = F, width = 5.5, height = 4.5)
  
  pr1 <- plotProjection(
    object = projection, 
    proj.name = proj.name, 
    covariate = 'group', 
    plot.title = plot.name1, 
    color.palette = glial.colors.atac$color
  )
  pr2 <- plotProjection(
    object = projection, 
    proj.name = proj.name, 
    covariate = 'is.Branch', 
    plot.title = plot.name2, 
    color.palette = c(glial.colors.atac$color[c(6,5,3)], 'grey67')
  )
  
  # Arrows
  pr3 <- ggplot(plot.df, aes(x = UMAP1.y, y = UMAP2.y)) + 
    geom_point(color = 'grey67') +
    geom_segment(data = plot.df %>% filter(is.Branch != 'Not Branch'),
                 arrow = arrow(length = unit(.1, 'inches')),
                 aes(x = UMAP1.y, xend = UMAP1.x, y = UMAP2.y, yend = UMAP2.x, color = is.Branch),
                 alpha = .5) +
    geom_point(data = plot.df %>% filter(is.Branch != 'Not Branch'), 
               aes(x = UMAP1.x, y = UMAP2.x, color = is.Branch)) +
    theme_classic() +
    theme(aspect.ratio = 1) +
    scale_color_manual(values = glial.colors.atac$color[c(6,5,3)]) +
    labs(title = plot.name1)
  
  sampled.plot.df <- plot.df %>% 
    sample_n(ceiling(nrow(plot.df) * plot.density), replace = F)
  
  pr4 <- ggplot(sampled.plot.df, aes(x = UMAP1.y, y = UMAP2.y)) + 
    geom_point(color = 'grey67') +
    geom_segment(data = sampled.plot.df,
                 arrow = arrow(length = unit(.1, 'inches')),
                 aes(x = UMAP1.y, xend = UMAP1.x, y = UMAP2.y, yend = UMAP2.x, color = group),
                 alpha = .5) +
    theme_classic() +
    theme(aspect.ratio = 1) +
    scale_color_manual(values = glial.colors.atac$color) +
    labs(title = plot.name2)
  
  print(pr1)
  print(pr2)
  print(pr3)
  print(pr4)
  
  dev.off()
  
  return(plot.df)
  
}

###############################################################################
#' Quick save function
#' 
#' @param object FCM object
#' @param file A character file name. Defaults to trial name
saveFCM <- function(object, file = NULL) {
  
  if (is.null(file)) {
    trial.name <- object[['Trial.Name']]
    file <- sprintf('fuzzy_clusters/%s/FCM_Object.RDS', trial.name)
  }
  saveRDS(object, file)
}


###############################################################################
#' Get genes in a particular module
#' 
#' Quickly access genes from a module. Takes integer input, where integers
#' correspond to module names as well (1 = m1)
getModules <- function(module, object = fuzz) {
  return(getMemberGeneList(object)[[module]])
}


###############################################################################
#' Which genes are transcription factors in a list?
#' 
#' @param genes A vector of genes
#' @param return What to return. Logical, or names
transcriptionFactors <- function(genes, return = 'names') {
  if (return == 'names') {
    return(genes[genes %in% tfs])
  } else {
    return(genes %in% tfs)
  }
}


###############################################################################
#' Plot gene expression in the FCM object.
#' 
#' Convenience function for quick access to data.
#' Parameters follow geneSetAveragePlot so they won't be listed here.
plotExpression <- function(genes,
                           project = fuzz,
                           plot.object = fuzz[['Data.Matrix']], 
                           plot.object.class = 'matrix', 
                           proj.name = 'UMAP', 
                           scale.expression = F,
                           plot.type = NULL,
                           convert = T,
                           point.size = 2,
                           aspectratio = 1,
                           trim = NULL,
                           color.palette = plant_palettes[['Iris']],
                           print = T,
                           facet.rows = NULL) {
  
  projection <- fuzz[['Projections']][[proj.name]][['DF']]
  
  if (is.null(plot.type)) {
    if (length(genes) == 1) {
      plot.type = 'single'
    } else {
      plot.type <- 'panels'
      if (is.null(facet.rows)) {
        facet.rows <- floor(length(genes)/3)
      }
    }
  }
  
  geneSetAveragePlot(
    genes = genes, 
    object = plot.object, 
    object.class = plot.object.class,
    projection = projection,
    scaled = scale.expression,
    trim = trim,
    plot.type = plot.type,
    convert = convert,
    point.size = point.size,
    aspectratio = aspectratio,
    color.palette = color.palette,
    print = print, 
    num.panel.rows = facet.rows
  )
  
  
}

###############################################################################
#' Plot gene score in the FCM object.
#' 
#' Convenience function for quick access to data.
#' Parameters follow geneSetAveragePlot so they won't be listed here.
plotActivity <- function(genes,
                         project = fuzz,
                         plot.object = fuzz[['Projections']][['ATAC.UMAP']][['Data.Matrix']], 
                         plot.object.class = 'matrix', 
                         proj.name = 'ATAC.UMAP', 
                         scale.expression = F,
                         plot.title = 'ATAC Gene Scores',
                         plot.type = NULL,
                         convert = T,
                         point.size = 2,
                         aspectratio = 1,
                         guide.name = 'Gene score',
                         color.palette = rev(brewer.pal(9, 'YlGnBu')),
                         trim = NULL) {
  
  projection <- fuzz[['Projections']][[proj.name]][['DF']]
  
  if (is.null(plot.type)) {
    if (length(genes) == 1) {
      plot.type = 'single'
      plot.title = genes
    } else {
      plot.type <- 'panels'
    }
  }
  
  geneSetAveragePlot(
    genes = genes, 
    object = plot.object, 
    object.class = plot.object.class,
    projection = projection,
    scaled = scale.expression,
    plot.title = plot.title,
    plot.type = plot.type,
    guide.name = guide.name,
    convert = convert,
    point.size = point.size,
    aspectratio = aspectratio,
    color.palette = color.palette,
    trim = trim
  )
}

###############################################################################
#' Plot gene score in the FCM object.
#' 
#' Convenience function for quick access to data.
#' Parameters follow geneSetAveragePlot so they won't be listed here.
plotCvar <- function(motifs,
                     project = fuzz,
                     plot.object = glial.pb.cvar.mat, 
                     plot.object.class = 'matrix', 
                     proj.name = 'ATAC.UMAP', 
                     scale.expression = T,
                     plot.title = 'ATAC ChromVAR deviations',
                     plot.type = NULL,
                     convert = F,
                     point.size = 2,
                     aspectratio = 1,
                     guide.name = 'Deviation',
                     color.palette = viridis_pal()(100),
                     trim = NULL) {
  
  projection <- fuzz[['Projections']][[proj.name]][['DF']]
  
  motifs <- toupper(motifs)
  motifs <- grep(motifs, toupper(rownames(plot.object)), value = T)
  
  if (is.null(plot.type)) {
    if (length(motifs) == 1) {
      plot.type = 'single'
      plot.title = motifs
    } else {
      plot.type <- 'panels'
    }
  }
  print(motifs)
  geneSetAveragePlot(
    genes = motifs, 
    object = plot.object, 
    object.class = plot.object.class,
    projection = projection,
    scaled = scale.expression,
    plot.title = plot.title,
    plot.type = plot.type,
    guide.name = guide.name,
    convert = convert,
    point.size = point.size,
    aspectratio = aspectratio,
    color.palette = color.palette,
    trim = trim
  )
  
  
}


###############################################################################
addModuleIgraph <- function(object, 
                            overlap.method = 'ratio', 
                            threshold = NULL, 
                            graph.mode = 'undirected') {
  
  if (is.null(threshold)) {
    threshold <- getThreshold(object)
  }
  
  adjacency.matrix <- getOverlaps(object, overlap.method)
  adjacency.matrix[adjacency.matrix < threshold] <- 0
  
  graph <- igraph::graph_from_adjacency_matrix(
    adjmatrix = adjacency.matrix, 
    mode = graph.mode, 
    weighted = T, 
    diag = F
  )
  
  set.vertex.attribute(graph, 'pseudotime', value = fuzz[['Module.Pseudotime']])
  
  object[['Graph']][[overlap.method]] <- graph
  plot.igraph(graph)
  return(object)
}





###############################################################################
#' Export lists of genes belonging to each module as TSV
#' 
#' @param object FCM object
#' @param table.name The name of the table. Defaults to Module_Gene_Lists
exportModuleGeneLists <- function(object, 
                                  core.only = F,
                                  table.name = NULL,
                                  return = F) {
  
  trial.name <- getTrialName(object)
  membership <- getMemberGeneList(object)
  
  df <- do.call(
    rbind,
    lapply(seq_along(membership), function(x) {
      
      genes <- membership[[x]]

      if (core.only == T) {
        other.mod.i <- seq_along(membership)[seq_along(membership) != x]
        other.mods <- Reduce('c', membership[other.mod.i])
        genes <- genes[!genes %in% other.mods]
      }
      
      ids <- names(genes)
      
      if (length(genes) == 0) {
        genes <- NA
        ids   <- NA
      }
      
      data.frame(
        cbind(
          module = paste0('m', as.character(x)), 
          gene.symbol = genes, 
          gene.id = ids
        )
      )
    })
  )
  
  if (is.null(table.name)) {
    if (core.only == T) {
      core <- "Core_"
    } else {
      core <- ""
    }
    table.name <- sprintf('fuzzy_clusters/%s/%sModule_Gene_Lists.txt', trial.name, core)
  }
  
  if (return == T) {
    return(df)
  } else {
    
    write.table(
      x = df, 
      file = table.name,
      quote = F,
      sep = "\t",
      row.names = F, 
      col.names = T
    )
    
  }
  
}

exportModuleOverlaps <- function(object, return = F) {
  
  trial.name = getTrialName(object)
  
  membership <- getMemberGeneList(object)
  
  df <- 
    lapply(seq_along(membership), function(i) {
      lapply(seq_along(membership), function(j) {
        
        intersections <- intersect(getModules(i), getModules(j))
        intersect.ids <- convertGeneIDs(intersections, name2ensembl)
      
        data.frame(
          cbind(
            module1 = paste0('m', as.character(i)), 
            module2 = paste0('m', as.character(j)),
            gene.symbol = intersections,
            gene.id = intersect.ids
          )
        )
      }) %>% Reduce('rbind', .)
    }) %>% Reduce('rbind', .) 
  
  df %<>% filter(module1 != module2)
  
  table.name <- sprintf('fuzzy_clusters/%s/Module_Gene_Overlaps.txt', trial.name)
  if (return == T) {
    return(df)
  } else {
    
    write.table(
      x = df, 
      file = table.name,
      quote = F,
      sep = "\t",
      row.names = F, 
      col.names = T
    )
    
  }
}

maxMemberhip <- function(object, threshold = NULL) {
  
  if (is.null(threshold)) {
    threshold <- getThreshold(object)
  }
  
  mm <- getMembership(object)
  
  apply(mm, 2, max)
  
}

overlap.plot.2 <- function(object,
                           overlap.method) {
  
  centers <- getCenterMatrix(object)
  fcm <- getMemberGeneList(object)
  fuzz.overlap <- object[['Overlaps']][[overlap.method]]
  threshold <- getThreshold(object)
  
  pca <- prcomp(centers, scale = TRUE)
  
  pcs <- t(t(pca[[2]]) %*% t(centers))
  pc.df <- data.frame(PC1 = pcs[ , 1], PC2 = pcs[ , 2])
  
  pc.df$n.genes <- do.call(c, lapply(fcm, length))
  pc.df$module <- paste0('m', seq_len(nrow(pc.df)))
  
  overlap.df <- fuzz.overlap %>%
    melt(value.name = 'overlap', varnames = c('from' ,'to')) %>%
    filter(from != to, overlap > threshold)
  
  overlap.df$x1 <- pc.df[overlap.df$from, 'PC1']
  overlap.df$x2 <- pc.df[overlap.df$to, 'PC1']
  overlap.df$y1 <- pc.df[overlap.df$from, 'PC2']
  overlap.df$y2 <- pc.df[overlap.df$to, 'PC2']
  
  gg <- ggplot(pc.df, aes(x = PC1, y = PC2, size = n.genes)) +
    geom_segment(data = overlap.df, inherit.aes = F, aes(x = x1, y = y1, xend = x2, yend = y2, alpha = overlap)) +
    geom_point(color = 'black', fill = '#FF3300', shape = 21) +
    theme_classic() +
    theme(aspect.ratio = 1) +
    geom_text_repel(aes(label = module), size = 3, point.padding = .3, min.segment.length = 1.5) +
    labs(title = 'Relationships between modules')
  print(gg)
}

plotModuleHeatmap <- function(object) {
  
  dat <- getDataMatrix(object)
  fcm <- getMemberGeneList(object)
  col.data <- object[['colData']]
  
  mod.order <- order(object[['Module.Pseudotime']])
  row.order <- names(fcm[mod.order])
  col.order <- order(col.data[['pseudotime']], decreasing = F)
  
  dat <- dat[] 
  
}

isInModule <- function(gene, object = fuzz) {
  
  fcm <- getMemberGeneList(object)
  
  is.in <- do.call(c, lapply(seq_along(fcm), function(x) { gene %in% fcm[[x]] } ))
  
  return(names(fcm)[is.in])
  
}


###############################################################################
#' Plot membership extent in modules
#' 
#' I use this to get color scales which can be arranged as a graph by module 
#' connectivity
#' @param genes
#' @param object
#' @param ORD the order of the factors/modules in the plot
plotGeneMembership <- function(gene,
                               object = fuzz,
                               membership.quantiles = c(0, 0.06, 0.12, 0.24, 0.48, 1),
                               quantile.color.vals = colorRampPalette(plant_palettes$`Jasminium polyanthum A`[1:8])(5)) {
  
  df <- cbind(fuzz$Module.Centroids, geneMembershipDF(genes = gene))
  
  df %<>%
    mutate(Membership.Score = cut(
      x = Membership, 
      breaks = membership.quantiles, labels = c(1,2,3,4,5), 
      ordered_result = T
    ))
  
  gg <- plotExpression(
    genes = gene, 
    color.palette = brewer.pal(9, "Greys")[2:8], 
    plot.object = fuzz$Counts.Matrix, 
    print = F,
    point.size = 2
  ) +
    theme(aspect.ratio = 1,
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank()) +
    geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2),
                 data = segments, 
                 inherit.aes = F, 
                 color = 'black') +
    geom_point(inherit.aes = F, size = 8, data = df, shape = 21, aes(x = UMAP1, y = UMAP2, fill = Membership.Score)) +
    geom_text(inherit.aes = F, size = 3, data = df, aes(x = UMAP1, y = UMAP2, label = Module)) +
    scale_fill_manual(values = quantile.color.vals, drop = F) +
    labs(title = gene) 
  print(gg)
  
}


###############################################################################
#' Plot motif enrichments between modules
#' 
#' I use this to get color scales which can be arranged as a graph by module 
#' connectivity
plotMotifsInModules <- function(motif, 
                                enr = mod.motif.enr, 
                                to.plot = 'pval', 
                                color.pal = brewer.pal(9, 'Blues'), 
                                ORD = c(6,9,11,5,14,2,10,3,8,12,1,7,4,13)) {
  
  titles <- paste0(toupper(motif), collapse = ' ')
  motif %<>% toupper(.)
  
  if (grepl('pval|p.val', to.plot)) {
    to.plot <- 'neg.log10.p'
  } else {
    to.plot <- 'log2FoldEnrichment'
  }
  
  df <- mod.motif.enr %>%
    filter(short.name == motif) %>%
    mutate(neg.log10.p = ifelse(log2FoldEnrichment < 0, 0, neg.log10.p))
  
  if (! is.null(ORD)) {
    df$group <- factor(df$group, levels = paste0('m', as.character(ORD)))
  }
  
  gg <- ggplot(df, aes_string(x = 'group', y = '1', fill = to.plot)) +
    geom_point(size = 9, shape = 21) +
    scale_fill_gradientn(colors = color.pal)
  
  # pdf(sprintf("fuzzy_clusters/redo_0813_14cl_1.25_noCC/Motif_Module_Plots/%s.pdf", titles),
  #     height = 1,
  #     width = 4,
  #     useDingbats = F)
  print(gg)
  # dev.off()
}


###############################################################################
#' Return gene membership in modules as DF
#' 
#' I use this to get color scales which can be arranged as a graph by module 
#' connectivity
#' @param genes
#' @param object
geneMembershipDF <- function(genes, object = fuzz) {
  
  titles <- paste0(toupper(genes), collapse = ' ')
  membership <- getMembership(object)
  gene.ids <- convertGeneIDs(toupper(genes), name2ensembl)
  
  if (is.null(membership)) {
    stop('No membership matrix detected.')
  }
  if (! all(gene.ids %in% rownames(membership))) {
    stop('Genes not found in membership matrix')
  }
  mmat <- matrix(membership[gene.ids, ], ncol = ncol(membership))
  m <- apply(mmat, 2, mean)
  
  df <- data.frame(
    Module = paste0('m', colnames(membership)), 
    Membership = m
  )
  
  return(df)
}


###############################################################################
#' Return motif enrichment in modules as DF
#' 
#' I use this to get color scales which can be arranged as a graph by module 
#' connectivity
#' @param motif
#' @param object
motifModuleDF <- function(motif, object = fuzz, to.plot = 'pval') {
  titles <- paste0(toupper(motif), collapse = ' ')
  motif %<>% toupper(.)
  
  df <- mod.motif.enr %>%
    filter(motif.name == motif) %>%
    mutate(neg.log10.p = ifelse(log2FoldEnrichment < 0, 0, 2 * neg.log10.p))
  
  return(df) 
}






















