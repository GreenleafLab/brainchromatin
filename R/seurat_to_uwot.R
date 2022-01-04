
# assuming RunUMAP was with 'dims!=NULL', i.e. the umap was based on another reduction
getUwotModelFromSeurat <- function(seu, assayUsed=DefaultAssay(seu), reductionUsed="pca", ...){
	require(uwot)
	paramL <- seu@commands[[paste0("RunUMAP.",assayUsed,".",reductionUsed)]]@params
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

	X <- Embeddings(seu[[reductionUsed]])[, paramL[["dims"]]]

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

projectMatrix_SeuratUMAP <- function(X, seu, assayUsed=DefaultAssay(seu), dataType="counts") {
	# X <- GetAssayData(object=seu[["RNA"]], slot="counts") # for testing
	# colnames(X) <- paste0("cell", 1:ncol(X))
	if (!is.element(dataType, c("counts", "normalized", "scaled", "PCA"))){
		logger.error(c("Invalid dataType:", dataType))
	}

	paramL_norm <- seu@commands[[paste0("NormalizeData.",assayUsed)]]@params
	paramL_scale <- seu@commands[[paste0("ScaleData.",assayUsed)]]@params
	if (!paramL_scale[["do.scale"]]) logger.error("Don't know how to deal with unscaled data yet")
	paramL_pca <- seu@commands[[paste0("RunPCA.",assayUsed)]]@params
	paramL_umap <- seu@commands[[paste0("RunUMAP.",assayUsed,".",paramL_pca[["reduction.name"]])]]@params

	seuQ <- X # query seurat dataset
	if (class(X)!="Seurat"){
		seuQ <- CreateSeuratObject(counts=X, project="SEUQ", assay=assayUsed)
	}
	if (is.element(dataType, c("counts"))){
		logger.status("Normalizing ...")
		seuQ <- NormalizeData(
			seuQ,
			assay=paramL_norm$assay,
			normalization.method=paramL_norm$normalization.method,
			scale.factor=paramL_norm$scale.factor,
			verbose=FALSE
		)
	}
	if (!is.element(dataType, c("scaled", "PCA"))){
		logger.status("Scaling ...")
		seuQ <- ScaleData(
			seuQ,
			assay=paramL_scale$assay,
			features=paramL_scale$features,
			do.scale=paramL_scale$do.scale,
			do.center=paramL_scale$do.center,
			model.use=paramL_scale$model.use,
			verbose=FALSE
		)
	}

	if (!is.element(dataType, c("PCA"))){
		logger.status("PC projection ...")
		X <- GetAssayData(seuQ[[assayUsed]], slot="scale.data")
		if (!all(paramL_pca[["features"]] %in% rownames(X))){
			logger.error("Could not find all features in X for the PC projection")
		}
		X <- t(X[paramL_pca[["features"]],])
	}

	projM <- Loadings(seu, reduction=paramL_pca[["reduction.name"]])
	pcaCoord_proj <- X %*% projM
	# pcaCoord_orig <- Embeddings(seu[[paramL_pca[["reduction.name"]]]]) # compare the original

	logger.status("Retrieving UMAP model ...")
	umapRes <- getUwotModelFromSeurat(seu, assayUsed=assayUsed, reductionUsed=paramL_pca[["reduction.name"]])

	logger.status("UMAP projection ...")
	umapCoord_orig <- Embeddings(seu[[paramL_umap[["reduction.name"]]]])
	umapCoord_proj <- uwot::umap_transform(pcaCoord_proj[,paramL_umap[["dims"]]], umapRes)
	rownames(umapCoord_proj) <- rownames(pcaCoord_proj)
	colnames(umapCoord_proj) <- colnames(umapCoord_orig)
	res <- list(
		pcaCoord=pcaCoord_proj,
		pcsUsed=paramL_umap[["dims"]],
		umapCoord=umapCoord_proj
	)
	class(res) <- "SeuUMAPProjection"
	return(res)
}
