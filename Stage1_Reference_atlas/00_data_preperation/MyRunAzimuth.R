MyRunAzimuth = function (query, reference, query.modality = "RNA", annotation.levels = NULL, 
          umap.name = "ref.umap", do.adt = FALSE, verbose = TRUE, 
          assay = NULL, k.weight = 50, n.trees = 20, mapping.score.k = 100, 
          ...) 
{
  CheckDots(...)
  assay <- assay %||% DefaultAssay(query)
  if (query.modality == "ATAC") {
    query <- RunAzimuthATAC(query = query, reference = reference, 
                            annotation.levels = annotation.levels, umap.name = umap.name, 
                            verbose = verbose, assay = assay, k.weight = k.weight, 
                            n.trees = n.trees, mapping.score.k = mapping.score.k, 
                            ...)
  }
  else {
    if (dir.exists(reference)) {
      reference <- LoadReference(reference)$map
    }
    else {
      reference <- tolower(reference)
      if (reference %in% InstalledData()$Dataset) {
        reference <- LoadData(reference, type = "azimuth")$map
      }
      else if (reference %in% AvailableData()$Dataset) {
        InstallData(reference)
        reference <- LoadData(reference, type = "azimuth")$map
      }
      else {
        possible.references <- AvailableData()$Dataset[grepl("*ref", 
                                                             AvailableData()$Dataset)]
        print("Choose one of:")
        print(possible.references)
        stop(paste("Could not find a reference for", 
                   reference))
      }
      if (!"num_precomputed_nns" %in% names(Misc(reference[["refUMAP"]])$model)) {
        Misc(reference[["refUMAP"]], slot = "model")$num_precomputed_nns <- 1
      }
      key.pattern = "^[^_]*_"
      new.colnames <- gsub(pattern = key.pattern, replacement = Key(reference[["refDR"]]), 
                           x = colnames(Loadings(object = reference[["refDR"]], 
                                                 projected = FALSE)))
      colnames(Loadings(object = reference[["refDR"]], 
                        projected = FALSE)) <- new.colnames
    }
    dims <- as.double(slot(reference, "neighbors")$refdr.annoy.neighbors@alg.info$ndim)
    if (isTRUE(do.adt) && !("ADT" %in% Assays(reference))) {
      warning("Cannot impute an ADT assay because the reference does not have antibody data")
      do.adt = FALSE
    }
    reference.version <- ReferenceVersion(reference)
    azimuth.version <- as.character(packageVersion(pkg = "Azimuth"))
    seurat.version <- as.character(packageVersion(pkg = "Seurat"))
    meta.data <- names(slot(reference, "meta.data"))
    if (is.null(annotation.levels)) {
      annotation.levels <- names(slot(object = reference, 
                                      name = "meta.data"))
      annotation.levels <- annotation.levels[!grepl(pattern = "^nCount", 
                                                    x = annotation.levels)]
      annotation.levels <- annotation.levels[!grepl(pattern = "^nFeature", 
                                                    x = annotation.levels)]
      annotation.levels <- annotation.levels[!grepl(pattern = "^ori", 
                                                    x = annotation.levels)]
    }
    query <- ConvertGeneNames(object = query, reference.names = rownames(x = reference), 
                              homolog.table = "./00_code_ref_atlas/0_ref_data_preperation/homologs.rds")
    if (!all(c("nCount_RNA", "nFeature_RNA") %in% c(colnames(x = query[[]])))) {
      calcn <- as.data.frame(x = Seurat:::CalcN(object = query[[assay]]))
      colnames(x = calcn) <- paste(colnames(x = calcn), 
                                   assay, sep = "_")
      query <- AddMetaData(object = query, metadata = calcn)
      rm(calcn)
    }
    if (any(grepl(pattern = "^MT-", x = rownames(x = query)))) {
      query <- PercentageFeatureSet(object = query, pattern = "^MT-", 
                                    col.name = "percent.mt", assay = assay)
    }
    anchors <- FindTransferAnchors(reference = reference, 
                                   query = query, k.filter = NA, reference.neighbors = "refdr.annoy.neighbors", 
                                   reference.assay = "refAssay", query.assay = assay, 
                                   reference.reduction = "refDR", normalization.method = "SCT", 
                                   features = rownames(Loadings(reference[["refDR"]])), 
                                   dims = 1:dims, n.trees = n.trees, mapping.score.k = mapping.score.k, 
                                   verbose = verbose)
    refdata <- lapply(X = annotation.levels, function(x) {
      reference[[x, drop = TRUE]]
    })
    names(x = refdata) <- annotation.levels
    if (isTRUE(do.adt)) {
      refdata[["impADT"]] <- GetAssayData(object = reference[["ADT"]], 
                                          layer = "data")
    }
    query <- TransferData(reference = reference, query = query, 
                          query.assay = assay, dims = 1:dims, anchorset = anchors, 
                          refdata = refdata, n.trees = 20, store.weights = TRUE, 
                          k.weight = k.weight, verbose = verbose)
    query <- IntegrateEmbeddings(anchorset = anchors, reference = reference, 
                                 query = query, query.assay = assay, reductions = "pcaproject", 
                                 reuse.weights.matrix = TRUE, verbose = verbose)
    query[["query_ref.nn"]] <- FindNeighbors(object = Embeddings(reference[["refDR"]]), 
                                             query = Embeddings(query[["integrated_dr"]]), return.neighbor = TRUE, 
                                             l2.norm = TRUE, verbose = verbose)
    query <- Azimuth:::NNTransform(object = query, meta.data = reference[[]])
    query[[umap.name]] <- RunUMAP(object = query[["query_ref.nn"]], 
                                  reduction.model = reference[["refUMAP"]], reduction.key = "UMAP_", 
                                  verbose = verbose)
    query <- AddMetaData(object = query, metadata = MappingScore(anchors = anchors, 
                                                                 ndim = dims), col.name = "mapping.score")
  }
  return(query)
}