### Functions are from decontX package, modified in a way to generate not plots right away but data frames, to the merge differnet 
# decontamination assays and to plot the levels of marker expression percentage 
#https://rdrr.io/bioc/celda/src/R/plot_decontx.R

plotDecontXMarkerPercentage_df <- function(x, markers, groupClusters = NULL,
                                        assayName = c(
                                          "counts",
                                          "decontXcounts"
                                        ),
                                        z = NULL, threshold = 1,
                                        exactMatch = TRUE, by = "rownames",
                                        ncol = round(sqrt(length(markers))),
                                        labelBars = TRUE, labelSize = 3) {
  cellTypeLabels <- percent <- NULL # fix check note
  
  legend <- "none"
  # Check that list arguments are named
  if (!is(markers, "list") || is.null(names(markers))) {
    stop("'markers' needs to be a named list.")
  }
  
  temp <- .processPlotDecontXMarkerInupt(
    x = x,
    z = z,
    markers = markers,
    groupClusters = groupClusters,
    by = by,
    exactMatch = exactMatch
  )
  x <- temp$x
  z <- temp$z
  geneMarkerIndex <- temp$geneMarkerIndex
  geneMarkerCellTypeIndex <- temp$geneMarkerCellTypeIndex
  groupClusters <- temp$groupClusters
  xlab <- temp$xlab
  
  if (inherits(x, "SingleCellExperiment")) {
    # If 'x' is SingleCellExperiment, then get percentage
    # for each matrix in 'assayName'
    df.list <- list()
    for (i in seq_along(assayName)) {
      counts <- SummarizedExperiment::assay(
        x[geneMarkerIndex, ],
        assayName[i]
      )
      df <- .calculateDecontXBarplotPercent(
        counts,
        z,
        geneMarkerCellTypeIndex,
        threshold
      )
      df.list[[i]] <- cbind(df, assay = assayName[i])
    }
    df <- do.call(rbind, df.list)
    assay <- as.factor(df$assay)
    if (length(assayName) > 1) {
      legend <- "right"
    }
  } else {
    ## If 'x' is matrix, then calculate percentages directly
    counts <- x[geneMarkerIndex, ]
    df <- .calculateDecontXBarplotPercent(
      counts,
      z,
      geneMarkerCellTypeIndex,
      threshold
    )
    assay <- "red3"
    legend <- "none"
  }
  
  # Build data.frame for ggplots
  df <- cbind(df, cellTypeLabels = names(groupClusters)[df$cellType])
  df$cellTypeLabels <- factor(df$cellTypeLabels,
                              levels = names(groupClusters)
  )
  df <- cbind(df, markerLabels = names(markers)[df$markers])
  df$markerLabels <- factor(df$markerLabels, levels = names(markers))
  return(df)
}


.processPlotDecontXMarkerInupt <- function(x, z, markers, groupClusters,
                                           by, exactMatch) {
  
  # Process z and convert to a factor
  if (is.null(z) & inherits(x, "SingleCellExperiment")) {
    cn <- colnames(SummarizedExperiment::colData(x))
    if (!("decontX_clusters" %in% cn)) {
      stop("'decontX_clusters' not found in 'colData(x)'. Make sure you have
           run 'decontX' or supply 'z' directly.")
    }
    z <- SummarizedExperiment::colData(x)$decontX_clusters
  } else if (length(z) == 1 & inherits(x, "SingleCellExperiment")) {
    if (!(z %in% colnames(SummarizedExperiment::colData(x)))) {
      stop("'", z, "' not found in 'colData(x)'.")
    }
    z <- SummarizedExperiment::colData(x)[, z]
  } else if (length(z) != ncol(x)) {
    stop("If 'x' is a SingleCellExperiment, then 'z' needs to be",
         " a single character or integer specifying the column in",
         " 'colData(x)'. Alternatively to specify the cell cluster",
         " labels directly as a vector, the length of 'z' needs to",
         " be the same as the number of columns in 'x'. This is",
         " required if 'x' is a matrix.")
  }
  if (!is.factor(z)) {
    z <- as.factor(z)
  }
  
  if (!is.null(groupClusters)) {
    if (!is(groupClusters, "list") || is.null(names(groupClusters))) {
      stop("'groupClusters' needs to be a named list.")
    }
    
    # Check that groupClusters are found in 'z'
    cellMappings <- unlist(groupClusters)
    if (any(!(cellMappings %in% z))) {
      missing <- cellMappings[!(cellMappings %in% z)]
      stop(
        "'groupClusters' not found in 'z': ",
        paste(missing, collapse = ",")
      )
    }
    
    labels <- rep(NA, ncol(x))
    for (i in seq_along(groupClusters)) {
      labels[z %in% groupClusters[i]] <- names(groupClusters)[i]
    }
    na.ix <- is.na(labels)
    labels <- labels[!na.ix]
    x <- x[, !na.ix]
    z <- as.integer(factor(labels, levels = names(groupClusters)))
    xlab <- "Cell types"
  } else {
    labels <- as.factor(z)
    groupClusters <- levels(labels)
    names(groupClusters) <- levels(labels)
    xlab <- "Clusters"
  }
  
  # Find index of each feature in 'x'
  geneMarkerCellTypeIndex <- rep(
    seq(length(markers)),
    lapply(markers, length)
  )
  geneMarkerIndex <- retrieveFeatureIndex(unlist(markers),
                                          x,
                                          by = by,
                                          removeNA = FALSE,
                                          exactMatch = exactMatch
  )
  
  # Remove genes that did not match
  na.ix <- is.na(geneMarkerIndex)
  geneMarkerCellTypeIndex <- geneMarkerCellTypeIndex[!na.ix]
  geneMarkerIndex <- geneMarkerIndex[!na.ix]
  
  return(list(
    x = x,
    z = z,
    geneMarkerIndex = geneMarkerIndex,
    geneMarkerCellTypeIndex = geneMarkerCellTypeIndex,
    groupClusters = groupClusters,
    xlab = xlab
  ))
}


.calculateDecontXBarplotPercent <- function(counts,
                                            z,
                                            geneMarkerCellTypeIndex,
                                            threshold) {
  
  # Get counts matrix and convert to DelayedMatrix
  counts <- DelayedArray::DelayedArray(counts)
  
  # Convert to boolean matrix and sum markers in same cell type
  # The "+ 0" is to convert boolean to numeric
  counts <- counts >= threshold
  countsByMarker <- DelayedArray::rowsum(counts + 0, geneMarkerCellTypeIndex)
  countsByCellType <- DelayedArray::colsum((countsByMarker > 0) + 0, z)
  
  # Calculate percentages within each cell cluster
  zTotals <- tabulate(z)
  percentByCellType <- round(sweep(countsByCellType, 2, zTotals, "/") * 100)
  df <- reshape2::melt(percentByCellType,
                       varnames = c("markers", "cellType"),
                       value.name = "percent"
  )
  
  return(df)
}

