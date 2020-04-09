
#' Create Dendrogram
#'
#' Creates dendrogram object from bipartite hierarchical clustering for plotting
#'
#' @param labels vector of group names to label dendrogram
#' @param mergeMat matrix of dimension [n - 1, 2] that indicates merge events, where n is number of starting groups (see documentation for "merge" matrix contained in hclust class for details)
#' @param nodeHeights vector of dendrogram node heights, in order of merge
#' @return object of class dendrogram
#' @export
createHclust <- function(labels, mergeMat, nodeHeights) {
  # Creates dendrogram object from bipartite hierarchical clustering for plotting
  # Input:
  #   labels: vector of group names to label dendrogram
  #   mergeMat: matrix of dimension [n - 1, 2] that indicates merge events, where n is number of
  #             starting groups (see documentation for "merge" matrix contained in hclust class for details)
  #   nodeHeights: vector of dendrogram node heights, in order of merge
  # Output:
  #   hclustObj: object of class dendrogram
  
  # generate vector of order by which leaves are merged given matrix of merge events
  dendro_order <- dendroOrder(mergeMat)

  hclustLst <- list("merge" = mergeMat, "height" = nodeHeights, "order" = dendro_order, "labels" = labels)
  class(hclustLst) <- "hclust"
  hclustObj <- as.dendrogram(hclustLst)
  return(hclustObj)
}

#' Merge Order
#'
#' Creates vector of order by which leaves are merged given matrix of merge events
#'
#' @param mergeMat matrix of dimension [n - 1, 2] that indicates merge events, where n is number of starting groups (see documentation for "merge" matrix contained in hclust class for details)
#' @return vector of leaf (group) indices, in order of merge
#' @export
dendroOrder <- function(mergeMat) {
  # creates vector of order by which leaves are merged given matrix of merge events
  # Input:
  #   mergeMat: matrix of dimension [n - 1, 2] that indicates merge events, where n is number of
  #             starting groups (see documentation for "merge" matrix contained in hclust class for details)
  # Output:
  #   dendro_order: vector of leaf (group) indices, in order of merge
  mergeMatT <- t(mergeMat)
  dendro_order <- c()
  for (i in seq(length(mergeMat))) {
    if (mergeMatT[i] < 0) {
      dendro_order <- c(dendro_order, -mergeMatT[i])
    }
  }
  return(dendro_order)
}

#' Merge groups
#'
#' Chooses which entry (row, column) of dissimilarity matrix has minimum dissimilarity
#'
#' @param dissMat square matrix of nucleus-norm based dissimilarity measure
#' @return list containing column and row number of entry with minimum dissimilarity
#' @export
whichMerge <- function(dissMat) {
  # Chooses which entry (row, column) of dissimilarity matrix has minimum dissimilarity
  # Input:
  #   dissMat: square matrix of nucleus-norm based dissimilarity measure
  # Output:
  #   list containing column and row number of entry with minimum dissimilarity
  min_idx <- which.min(dissMat)
  col <- ceiling(min_idx / sqrt(length(dissMat)))
  row <- ((min_idx - 1) %% sqrt(length(dissMat))) + 1
  return(list("col" = col, "row" = row))
}

#' Bipartite Hierarchical Clustering
#'
#' Main bipartite hierarchial clustering algorithm. Execute browseVignettes("hierBipartite") for vignette on using the
#' hierBipartite package.
#'
#' @param X an n x p matrix (e.g. for gene expression)
#' @param Y an n x q matrix (e.g. for drug sensitivity)
#' @param groups a list of starting group membership (e.g. list("1" = c(1,2,3), "2" = c(4,5,6)) means group 1 has samples 1, 2, 3, and group 2 has samples 4, 5, 6.
#' @param n_subsample number of subsampling to generate matrix B (see paper)
#' @param subsampling_ratio fraction of samples to sample for subsampling to generate matrix B (see paper)
#' @param p.value boolean for whether to generate p-values for each merge
#' @param n_perm number of permutations for generating p-values. Ignored if p.value = FALSE
#' @param parallel boolean for whether to parallelize subsampling and p-value generation step
#' @return list of results from bipartite hierarchical clustering, containing 
#'         (1) nodeMemberships: list of samples (in terms of indices of X, Y) of each new merged group, in order of merge,
#'         (2) hclustObj: dendrogram class of resulting dendrogram, 
#'         (3) nodeSCCA: list of SCCA output for each new merged group, in order of merge, 
#'         (4) nodeGroups: list of groups for each merge, in order fo merge,
#'         (5) nodePvals: list of p-value of each new merge, in order of merge if p.value = TRUE.
#' @export
hierBipartite <- function(X, Y, groups, n_subsample = 1, subsampling_ratio = 1, p.value = FALSE, n_perm = 100, parallel = TRUE) {
  # Main bipartite hierarchial clustering algorithm
  # Input:
  #   X: an n x p matrix (e.g. for gene expression)
  #   Y: an n x q matrix (e.g. for drug sensitivity)
  #   groups: a list of starting group membership (e.g. list("1" = c(1,2,3), "2" = c(4,5,6)) means group 1 has 
  #           samples 1, 2, 3, and group 2 has samples 4, 5, 6.
  #   n_subsample: number of subsampling to generate matrix B (see paper)
  #   subsampling_ratio: fraction of samples to sample for subsampling to generate matrix B (see paper)
  #   p.value: boolean for whether to generate p-values for each merge
  #   n_perm: number of permutations for generating p-values. Ignored if p.value = FALSE
  #   parallel: boolean for whether to parallelize subsampling and p-value generation step
  # Output:
  #   retLst: list of results from bipartite hierarchical clustering, containing
  #     nodeMembership: list of samples (in terms of indices of X, Y) of each new merged group, in order of merge
  #     hclustObj: dendrogram class of resulting dendrogram
  #     nodeSCCA: list of SCCA output for each new merged group, in order of merge
  #     mergeGroups: list of groups for each merge, in order fo merge
  #     nodePvals: list p-value of each new merge, in order of merge if p.value = TRUE
  groupNames <- names(groups)
  
  # construct matrix B representing bipartite relationship
  bgraphs <- lapply(groups, function(gid) {
    mat1 <- X[gid, ]
    mat2 <- Y[gid, ]
    
    constructBipartiteGraph(mat1, mat2, n_subsample, subsampling_ratio, parallel = parallel)
  })
  
  # construct starting dissimilarity matrix
  n_groups <- length(groups)
  dissMat <- matrix(NA, nrow = n_groups, ncol = n_groups)
  dissMatNames <- seq(n_groups)
  for (i in seq(1, n_groups - 1, 1)) {
    for (j in seq(i + 1, n_groups, 1)) {
      dissMat[i, j] <- matrixDissimilarity(bgraphs[[i]], bgraphs[[j]])
      dissMat[j, i] <- dissMat[i, j]
    }
  }
  
  # keep track of (1) group index for creating hclust object for dendrogram plotting,
  # (2) height of merged groups (starts at 0 for each group)
  groupIdx <- list()
  nodeHeightsLst <- list()
  mergeGroups <- list()
  for (i in seq(n_groups)) {
    groupIdx[[i]] <- -i
    nodeHeightsLst[[i]] <- 0
    mergeGroups[[i]] <- groupNames[i]
  }
  mergeMat <- matrix(0, nrow = n_groups - 1, ncol = 2)
  
  # hierarchical clustering
  n_nodes <- n_groups - 1
  nodeMemberships <- list()
  nodePvals <- vector(length = n_nodes)
  nodeHeights <- vector(length = n_nodes)
  nodeSCCA <- list()
  nodeGroups <- list()
  index <- 1
  
  while (length(dissMat) > 0) {
    rowCol <- whichMerge(dissMat)
    row <- rowCol[["row"]]
    col <- rowCol[["col"]]
    
    # merge: new group, new graph
    newGroupMembers <- c(groups[[row]], groups[[col]])
    nodeMemberships[[index]] <- newGroupMembers
    nodeHeightsLst[[n_groups + index]] <- max(nodeHeightsLst[[dissMatNames[row]]], nodeHeightsLst[[dissMatNames[col]]]) + dissMat[row, col]
    nodeHeights[index] <- nodeHeightsLst[[n_groups + index]] 
    newGraph <- constructBipartiteGraph(X[newGroupMembers, ], Y[newGroupMembers, ], 
      n_subsample, subsampling_ratio, parallel = parallel)

    # standardize columns of X, Y
    Xnew = X[newGroupMembers, ]
    Ynew = Y[newGroupMembers, ]

    Xnew <- scale(Xnew)
    Ynew <- scale(Ynew)
      
    Xnew[is.nan(Xnew)] = 0
    Ynew[is.nan(Ynew)] = 0

    nodeSCCA[[index]] = scca::scca(Xnew, Ynew, penalty = "LASSO")
    nodeGroups[[index]] <- c(mergeGroups[[row]], mergeGroups[[col]])
    
    # calculate p-values
    if (p.value) {
      null_distribution <- null_distri(X[groups[[row]], ], Y[groups[[row]], ], 
        X[groups[[col]], ], Y[groups[[col]], ], n.perm = n_perm, parallel = parallel)
      nodePvals[index] <- p_value(dissMat[row, col], null_distribution)
    }

    print(paste0("Merge ", index, "/", n_nodes, " complete"))
    
    # update: dissimilarity matrix, group membership, and list of weighted bipartite graphs, groupIdx
    groups[[row]] <- NULL
    groups[[col]] <- NULL
    groups[[length(groups) + 1]] <- newGroupMembers
    names(groups) <- seq(length(groups))
    
    mergeMat[index, 1] <- groupIdx[[row]]
    mergeMat[index, 2] <- groupIdx[[col]]
    groupIdx[[row]] <- NULL
    groupIdx[[col]] <- NULL
    groupIdx[[length(groupIdx) + 1]] <- index
    names(groupIdx) <- seq(length(groupIdx))

    mergeGroups[[length(mergeGroups) + 1]] <- c(mergeGroups[[row]], mergeGroups[[col]])
    mergeGroups[[row]] <- NULL
    mergeGroups[[col]] <- NULL
    names(mergeGroups) <- seq(length(mergeGroups))
    
    bgraphs[[row]] <- NULL
    bgraphs[[col]] <- NULL
    bgraphs[[length(bgraphs) + 1]] <- newGraph
    names(bgraphs) <- seq(length(bgraphs))
    dissMatNames <- c(dissMatNames[-c(row, col)], n_groups + index)
    dissMat <- dissMat[-c(row, col), -c(row, col)]

    if (length(dissMat) > 0) {
      newDist <- vector(length = sqrt(length(dissMat)))
      for (i in seq(sqrt(length(dissMat)))) {
        newDist[i] <- matrixDissimilarity(bgraphs[[i]], newGraph)
      }
      dissMat <- cbind(dissMat, newDist)
      dissMat <- rbind(dissMat, c(newDist, NA))
      
      index <- index + 1
    }
  }
  
  # create hclust object
  hclustObj <- createHclust(groupNames, mergeMat, nodeHeights)
  retLst <- list(nodeMemberships = nodeMemberships, hclustObj = hclustObj, nodeSCCA = nodeSCCA, nodeGroups = nodeGroups)

  if (p.value) {
    retLst[["nodePvals"]] <- nodePvals
  }
  
  return(retLst)
}


#' Construct Bipartite Graph
#'
#' Constructs matrix B containing information of bipartite relationship between mat1 and mat2 (see paper)
#'
#' @param mat1 an n x p matrix (e.g. for gene expression)
#' @param mat2 an n x q matrix (e.g. for drug sensitivity)
#' @param n_subsample number of times to perform subsampling to generate B
#' @param subsampling_ratio fraction of samples to subsample each time
#' @param parallel boolean for whether to parallelize subsampling
#' @return a p x q matrix containing information of bipartite relationship
#' @export
constructBipartiteGraph <- function(mat1, mat2, n_subsample = 1, subsampling_ratio = 1, parallel = TRUE) {
  # Constructs matrix B containing information of bipartite relationship between mat1 and mat2 (see paper)
  # Input:
  #   mat1: an n x p matrix (e.g. for gene expression)
  #   mat2: an n x q matrix (e.g. for drug sensitivity)
  #   n_subsample: number of times to perform subsampling to generate B
  #   subsampling_ratio: fraction of samples to subsample each time
  #   parallel: boolean for whether to parallelize subsampling
  # Output:
  #   B: a p x q matrix containing information of bipartite relationship

  n_samples <- nrow(mat1)
  p <- ncol(mat1)
  q <- ncol(mat2)

  n_cores <- parallel::detectCores() - 1
  if (parallel && n_cores > 0 && n_subsample > 1) {
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterExport(cl, c("mat1", "mat2", "n_samples", "p", "q", 
      "subsampling_ratio"), envir = environment())
    parallel::clusterSetRNGStream(cl = cl)

    Bmatrices <- parallel::parLapply(cl = cl, seq(n_subsample), function (x) {
      sample_ids <- sample(seq(n_samples), round(n_samples * subsampling_ratio), replace = FALSE)

      # standardize columns of X, Y
      mat1_sub <- scale(mat1[sample_ids, ])
      mat2_sub <- scale(mat2[sample_ids, ])
      
      mat1_sub[is.nan(mat1_sub)] = 0
      mat2_sub[is.nan(mat2_sub)] = 0

      # don't use center, scale in scca() because if an entire column is 0, NaNs will be produced
      scca_rlt <- scca::scca(mat1_sub, mat2_sub, penalty = "LASSO") # might need to tune parameters
      
      vec_p <- scca_rlt$A[, 1]
      vec_q <- scca_rlt$B[, 1]

      return(outer(vec_p, vec_q))
    })
    Bsum <- Reduce("+", Bmatrices)
    B <- Bsum / length(Bmatrices)
    parallel::stopCluster(cl)
  } else {
    # sub-sampling the sampling in each group to get a weighted bipartite graph
    B <- replicate(n_subsample, {
      sample_ids <- sample(seq(n_samples), round(n_samples * subsampling_ratio), replace = FALSE)
      
      # standardize columns of X, Y
      mat1_sub <- scale(mat1[sample_ids, ])
      mat2_sub <- scale(mat2[sample_ids, ])
      
      mat1_sub[is.nan(mat1_sub)] = 0
      mat2_sub[is.nan(mat2_sub)] = 0
      
      # don't use center, scale in scca() because if an entire column is 0, NaNs will be produced
      scca_rlt <- scca::scca(mat1_sub, mat2_sub, penalty = "LASSO") # might need to tune parameters
      
      vec_p <- scca_rlt$A[, 1]
      vec_q <- scca_rlt$B[, 1]
      
      return(c(outer(vec_p, vec_q)))
      }) %>% rowMeans() %>% matrix(nrow = p, ncol = q)
  }
  return(B)
}

#' Matrix dissimilarity
#'
#' Computes nucleus norm-based dissimilarity measure between two matrices
#'
#' @param B1 a p x q matrix containing information of bipartite relationship
#' @param B2 a p x q matrix containing information of bipartite relationship
#' @return nucleus norm-based dissimilarity
#' @export
matrixDissimilarity <- function(B1, B2){
  # Computes nucleus norm-based dissimilarity measure between two matrices
  # Input:
  #   B1: a p x q matrix containing information of bipartite relationship
  #   B2: a p x q matrix containing information of bipartite relationship
  # Output:
  #   nucleus norm-based dissimilarity
  return(sum(svd(B1 - B2)$d) / (sum(svd(B1)$d) + sum(svd(B2)$d)))
}

#' Null distribution of dissimilarity measures
#'
#' Generates null distribution of dissimilarity measures between (X1, Y1) and (X2, Y2)
#'
#' @param X1 an n x p matrix (e.g. for gene expression) for group 1
#' @param Y1 an n x q matrix (e.g. for drug sensitivity) for group 1
#' @param X2 an n x p matrix (e.g. for gene expression) for group 2
#' @param Y2 an n x q matrix (e.g. for drug sensitivity) for group 2
#' @param n.perm number of null dissimilarity measures to generate
#' @param parallel boolean for whether to parallelize subsampling
#' @return vector of length n.perm of null dissimilarity measures
#' @export
null_distri <- function(X1, Y1, X2, Y2, n.perm = 100, parallel = TRUE) {
  # Generates null distribution of dissimilarity measures between (X1, Y1) and (X2, Y2)
  # Input:
  #   X1: an n x p matrix (e.g. for gene expression) for group 1
  #   Y1: an n x q matrix (e.g. for drug sensitivity) for group 1
  #   X1: an n x p matrix (e.g. for gene expression) for group 2
  #   Y1: an n x q matrix (e.g. for drug sensitivity) for group 2
  #   n.perm: number of null dissimilarity measures to generate
  #   parallel: boolean for whether to parallelize subsampling
  # Output:
  #   dissimilarities: vector of length n.perm of null dissimilarity measures
  X1 <- scale(X1)
  X1[is.nan(X1)] = 0
  Y1 <- scale(Y1)
  Y1[is.nan(Y1)] = 0
  
  X2 <- scale(X2)
  X2[is.nan(X2)] = 0
  Y2 <- scale(Y2)
  Y2[is.nan(Y2)] = 0

  n_cores <- parallel::detectCores() - 1
  if (parallel && n_cores > 0) {
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterExport(cl, c("X1", "Y1", "X2", "Y2", "matrixDissimilarity"),
      envir = environment())
    parallel::clusterSetRNGStream(cl = cl)
  
    dissimilarities <- parallel::parSapply(cl = cl, seq(n.perm), function (x) {
      sampleIdx1 <- sample(seq(nrow(X1)), size = nrow(X1), replace = FALSE)
      sampleIdx2 <- sample(seq(nrow(X2)), size = nrow(X2), replace = FALSE)
      
      X1.b <- X1[sampleIdx1, ]
      X2.b <- X2[sampleIdx2, ]
      
      scca1 <- scca::scca(X1.b, Y1, penalty = "LASSO")
      scca2 <- scca::scca(X2.b, Y2, penalty = "LASSO")
      
      B1 <- outer(scca1$A, scca1$B)
      B2 <- outer(scca2$A, scca2$B)
      
      return(matrixDissimilarity(B1, B2))
    }, simplify = TRUE)
    parallel::stopCluster(cl)
  } else {
    dissimilarities <- replicate(n.perm, {
      sampleIdx1 <- sample(seq(nrow(X1)), size = nrow(X1), replace = FALSE)
      sampleIdx2 <- sample(seq(nrow(X2)), size = nrow(X2), replace = FALSE)
    
      X1.b <- X1[sampleIdx1, ]
      X2.b <- X2[sampleIdx2, ]
      
      scca1 <- scca::scca(X1.b, Y1, penalty = "LASSO")
      scca2 <- scca::scca(X2.b, Y2, penalty = "LASSO")
      
      B1 <- outer(scca1$A, scca1$B)
      B2 <- outer(scca2$A, scca2$B)
      
      return(matrixDissimilarity(B1, B2))
      }, simplify = TRUE)
  }
  return(dissimilarities)
}

#' P-value
#'
#' Computes p-value as number of null dissimilarities less than or equal to observed dissimilarity
#'
#' @param dissimilarity observed dissimilarity
#' @param dissimilarities null distribution of dissimilarities
#' @return p-value
#' @export
p_value <- function(dissimilarity, dissimilarities) {
  # Computes p-value as number of null dissimilarities less than or equal to observed dissimilarity
  # Input:
  #   dissimilarity: observed dissimilarity
  #   dissimilarities: null distribution of dissimilarities
  # Output:
  #   p-value
  return(mean(dissimilarities <= dissimilarity))
}

#' Significant Merges
#'
#' Filters bipartite hierarchical clustering merged groups by p-value
#'
#' @param results list of results from bipartite hierarchical clustering
#' @param p p-value to filter merged groups by
#' @return list of results from bipartite hierarchical clustering filtered by p-value
#' @export
getSignificantMergedGroups <- function(results, p = 0.05) {
  # filters bipartite hierarchical clustering merged groups by p-value
  # Input:
  #   results: list of results from bipartite hierarchical clustering
  #   p: p-value to filter merged groups by
  # Output:
  #   retLst; list of results from bipartite hierarchical clustering filtered by p-value
  if (!"nodePvals" %in% names(results)) {
    print("p-value must be computed first!")
  } else {
    nodeMemberships <- results[["nodeMemberships"]]
    nodeSCCA <- results[["nodeSCCA"]]
    nodePvals <- results[["nodePvals"]]
    nodeGroups <- results[["nodeGroups"]]

    retLst <- list()
    nodeMembershipsFiltered <- list()
    nodeSCCAFiltered <- list()
    nodePvalsFiltered <- list()
    nodeGroupsFiltered <- list()

    n <- length(nodePvals)
    index <- 1
    for (i in seq(n)) {
      if (nodePvals[i] <= p) {
        nodeMembershipsFiltered[[toString(index)]] <- nodeMemberships[[i]]
        nodeSCCAFiltered[[toString(index)]] <- nodeSCCA[[i]]
        nodePvalsFiltered[[toString(index)]] <- nodePvals[[i]]
        nodeGroupsFiltered[[toString(index)]] <- nodeGroups[[i]]
      }
      index <- index + 1
    }
    retLst <- list("nodeMemberships" = nodeMembershipsFiltered, "nodeSCCA" = nodeSCCAFiltered,
      "nodePvals" = nodePvalsFiltered, "hclustObj" = results[["hclustObj"]], "nodeGroups" = nodeGroupsFiltered)
    return(retLst)
  }
}



