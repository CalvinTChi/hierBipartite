#' Given index from hclust merge matrix, return X, Y row indices corresponding to groups specified by index
#'
#' @param idx index from hclust merge
#' @param groups list of starting group membership (input to hierBipartite())
#' @param groupMerges list of merged groups thus far (groupMerge[[3]] is a vector of group names from cluster created at third merge)
#' @return vector of row indices corresponding to X, Y
#' @export
getMergeGroupRows <- function(idx, groups, groupMerges) {
  # Given index from hclust merge matrix, return X, Y row indices corresponding to groups specified by index
  # Input:
  #   idx: index from hclust merge
  #   groups: list of starting group membership (input to hierBipartite())
  #   groupMerges: list of merged groups thus far (groupMerge[[3]] is a vector of group names from cluster created at third merge)
  # Output:
  #   rows: vector of row indices corresponding to X, Y
  if (idx < 0) {
    return(groups[[-idx]])
  } 
  rows = c()
  for (group in groupMerges[[idx]]) {
    rows = c(rows, groups[[group]])
  }
  return(rows)
}

#' Given group indices for a merge from hclust merge matrix, return group names of merged cluster
#'
#' @param idx1 group index of first group from hclust merge
#' @param idx2 group index of second group from hclust merge
#' @param groups list of starting group membership (input to hierBipartite())
#' @param groupMerges list of merged groups thus far (groupMerge[[3]] is a vector of group names from cluster created at third merge)
#' @return vector of group names after merge of group indicated by idx1 and idx2
#' @export
newMergedGroup <- function(idx1, idx2, groups, groupMerges) {
  # Given group indices for a merge from hclust merge matrix, return group names of merged cluster
  # Input:
  #   idx1: group index of first group from hclust merge
  #   idx2: group index of second group from hclust merge
  #   groups: list of starting group membership (input to hierBipartite())
  #   groupMerges: list of merged groups thus far (groupMerge[[3]] is a vector of group names from cluster created at third merge)
  # Output:
  #   newestGroup: vector of group names after merge of group indicated by idx1 and idx2

  groupNames = names(groups)

  newestGroup = c()
  # update first group
  if (idx1 < 0) {
    newestGroup = c(newestGroup, groupNames[-idx1])
  } else {
    newestGroup = c(newestGroup, groupMerges[[idx1]])
  }

  # update second group
  if (idx2 < 0) {
    newestGroup = c(newestGroup, groupNames[-idx2])
  } else {
    newestGroup = c(newestGroup, groupMerges[[idx2]])
  }

  return(newestGroup)
}

#' Bipartite Hierarchical Clustering
#'
#' Main bipartite hierarchial clustering algorithm. Execute browseVignettes("hierBipartite") or visit \href{https://calvintchi.github.io/html/hierBipartite}{here} for vignette on using the
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
#'         (1) hclustObj: dendrogram class of resulting dendrogram, 
#'         (2) groupMerges: list of groups for each merge, in order of merge,
#'         (3) nodePvals: list of p-value of each new merge, in order of merge if p.value = TRUE.
#' @export
hierBipartite <- function(X, Y, groups, link = "ward.D2", n_subsample = 1, subsampling_ratio = 1, p.value = FALSE, n_perm = 100, parallel = FALSE) {
  # Main bipartite hierarchial clustering algorithm
  # Input:
  #   X: an n x p matrix (e.g. for gene expression)
  #   Y: an n x q matrix (e.g. for drug sensitivity)
  #   link: string indicating link function as input to hclust. One of "ward.D", "ward.D2", "single", "complete", "average",
  #         "mcquitty", "median", "centroid".
  #   groups: a list of starting group membership (e.g. list("1" = c(1,2,3), "2" = c(4,5,6)) means group 1 has 
  #           samples 1, 2, 3, and group 2 has samples 4, 5, 6.
  #   n_subsample: number of subsampling to generate matrix B (see paper)
  #   subsampling_ratio: fraction of samples to sample for subsampling to generate matrix B (see paper)
  #   p.value: boolean for whether to generate p-values for each merge
  #   n_perm: number of permutations for generating p-values. Ignored if p.value = FALSE
  #   parallel: boolean for whether to parallelize subsampling and p-value generation step
  # Output:
  #   retLst: list of results from bipartite hierarchical clustering, containing
  #     hclustObj: dendrogram class of resulting dendrogram
  #     groupMerges: list of groups for each merge, in order of merge
  #     nodePvals: list p-value of each new merge, in order of merge if p.value = TRUE
  groupNames <- names(groups)
  
  print("Computing starting dissimilarity matrix...")
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
  hclustObj = hclust(as.dist(dissMat), method = link)
  hclustObj$labels = names(groups)
  retLst <- list(hclustObj = hclustObj)

  # create groupMerges
  groupMerges = list()
  merge = hclustObj$merge
  n_merges = nrow(merge)
  for (i in seq(n_merges)) {
    groupMerges[[length(groupMerges) + 1]] = newMergedGroup(merge[i, 1], merge[i, 2], groups, groupMerges)
  }
  retLst[["groupMerges"]] = groupMerges

  # p-value
  if (p.value) {
    nodePvals = list()
    for (i in seq(n_merges)) {
      rows1 = getMergeGroupRows(merge[i, 1], groups, groupMerges)
      rows2 = getMergeGroupRows(merge[i, 2], groups, groupMerges)

      X1 = X[rows1, ]
      Y1 = Y[rows1, ]
      X2 = X[rows2, ]
      Y2 = Y[rows2, ]

      if (merge[i, 1] < 0 && merge[i, 2] < 0) {
        d = dissMat[-merge[i, 1], -merge[i, 2]]
      } else {
        B1 = constructBipartiteGraph(X1, Y1, n_subsample = 1, 
          subsampling_ratio = 1, parallel = FALSE)
        B2 = constructBipartiteGraph(X2, Y2, n_subsample = 1, 
          subsampling_ratio = 1, parallel = FALSE)
        d = matrixDissimilarity(B1, B2)
      }

      dissimilarities = null_distri(X1, Y1, X2, Y2)
      nodePvals[[i]] = p_value(d, dissimilarities)

      print(paste0("P-value for merge ", i, "/", n_merges, " complete"))
    }
    retLst[["nodePvals"]] = nodePvals
  }
  
  return(retLst)
}

#' Only scales features (columns) with positive variance
#'
#' @param mat an n x p matrix (gene expression or drug sensitivity)
#' @return an n x p matrix with columns with positive variance scaled 
#' @export
scale_features <- function(mat) {
  # Only scales features (columns) with positive variance
  # Input:
  #   mat: an n x p matrix (gene expression or drug sensitivity)
  # Output:
  #   mat: an n x p matrix with columns with positive variance scaled
  vars <- apply(mat, 2, var)
  idx <- which(vars > 0)
  mat[, idx] <- scale(mat[, idx])
  idx <- which(vars == 0)
  mat[, idx] <- 0
  return(mat)
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

  # if n_samples <= 5, cannot perform CV from scca(), so do not subsample
  if (n_samples <= 5) {
    mat1 <- scale_features(mat1)
    mat2 <- scale_features(mat2)
    scca_rlt <- scca::scca(mat1, mat2, penalty = "LASSO")
    return(outer(c(scca_rlt$A), c(scca_rlt$B)))
  }

  n_cores <- parallel::detectCores() - 1
  if (parallel && n_cores > 0) {
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterExport(cl, c("mat1", "mat2", "n_samples", "p", "q", 
      "subsampling_ratio"), envir = environment())
    parallel::clusterSetRNGStream(cl = cl)

    Bmatrices <- parallel::parLapply(cl = cl, seq(n_subsample), function (x) {
      sample_ids <- sample(seq(n_samples), round(n_samples * subsampling_ratio), replace = FALSE)

      mat1_sub <- mat1[sample_ids, ]
      mat2_sub <- mat2[sample_ids, ]

      # only standardize features with positive variance
      mat1_sub <- scale_features(mat1_sub)
      mat2_sub <- scale_features(mat2_sub)

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

      mat1_sub <- mat1[sample_ids, ]
      mat2_sub <- mat2[sample_ids, ]

      # only standardize features with positive variance
      mat1_sub <- scale_features(mat1_sub)
      mat2_sub <- scale_features(mat2_sub)
      
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
  return(sum(irlba::irlba(B1 - B2)$d) / (sum(irlba::irlba(B1)$d) + sum(irlba::irlba(B2)$d)))
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
  X1 = scale_features(X1)
  Y1 = scale_features(Y1)
  
  X2 = scale_features(X2)
  Y2 = scale_features(Y2)

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
      
      B1 <- outer(c(scca1$A), c(scca1$B))
      B2 <- outer(c(scca2$A), c(scca2$B))
      
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
      
      B1 <- outer(c(scca1$A), c(scca1$B))
      B2 <- outer(c(scca2$A), c(scca2$B))
      
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
    nodePvals <- results[["nodePvals"]]
    groupMerges <- results[["groupMerges"]]

    retLst <- list()
    nodePvalsFiltered <- list()
    groupMergesFiltered <- list()

    n <- length(nodePvals)
    index <- 1
    for (i in seq(n)) {
      if (nodePvals[i] <= p) {
        nodePvalsFiltered[[toString(index)]] <- nodePvals[[i]]
        groupMergesFiltered[[toString(index)]] <- groupMerges[[i]]
      }
      index <- index + 1
    }
    retLst <- list("nodePvals" = nodePvalsFiltered, "hclustObj" = results[["hclustObj"]], "groupMerges" = groupMergesFiltered)
    return(retLst)
  }
}



