#' Bipartite Graph-based Hierarchical Clustering
#'
#' Main bipartite graph-based hierarchial clustering algorithm. Visit \href{https://calvintchi.github.io/html/hierBipartite}{here} for vignette on using the
#' hierBipartite package.
#' @importFrom stats var hclust as.dist median sd
#'
#' @param X an n x p matrix of variable set 1 (e.g. gene expression)
#' @param Y an n x q matrix of variable set 2 (e.g. drug sensitivity)
#' @param groups a list of starting group membership (e.g. list("1" = c(1,2,3), "2" = c(4,5,6)) means group 1 has samples 1, 2, 3, and group 2 has samples 4, 5, 6.
#' @param link string indicating link function as input to hclust(). One of "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid".
#' @param n_subsample number of subsampling to generate matrix B (see paper)
#' @param subsampling_ratio fraction of samples to sample for subsampling to generate matrix B (see paper)
#' @param p.value boolean for whether to generate p-values for each merge
#' @param n_perm number of permutations for generating p-values. Ignored if p.value = FALSE
#' @param parallel boolean for whether to parallelize subsampling and p-value generation step
#' @param maxCores maximum number of cores to use (only applicable when parallel = TRUE)
#' @param p_cutoff p-value cutoff that determines whether merge is significant. If p-value > p_cutoff, p-values will not be calculated for future merges involving current group. Ignored if p.value = FALSE.
#' @return list of results from bipartite graph-based hierarchical clustering, containing up to
#' \itemize{
#'   \item hclustObj: hclust object
#'   \item groupMerges: list of clusters after each merge, in order of merge.
#'   Each cluster is indicated by a vector of cell line groups
#'   \item nodePvals: list of p-value of each new merge,
#'   in order of merge. Only available if p.value = TRUE
#'   \item D: dissimilarity matrix
#' }
#' @examples
#' # Get a small subset of the test dataset
#' data(ctrp2)
#'
#' groups = ctrp2$groups
#' X = ctrp2$X
#' Y = ctrp2$Y
#'
#' groupNames = names(groups)
#' groupSmall = groups[groupNames[1:3]]
#'
#' \dontrun{
#' # Basic call of hierBipartite() on small test dataset
#' result0 = hierBipartite(X, Y, groupSmall)
#'
#' # Calling hierBipartite() with subsampling
#' result1 = hierBipartite(X, Y, groupSmall, n_subsample = 100, subsampling_ratio = 0.90)
#'
#' # Calling hierBipartite() with p-value generation
#' result2 = hierBipartite(X, Y, groupSmall, n_perm = 100, p.value = TRUE, p_cutoff = 0.10)
#'
#' # Calling hierBipartite() with both subsampling and p-value generation (expensive)
#' result3 = hierBipartite(X, Y, groupSmall, n_subsample = 100, subsampling_ratio = 0.90,
#'                         n_perm = 100, p.value = TRUE, p_cutoff = 0.10)
#' }
#'
#' @export
hierBipartite <- function(X, Y, groups, link = "ward.D2", n_subsample = 1, subsampling_ratio = 1, p.value = FALSE,
                          n_perm = 100, parallel = FALSE, maxCores = 7, p_cutoff = 0.10) {
  # Main bipartite graph-based hierarchial clustering algorithm
  # Input:
  #   X: an n x p matrix of variable set 1 (e.g. gene expression)
  #   Y: an n x q matrix of variable set 2 (e.g. drug sensitivity)
  #   link: string indicating link function as input to hclust. One of "ward.D", "ward.D2", "single", "complete", "average",
  #         "mcquitty", "median", "centroid".
  #   groups: a list of starting group membership (e.g. list("1" = c(1,2,3), "2" = c(4,5,6)) means group 1 has
  #           samples 1, 2, 3, and group 2 has samples 4, 5, 6.
  #   n_subsample: number of subsampling to generate matrix B (see paper)
  #   subsampling_ratio: fraction of samples to sample for subsampling to generate matrix B (see paper)
  #   p.value: boolean for whether to generate p-values for each merge
  #   n_perm: number of permutations for generating p-values. Ignored if p.value = FALSE
  #   parallel: boolean for whether to parallelize subsampling and p-value generation step
  #   maxCores: maximum number of cores to use (only applicable when parallel = TRUE)
  #   p_cutoff: p-value cutoff that determines whether merge is significant. If p-value > p_cutoff, p-values will not be
  #             calculated for future merges involving current group.
  # Output:
  #   retLst: list of results from bipartite graph-based hierarchical clustering, containing up to
  #     hclustObj: hclust object
  #     groupMerges: list of groups for each merge, in order of merge
  #     nodePvals: list p-value of each new merge, in order of merge. Only available if p.value = TRUE
  #     D: dissimilarity matrix
  groupNames <- names(groups)

  print("Computing starting dissimilarity matrix...")
  # construct matrix B representing bipartite relationship
  bgraphs <- lapply(groups, function(gid) {
    mat1 <- X[gid, ]
    mat2 <- Y[gid, ]

    constructBipartiteGraph(mat1, mat2, n_subsample, subsampling_ratio, parallel = parallel, maxCores = maxCores)
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
  hclustObj <- hclust(as.dist(dissMat), method = link)
  hclustObj$labels <- names(groups)
  retLst <- list(hclustObj = hclustObj)
  retLst[["D"]] <- dissMat

  # create groupMerges
  groupMerges <- list()
  merge <- hclustObj$merge
  n_merges <- nrow(merge)
  for (i in seq(n_merges)) {
    groupMerges[[length(groupMerges) + 1]] <- newMergedGroup(merge[i, 1], merge[i, 2], groups, groupMerges)
  }
  retLst[["groupMerges"]] <- groupMerges

  # p-value
  if (p.value) {
    nodePvals <- list()
    for (i in seq(n_merges)) {
      print(paste0("Compute p-value for merge ", i, "/", n_merges))
      m1 <- merge[i, 1]
      m2 <- merge[i, 2]
      if ((m1 > 0 && (nodePvals[[m1]] > p_cutoff || is.na(nodePvals[[m1]]))) || (m2 > 0 && (nodePvals[[m2]] > p_cutoff || is.na(nodePvals[[m2]])))) {
        nodePvals[[i]] <- NA
      } else {
        rows1 <- getMergeGroupRows(m1, groups, groupMerges)
        rows2 <- getMergeGroupRows(m2, groups, groupMerges)

        X1 <- X[rows1, ]
        Y1 <- Y[rows1, ]
        X2 <- X[rows2, ]
        Y2 <- Y[rows2, ]

        if (merge[i, 1] < 0 && merge[i, 2] < 0) {
          d <- dissMat[-m1, -m2]
        } else {
          B1 <- constructBipartiteGraph(X1, Y1, n_subsample = 1, subsampling_ratio = 1, parallel = FALSE)
          B2 <- constructBipartiteGraph(X2, Y2, n_subsample = 1, subsampling_ratio = 1, parallel = FALSE)
          d <- matrixDissimilarity(B1, B2)
        }

        dissimilarities <- null_distri(X1, Y1, X2, Y2, parallel = parallel, maxCores = maxCores)
        nodePvals[[i]] <- p_value(d, dissimilarities)
      }
      print(paste0("P-value for merge ", i, "/", n_merges, " complete"))
    }
    retLst[["nodePvals"]] = nodePvals
  }

  return(retLst)
}

#' Construct Bipartite Graph Edge Weight Matrix of Gene-drug Association Patterns
#'
#' Constructs edge weight matrix B representing association between set of variables in mat1 and set of variables in mat2 (see paper).
#'
#' @importFrom magrittr %>%
#' @param mat1 an n x p matrix of variable set 1 (e.g. gene expression)
#' @param mat2 an n x q matrix of variable set 2 (e.g. drug sensitivity)
#' @param n_subsample number of times to perform subsampling to generate B
#' @param subsampling_ratio fraction of samples to subsample each time
#' @param parallel boolean for whether to parallelize subsampling
#' @param maxCores maximum number of cores to use (only applicable when parallel = TRUE)
#' @return a p x q matrix of bipartite graph edge weights
#'
#' @examples
#' # Extract bipartite edge weight matrix B for cell lines from the
#' # squamous cell carcinoma, esophagus group
#' data(ctrp2)
#'
#' groups = ctrp2$groups
#' X = ctrp2$X
#' Y = ctrp2$Y
#'
#' x = X[groups[["squamous_cell_carcinoma_esophagus"]], ]
#' y = Y[groups[["squamous_cell_carcinoma_esophagus"]], ]
#'
#' # Extract bipartite edge weight matrix B with subsampling
#' \dontrun{
#' B = constructBipartiteGraph(x, y, n_subsample = 100,
#'                             subsampling_ratio = 0.90,
#'                             parallel = TRUE, maxCores = 2)
#' }
#'
#' @export
constructBipartiteGraph <- function(mat1, mat2, n_subsample = 1,
                                    subsampling_ratio = 1,
                                    parallel = FALSE,
                                    maxCores = 7) {
  # Constructs matrix B containing information of bipartite relationship between mat1 and mat2 (see paper)
  # Input:
  #   mat1: an n x p matrix of variable set 1 (e.g. gene expression)
  #   mat2: an n x q matrix of variable set 2 (e.g. drug sensitivity)
  #   n_subsample: number of times to perform subsampling to generate B
  #   subsampling_ratio: fraction of samples to subsample each time
  #   parallel: boolean for whether to parallelize subsampling
  #   maxCores: maximum number of cores to use (only applicable when parallel = TRUE)
  # Output:
  #   B: a p x q matrix of bipartite graph edge weights
  n_samples <- nrow(mat1)
  p <- ncol(mat1)
  q <- ncol(mat2)

  # if n_samples <= 5, cannot perform CV from scca(), so do not subsample
  if (n_samples <= 5) {
    mat1 <- scale_features(mat1)
    mat2 <- scale_features(mat2)
    scca_rlt <- scca(mat1, mat2, penalty = "LASSO")
    return(outer(c(scca_rlt$A), c(scca_rlt$B)))
  }

  n_cores <- min(parallel::detectCores() - 1, maxCores)
  if (parallel && n_cores > 1 && n_subsample > 1) {
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
      scca_rlt <- scca(mat1_sub, mat2_sub, penalty = "LASSO") # might need to tune parameters

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
      scca_rlt <- scca(mat1_sub, mat2_sub, penalty = "LASSO") # might need to tune parameters

      vec_p <- scca_rlt$A[, 1]
      vec_q <- scca_rlt$B[, 1]

      return(c(outer(vec_p, vec_q)))
      }) %>% rowMeans() %>% matrix(nrow = p, ncol = q)
  }
  return(B)
}

#' Matrix dissimilarity
#'
#' Computes nuclear norm-based dissimilarity measure between two matrices.
#'
#' @param B1 first p x q bipartite graph edge weight matrix
#' @param B2 second p x q bipartite graph edge weight matrix
#' @return nuclear norm-based dissimilarity
#' @examples
#' # Compute matrix dissimilarity in edge weight matrix between squamous cell
#' # carcinoma, esophagus and squamous cell carcinoma, upper aerodigestive
#' data(ctrp2)
#'
#' groups = ctrp2$groups
#' X = ctrp2$X
#' Y = ctrp2$Y
#'
#' x1 = X[groups[["squamous_cell_carcinoma_esophagus"]], ]
#' y1 = Y[groups[["squamous_cell_carcinoma_esophagus"]], ]
#'
#' \dontrun{
#' B1 = constructBipartiteGraph(x1, y1)
#' }
#'
#' x2 = X[groups[["squamous_cell_carcinoma_upper_aerodigestive"]], ]
#' y2 = Y[groups[["squamous_cell_carcinoma_upper_aerodigestive"]], ]
#'
#' \dontrun{
#' B2 = constructBipartiteGraph(x2, y2)
#' matrixDissimilarity(B1, B2)
#' }
#'
#' @export
matrixDissimilarity <- function(B1, B2){
  # Computes nuclear norm-based dissimilarity measure between two bipartite graph edge weight matrices
  # Input:
  #   B1: first p x q bipartite graph edge weight matrix
  #   B2: second p x q bipartite graph edge weight matrix
  # Output:
  #   nuclear norm-based dissimilarity
  return(sum(irlba::irlba(B1 - B2, nv = 2)$d) / (sum(irlba::irlba(B1, nv = 1)$d) + sum(irlba::irlba(B2, nv = 1)$d)))
}

#' Null distribution of dissimilarity measures
#'
#' Generates null distribution of dissimilarity measures between group 1 (X1, Y1) and group 2 (X2, Y2).
#'
#' @param X1 an n x p matrix of variable set 1 (e.g. gene expression) from group 1
#' @param Y1 an n x q matrix of variable set 2 (e.g. drug sensitivity) from group 1
#' @param X2 an n x p matrix of variable set 1 (e.g. gene expression) from group 2
#' @param Y2 an n x q matrix of varaible set 2 (e.g. drug sensitivity) from group 2
#' @param n.perm number of null dissimilarity measures to generate
#' @param parallel boolean for whether to parallelize permutation
#' @param maxCores maximum number of cores to use (only applicable when parallel = TRUE)
#' @return vector of length n.perm of null dissimilarity measures
#' @examples
#' # Get data for group squamous cell carcinoma, esophagus and for group
#' # squamous cell carcinoma, upper aerodigestive
#' data(ctrp2)
#'
#' groups = ctrp2$groups
#' X = ctrp2$X
#' Y = ctrp2$Y
#'
#' x1 = X[groups[["squamous_cell_carcinoma_esophagus"]], ]
#' y1 = Y[groups[["squamous_cell_carcinoma_esophagus"]], ]
#'
#' x2 = X[groups[["squamous_cell_carcinoma_upper_aerodigestive"]], ]
#' y2 = Y[groups[["squamous_cell_carcinoma_upper_aerodigestive"]], ]
#'
#' \dontrun{
#' dissimilarities = null_distri(x1, y1, x2, y2, n.perm = 100)
#' }
#'
#' @export
null_distri <- function(X1, Y1, X2, Y2, n.perm = 100, parallel = FALSE, maxCores = 7) {
  # Generates null distribution of dissimilarity measures between (X1, Y1) and (X2, Y2)
  # Input:
  #   X1: an n x p matrix of variable set 1 (e.g. gene expression) from group 1
  #   Y1: an n x q matrix of variable set 2 (e.g. drug sensitivity) from group 1
  #   X2: an n x p matrix of variable set 1 (e.g. gene expression) from group 2
  #   Y2: an n x q matrix of varaible set 2 (e.g. drug sensitivity) from group 2
  #   n.perm: number of null dissimilarity measures to generate
  #   parallel: boolean for whether to parallelize permutation
  #   maxCores: maximum number of cores to use (only applicable when parallel = TRUE)
  # Output:
  #   dissimilarities: vector of length n.perm of null dissimilarity measures
  X1 = scale_features(X1)
  Y1 = scale_features(Y1)

  X2 = scale_features(X2)
  Y2 = scale_features(Y2)

  n_cores <- min(parallel::detectCores() - 1, maxCores)
  if (parallel && n_cores > 1) {
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterExport(cl, c("X1", "Y1", "X2", "Y2", "matrixDissimilarity"),
      envir = environment())
    parallel::clusterSetRNGStream(cl = cl)

    dissimilarities <- parallel::parSapply(cl = cl, seq(n.perm), function (x) {
      sampleIdx1 <- sample(seq(nrow(X1)), size = nrow(X1), replace = FALSE)
      sampleIdx2 <- sample(seq(nrow(X2)), size = nrow(X2), replace = FALSE)

      X1.b <- X1[sampleIdx1, ]
      X2.b <- X2[sampleIdx2, ]

      scca1 <- scca(X1.b, Y1, penalty = "LASSO")
      scca2 <- scca(X2.b, Y2, penalty = "LASSO")

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

      scca1 <- scca(X1.b, Y1, penalty = "LASSO")
      scca2 <- scca(X2.b, Y2, penalty = "LASSO")

      B1 <- outer(c(scca1$A), c(scca1$B))
      B2 <- outer(c(scca2$A), c(scca2$B))

      return(matrixDissimilarity(B1, B2))
      }, simplify = TRUE)
  }
  return(dissimilarities)
}

#' P-value of Similarity in Gene-drug Associations
#'
#' Computes p-value as number of null dissimilarities less than or equal to observed dissimilarity.
#'
#' @param dissimilarity observed dissimilarity
#' @param dissimilarities null distribution of dissimilarities
#' @return p-value
#' @examples
#' # simulate null distribution of dissimilarities
#' dissimilarities = runif(100, min = 0, max = 1)
#'
#' d = 0.10
#' p_value(d, dissimilarities)
#'
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

#' Select Significant Results from 'HierBipartite'
#'
#' Selects clusters from bipartite graph-based hierarchical clustering with p-value less than or equal to a p-value cutoff.
#'
#' @param results list of results from bipartite graph-based hierarchical clustering
#' @param p p-value cutoff
#' @return list of results from bipartite graph-based hierarchical clustering, but only with clusters with p-value at or below p-value cutoff
#' @examples
#' # sample bipartite graph-based hierarchical clustering of three groups
#' data(ctrp2)
#'
#' groups = ctrp2$groups
#' X = ctrp2$X
#' Y = ctrp2$Y
#'
#' groupNames = names(groups)
#' groupSmall = groups[groupNames[1:3]]
#'
#' \dontrun{
#' result = hierBipartite(X, Y, groupSmall)
#'
#' # set fictitious p-values, with one cluster with p-value less than the cutoff
#' # and the other not
#' result$nodePvals = list(0.03, 0.12)
#' getSignificantMergedGroups(result, p = 0.05)
#' }
#'
#' @export
getSignificantMergedGroups <- function(results, p = 0.05) {
  # filters bipartite hierarchical clustering merged groups by p-value
  # Input:
  #   results: list of results from bipartite hierarchical clustering
  #   p: p-value to filter merged groups by
  # Output:
  #   retLst; list of results from bipartite hierarchical clustering filtered by p-value
  if (!"nodePvals" %in% names(results)) {
    stop("p-value must be computed first!")
  } else {
    nodePvals <- results[["nodePvals"]]
    groupMerges <- results[["groupMerges"]]

    retLst <- list()
    nodePvalsFiltered <- list()
    groupMergesFiltered <- list()

    n <- length(nodePvals)
    index <- 1
    for (i in seq(n)) {
      if (!is.na(nodePvals[[i]]) && nodePvals[[i]] <= p) {
        nodePvalsFiltered[[toString(index)]] <- nodePvals[[i]]
        groupMergesFiltered[[toString(index)]] <- groupMerges[[i]]
      }
      index <- index + 1
    }
    retLst <- list("nodePvals" = nodePvalsFiltered, "hclustObj" = results[["hclustObj"]],
      "groupMerges" = groupMergesFiltered, "D" = results[["D"]])
    return(retLst)
  }
}
