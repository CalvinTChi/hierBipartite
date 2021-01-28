#' Given index from hclust merge matrix, return X, Y row indices
#' corresponding to groups specified by index
#'
#' @param idx index from hclust merge
#' @param groups list of starting group membership (input to hierBipartite())
#' @param groupMerges list of merged groups thus far (groupMerge[[3]] is a
#' vector of group names from cluster created at third merge)
#' @keywords internal
#' @return vector of row indices corresponding to X, Y
getMergeGroupRows <- function(idx, groups, groupMerges) {
  # Given index from hclust merge matrix, return X, Y row indices corresponding
  # to groups specified by index
  # Input:
  #   idx: index from hclust merge
  #   groups: list of starting group membership (input to hierBipartite())
  #   groupMerges: list of merged groups thus far (groupMerge[[3]] is a vector
  #                of group names from cluster created at third merge)
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
#' @keywords internal
#' @return vector of group names after merge of group indicated by idx1 and idx2
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

#' Only scales features (columns) with positive variance
#'
#' @param mat an n x p matrix (gene expression or drug sensitivity)
#' @keywords internal
#' @return an n x p matrix with columns with positive variance scaled
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

