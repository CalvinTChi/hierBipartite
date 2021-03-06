% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HierBipartite.R
\name{hierBipartite}
\alias{hierBipartite}
\title{Bipartite Graph-based Hierarchical Clustering}
\usage{
hierBipartite(
  X,
  Y,
  groups,
  link = "ward.D2",
  n_subsample = 1,
  subsampling_ratio = 1,
  p.value = FALSE,
  n_perm = 100,
  parallel = FALSE,
  maxCores = 7,
  p_cutoff = 0.1
)
}
\arguments{
\item{X}{an n x p matrix of variable set 1 (e.g. gene expression)}

\item{Y}{an n x q matrix of variable set 2 (e.g. drug sensitivity)}

\item{groups}{a list of starting group membership (e.g. list("1" = c(1,2,3), "2" = c(4,5,6)) means group 1 has samples 1, 2, 3, and group 2 has samples 4, 5, 6.}

\item{link}{string indicating link function as input to hclust(). One of "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid".}

\item{n_subsample}{number of subsampling to generate matrix B (see paper)}

\item{subsampling_ratio}{fraction of samples to sample for subsampling to generate matrix B (see paper)}

\item{p.value}{boolean for whether to generate p-values for each merge}

\item{n_perm}{number of permutations for generating p-values. Ignored if p.value = FALSE}

\item{parallel}{boolean for whether to parallelize subsampling and p-value generation step}

\item{maxCores}{maximum number of cores to use (only applicable when parallel = TRUE)}

\item{p_cutoff}{p-value cutoff that determines whether merge is significant. If p-value > p_cutoff, p-values will not be calculated for future merges involving current group. Ignored if p.value = FALSE.}
}
\value{
list of results from bipartite graph-based hierarchical clustering, containing up to
\itemize{
  \item hclustObj: hclust object
  \item groupMerges: list of clusters after each merge, in order of merge.
  Each cluster is indicated by a vector of cell line groups
  \item nodePvals: list of p-value of each new merge,
  in order of merge. Only available if p.value = TRUE
  \item D: dissimilarity matrix
}
}
\description{
Main bipartite graph-based hierarchial clustering algorithm. Visit \href{https://calvintchi.github.io/html/hierBipartite}{here} for vignette on using the
hierBipartite package.
}
\examples{
# Get a small subset of the test dataset
data(ctrp2)

groups = ctrp2$groups
X = ctrp2$X
Y = ctrp2$Y

groupNames = names(groups)
groupSmall = groups[groupNames[1:3]]

\dontrun{
# Basic call of hierBipartite() on small test dataset
result0 = hierBipartite(X, Y, groupSmall)

# Calling hierBipartite() with subsampling
result1 = hierBipartite(X, Y, groupSmall, n_subsample = 100, subsampling_ratio = 0.90)

# Calling hierBipartite() with p-value generation
result2 = hierBipartite(X, Y, groupSmall, n_perm = 100, p.value = TRUE, p_cutoff = 0.10)

# Calling hierBipartite() with both subsampling and p-value generation (expensive)
result3 = hierBipartite(X, Y, groupSmall, n_subsample = 100, subsampling_ratio = 0.90,
                        n_perm = 100, p.value = TRUE, p_cutoff = 0.10)
}

}
