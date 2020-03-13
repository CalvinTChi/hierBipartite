#' Preprocessed Cancer Cell Line Encyclopedia (CCLE) Dataset
#'
#' The dataset from CCLE (2012) that has been preprocessed for bipartite hierarhical clustering. 
#' The preprocessing steps removed samples and features based on missingness, filtered out
#' genes with expression variance less than 3, and removed cell lines from tissues with less than
#' 15 cell lines.
#'
#' @format A list with elements of gene expression, drug sensitivities, and group membership.
#' \describe{
#'   \item{X}{n x p gene expression matrix}
#'   \item{Y}{n x q drug sensitivities matrix}
#'   \item{groups}{List of starting groups. Each group is represented by a vector of row indices for X, Y.}
#' }
#'
#' @examples
#' data(CCLE)
#' X = CCLE[["X"]]
#' Y = CCLE[["Y"]]
#' groups = CCLE[["groups"]]
#'
#' @references Barretina, J., et al. (2012). The cancer cell line encyclopedia enables predictive modelling of anticancer drug sensitivity. Nature, 483(7391), 603–607.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/22460905}{PubMed})
"CCLE"

#' Simulation Dataset
#'
#' Dataset simulated as described in the paper.
#'
#' @format A list with the elements of simulated gene expression, simulated drug sensitivities, and group membership.
#' \describe{
#'   \item{X}{n x p gene expression matrix}
#'   \item{Y}{n x q drug sensitivities matrix}
#'   \item{groups}{List of starting groups. Each group is represented by a vector of row indices for X, Y.}
#' }
#'
#' @examples
#' data(simulation)
#' X = simulation[["X"]]
#' Y = simulation[["Y"]]
#' groups = simulation[["groups"]]
"simulation"