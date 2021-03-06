#' Processed Cancer Cell Line Encyclopedia (CCLE) and Cancer Therapeutics Response Portal (CTRP2) Dataset
#'
#' Smaller test dataset version of the "CTRP2" carcinoma dataset in the paper. Specifically,
#' only the top 1,000 transcripts by correlation with drug sensitivity are included instead of 5,000.
#' Otherwise the dataset has been processed exactly as described in the paper. Note the expression dataset is
#' provided by CCLE and the drug sensitivity dataset is provided by CTRP2, and the pharmacogenomic datasets in the paper are
#' are referred to by the resource providing the sensitivity data. The cell lines are grouped by carcinoma subtype and primary site (e.g. lung NSC).
#'
#' @usage data(ctrp2)
#'
#' @format A list with elements of gene expression, drug sensitivities, and group membership.
#' \describe{
#'   \item{X}{n x p gene expression matrix}
#'   \item{Y}{n x q drug sensitivities matrix}
#'   \item{groups}{List of starting groups. Each group is represented by a vector of row indices for X, Y.}
#' }
#'
#' @examples
#' data(ctrp2)
#' X = ctrp2[["X"]]
#' Y = ctrp2[["Y"]]
#' groups = ctrp2[["groups"]]
#'
#' @references
#' Barretina, J., et al. (2012). The cancer cell line encyclopedia enables predictive modelling of anticancer drug sensitivity.
#' Nature, 483(7391), 603–607. (\href{https://pubmed.ncbi.nlm.nih.gov/22460905}{PubMed})
#'
#' Seashore-Ludlow, B., et al. (2015). Harnessing connectivity in a large-scale small-molecule sensitivity dataset.
#' Cancer discovery, 5(11), 1210-1223. (\href{https://pubmed.ncbi.nlm.nih.gov/26482930}{PubMed})
"ctrp2"
