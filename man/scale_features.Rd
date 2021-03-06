% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{scale_features}
\alias{scale_features}
\title{Only scales features (columns) with positive variance}
\usage{
scale_features(mat)
}
\arguments{
\item{mat}{an n x p matrix (gene expression or drug sensitivity)}
}
\value{
an n x p matrix with columns with positive variance scaled
}
\description{
Only scales features (columns) with positive variance
}
\keyword{internal}
