
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hierBipartite

[![CRAN](https://www.r-pkg.org/badges/version/hierBipartite)](https://www.r-pkg.org/pkg/hierBipartite)
[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/hierBipartite)](https://CRAN.R-project.org/package=hierBipartite)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![MIT
license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](https://opensource.org/licenses/MIT)

> Bipartite Graph-based Hierarchical Clustering

**Author:** [Calvin Chi](https://calvintchi.github.io/)

-----

## About

Bipartite graph-based hierarchical clustering performs hierarchical
clustering of groups of samples based on association patterns between
two sets of variables. It is developed for pharmacogenomic datasets and
datasets sharing the same data structure. In the context of
pharmacogenomic datasets, the samples are cell lines, and the two sets
of variables are typically expression levels and drug sensitivity
values. For this method, sparse canonical correlation analysis from Lee,
W., Lee, D., Lee, Y. and Pawitan, Y. (2011) <doi:10.2202/1544-6115.1638>
is first applied to extract association patterns for each group of
samples. Then, a nuclear norm-based dissimilarity measure is used to
construct a dissimilarity matrix between groups based on the extracted
associations. Finally, hierarchical clustering is applied.

-----

## Installation

For standard use, we recommend installing the package from
[CRAN](https://CRAN.R-project.org/package=hierBipartite) via

``` r
install.packages("hierBipartite")
```

-----

## Example

Minimially sufficient example of using `hierBipartite`. Please refer to
the
[vignette](https://cran.r-project.org/web/packages/hierBipartite/vignettes/hierBipartite.html)
for detailed usage.

``` r
library(hierBipartite)
data(ctrp2)

groups = ctrp2$groups
X = ctrp2$X
Y = ctrp2$Y

result = hierBipartite(X, Y, groups)
```
