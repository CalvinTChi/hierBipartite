test_that("parallel_cores", {
  # null_distri
  X1 = matrix(rnorm(1000), nrow = 100)
  Y1 = matrix(rnorm(1000), nrow = 100)
  X2 = matrix(rnorm(1000), nrow = 100)
  Y2 = matrix(rnorm(1000), nrow = 100)

  d1 = null_distri(X1, Y1, X2, Y2, n.perm = 10, parallel = TRUE,
                   maxCores = 2)
  d2 = null_distri(X1, Y1, X2, Y2, n.perm = 10, parallel = FALSE)
  expect_equal(length(d1), length(d2))

  # constructBipartiteGraph
  B1 = constructBipartiteGraph(X1, Y1, n_subsample = 10, subsampling_ratio = 0.9,
                               parallel = TRUE, maxCores = 2)
  B2 = constructBipartiteGraph(X1, Y1, n_subsample = 10, subsampling_ratio = 0.9,
                               parallel = FALSE)
  expect_equal(dim(B1), dim(B2))
})

test_that("getSignificantMergedGroups", {
  data(ctrp2)

  groups = ctrp2$groups
  X = ctrp2$X
  Y = ctrp2$Y

  groupNames = names(groups)
  groupSmall = groups[groupNames[1:3]]

  result = hierBipartite(X, Y, groupSmall)
  expect_error(getSignificantMergedGroups(result),
               "p-value must be computed first!")
})
