test_that("generateRandomGLVadjMat returns a valid square GLV adjacency matrix", {
  set.seed(2)
  n <- 100
  M <- generateRandomGLVadjMat(n, 0.1, c(0.2, 0.3, 0.3, 0.1, 0.1))

  expect_true(is.matrix(M))
  expect_equal(dim(M), c(n, n))
  expect_true(all(M %in% c(-1, 0, 1)))
  expect_true(all(diag(M) == 0))
})

test_that("generateRandomGLVadjMat respects connectivity approximately", {
  set.seed(2)
  n <- 100
  M <- generateRandomGLVadjMat(n, 0.1, c(0.2, 0.3, 0.3, 0.1, 0.1))
  expect_equal(sum(M != 0) / n^2, 0.1, tolerance = 0.3)   # within 30% of C
})

test_that("generateRandomGLVadjMat respects the proportion of interactions approximately", {
  set.seed(2)
  n <- 100
  prop <- c(0.2, 0.3, 0.3, 0.1, 0.1)
  M <- generateRandomGLVadjMat(n, 0.1, prop)
  p <- calcPropInteractionsGLVadjMat(M, rep(1, n))
  expect_equal(p, prop, tolerance = 0.1, scale = 1)
})

test_that("predator-prey pairs have opposite signs", {
  set.seed(3)
  M <- generateRandomGLVadjMat(60, 0.2, c(0, 1, 0, 0, 0))   # only predator-prey
  nz <- M != 0
  expect_true(any(nz))
  expect_true(all(M[nz] == -t(M)[nz]))
})

test_that("generateRandomGLVadjMat is reproducible with the same seed", {
  set.seed(10); A <- generateRandomGLVadjMat(30, 0.1, rep(0.2, 5))
  set.seed(10); B <- generateRandomGLVadjMat(30, 0.1, rep(0.2, 5))
  expect_identical(A, B)
})
