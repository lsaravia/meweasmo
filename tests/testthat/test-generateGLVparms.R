make_adj <- function(n = 40, seed = 3) {
  set.seed(seed)
  generateRandomGLVadjMat(n, 0.15, c(0.2, 0.3, 0.3, 0.1, 0.1))
}

test_that("generateGLVparmsFromAdj returns the documented structure", {
  A <- make_adj()
  set.seed(1)
  p <- generateGLVparmsFromAdj(A, 0.1, 0.01, migrAvg = 0.5)

  expect_named(p, c("interM", "r", "m"))
  expect_equal(dim(p$interM), dim(A))
  expect_length(p$r, nrow(A))
  expect_length(p$m, nrow(A))
})

test_that("self-limitation is negative and migration is within [0, 2*migrAvg]", {
  A <- make_adj()
  set.seed(1)
  p <- generateGLVparmsFromAdj(A, 0.1, 0.01, migrAvg = 0.5)

  expect_true(all(diag(p$interM) <= 0))
  expect_true(all(p$m >= 0 & p$m <= 1))
})

test_that("migrAvg = 0 gives no migration", {
  A <- make_adj()
  set.seed(1)
  p <- generateGLVparmsFromAdj(A, 0.1, 0.01, migrAvg = 0)
  expect_true(all(p$m == 0))
})

test_that("off-diagonal signs and zeros of the adjacency matrix are preserved", {
  A <- make_adj()
  set.seed(1)
  p <- generateGLVparmsFromAdj(A, 0.1, 0.01)
  off <- row(A) != col(A)
  expect_identical(p$interM[off] == 0, A[off] == 0)
  expect_true(all(sign(p$interM[off]) == sign(A[off])))
})

test_that("predator gain equals efficiency times prey loss", {
  set.seed(4)
  A <- generateRandomGLVadjMat(40, 0.2, c(0, 1, 0, 0, 0))   # only predator-prey
  p <- generateGLVparmsFromAdj(A, 0.25, 0.01)
  off <- row(A) != col(A)
  prey_loss <- p$interM[off & A < 0]
  pred_gain <- t(p$interM)[off & A < 0]
  expect_equal(pred_gain, -0.25 * prey_loss)
})

test_that("uniform interactions have the requested mean", {
  A <- make_adj(150, 5)
  set.seed(1)
  p <- generateGLVparmsFromAdj(A, 0.1, 0.01, rndType = 0)
  x <- p$interM; diag(x) <- 0
  expect_equal(mean(abs(x[x != 0])), 0.01, tolerance = 0.2)
})

test_that("gamma interactions have the requested mean and sd", {
  A <- make_adj(150, 5)
  set.seed(1)
  p <- generateGLVparmsFromAdj(A, 0.1, 0.01, predIntSd = 0.001, rndType = 1)
  x <- p$interM; diag(x) <- 0
  x <- abs(x[x != 0])
  expect_equal(mean(x), 0.01, tolerance = 0.1)
  expect_equal(sd(x), 0.001, tolerance = 0.3)
})

test_that("preserveInt = 1 keeps the original interaction values", {
  A <- make_adj()
  A[A != 0] <- A[A != 0] * 0.005
  set.seed(1)
  p <- generateGLVparmsFromAdj(A, 0.1, 0.01, preserveInt = 1)
  off <- row(A) != col(A)
  expect_equal(p$interM[off], A[off])
})

test_that("the R and Rcpp implementations give the same results", {
  A <- make_adj()
  for (pres in 0:2) {
    set.seed(7); a <- generateGLVparmsFromAdj(A, 0.1, 0.01, preserveInt = pres, rndType = 0)
    set.seed(7); b <- generateGLVparmsFromAdjR(A, 0.1, 0.01, preserveInt = pres, rndType = 0)
    expect_equal(a, b, info = paste("preserveInt =", pres))
  }
})

test_that("the R implementation matches Rcpp with gamma-distributed interactions", {
  A <- make_adj()
  set.seed(7); a <- generateGLVparmsFromAdj(A, 0.1, 0.01, predIntSd = 0.001, rndType = 1)
  set.seed(7); b <- generateGLVparmsFromAdjR(A, 0.1, 0.01, predIntSd = 0.001, rndType = 1)
  expect_equal(a, b)
})

test_that("the R implementation gives gamma interactions with the requested mean and sd", {
  A <- make_adj(150, 5)
  set.seed(1)
  p <- generateGLVparmsFromAdjR(A, 0.1, 0.01, predIntSd = 0.001, rndType = 1)
  x <- p$interM; diag(x) <- 0
  x <- abs(x[x != 0])
  expect_equal(mean(x), 0.01, tolerance = 0.1)
  expect_equal(sd(x), 0.001, tolerance = 0.3)
})

test_that("generateGLVparmsFromAdj does not modify its input", {
  A <- make_adj() + 0   # double storage: Rcpp shares memory instead of copying
  A0 <- A + 0           # independent copy
  p <- generateGLVparmsFromAdj(A, 0.1, 0.01)
  expect_identical(A, A0)
})
