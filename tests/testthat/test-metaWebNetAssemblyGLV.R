test_that("metaWebNetAssemblyGLV returns the documented structure", {
  A <- matrix(c(-0.001, -0.08, 0.01, 0), nrow = 2, byrow = TRUE)
  set.seed(123)
  out <- metaWebNetAssemblyGLV(A, c(0.1, 0.2), c(1.5, -0.25), c(0, 0), 20, 0.1)

  expect_named(out, c("S", "L", "STime", "A"))
  expect_length(out$S, 20)
  expect_length(out$L, 20)
  expect_equal(dim(out$STime), c(2, 20))
  expect_equal(dim(out$A), c(2, 2))
  expect_true(all(out$STime >= 0))
})

test_that("no migration and no initial individuals gives an empty system", {
  A <- matrix(c(-0.001, -0.08, 0.01, 0), nrow = 2, byrow = TRUE)
  set.seed(1)
  out <- metaWebNetAssemblyGLV(A, c(0, 0), c(1.5, -0.25), c(0, 0), 10, 0.1)
  expect_true(all(out$S == 0))
  expect_true(all(out$L == 0))
  expect_true(all(out$STime == 0))
})

test_that("populations stay constant without interactions, growth or migration", {
  A <- matrix(0, 3, 3)
  set.seed(1)
  out <- metaWebNetAssemblyGLV(A, rep(0, 3), rep(0, 3), c(5, 6, 7), 10, 0.1)
  expect_true(all(out$S == 3))
  expect_equal(out$STime[, 1], c(5, 6, 7))
  expect_true(all(out$STime == c(5, 6, 7)))   # recycled down the columns
})

test_that("initial conditions are recorded at time 1", {
  A <- matrix(0, 2, 2)
  set.seed(1)
  out <- metaWebNetAssemblyGLV(A, c(0, 0), c(0, 0), c(10, 0), 5, 0.1)
  expect_equal(out$STime[, 1], c(10, 0))
  expect_equal(out$S[1], 1)
})

test_that("the predator-prey link is directed prey -> predator in A", {
  # metaW[1, 2] < 0 and metaW[2, 1] > 0: sp 2 preys on sp 1
  A <- matrix(c(0, -1e-6, 1e-6, 0), nrow = 2, byrow = TRUE)
  set.seed(1)
  out <- metaWebNetAssemblyGLV(A, c(0, 0), c(0, 0), c(1000, 1000), 3, 0.1)
  expect_true(all(out$S == 2))
  expect_true(out$A[1, 2])
  expect_false(out$A[2, 1])
  expect_equal(out$L[3], 1)
})

test_that("metaWebNetAssemblyGLV is reproducible with the same seed", {
  A <- matrix(c(-0.001, -0.08, 0.01, 0), nrow = 2, byrow = TRUE)
  set.seed(9); a <- metaWebNetAssemblyGLV(A, c(0.1, 0.2), c(1.5, -0.25), c(10, 10), 20, 0.1)
  set.seed(9); b <- metaWebNetAssemblyGLV(A, c(0.1, 0.2), c(1.5, -0.25), c(10, 10), 20, 0.1)
  expect_identical(a, b)
})

test_that("a single initial value is recycled to all species", {
  A <- matrix(0, 3, 3)
  set.seed(1)
  out <- metaWebNetAssemblyGLV(A, rep(0, 3), rep(0, 3), 4, 5, 0.1)
  expect_equal(out$STime[, 1], c(4, 4, 4))
  expect_true(all(out$S == 3))

  out0 <- metaWebNetAssemblyGLV(A, rep(0, 3), rep(0, 3), 0, 5, 0.1)
  expect_true(all(out0$STime == 0))
  expect_true(all(out0$S == 0))
})

test_that("wrong argument lengths raise an error", {
  A <- matrix(0, 3, 3)
  expect_error(metaWebNetAssemblyGLV(A, rep(0, 3), rep(0, 3), c(1, 2), 5, 0.1), "ini")
  expect_error(metaWebNetAssemblyGLV(A, rep(0, 2), rep(0, 3), 1, 5, 0.1), "length")
  expect_error(metaWebNetAssemblyGLV(A, rep(0, 3), rep(0, 4), 1, 5, 0.1), "length")
})
