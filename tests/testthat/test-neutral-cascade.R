test_that("neutralNetAssembly returns the documented structure", {
  set.seed(1)
  out <- neutralNetAssembly(20, 0.1, 0.2, 0.1, 50, BB = 3L)

  expect_named(out, c("S", "L", "B", "A"))
  expect_length(out$S, 50)
  expect_length(out$L, 50)
  expect_equal(dim(out$A), c(20, 20))
  expect_true(all(out$S >= 0 & out$S <= 20))
  expect_true(all(out$L >= 0))
})

test_that("neutralNetAssembly forbids self interactions", {
  set.seed(2)
  out <- neutralNetAssembly(30, 0.2, 0.3, 0.05, 100, BB = 3L)
  expect_false(any(diag(out$A)))
})

test_that("neutralNetAssembly forbids double interactions", {
  set.seed(2)
  out <- neutralNetAssembly(30, 0.2, 0.3, 0.05, 100, BB = 3L)
  expect_false(any(out$A & t(out$A)))
})

test_that("neutralNetAssembly with no migration and no new species keeps the initial basal species", {
  set.seed(1)
  out <- neutralNetAssembly(20, 0, 0.2, 0.1, 30, BB = 5L)
  expect_equal(out$S[1], 5)
  expect_true(all(out$S <= 5))
})

test_that("cascadeNetAssembly returns the documented structure", {
  set.seed(1)
  out <- cascadeNetAssembly(20, 0.1, 0.2, 0.1, 50, BB = 3L)

  expect_named(out, c("S", "L", "B", "A"))
  expect_length(out$S, 50)
  expect_equal(dim(out$A), c(20, 20))
})

test_that("cascadeNetAssembly interactions follow the cascade hierarchy", {
  set.seed(3)
  out <- cascadeNetAssembly(30, 0.2, 0.3, 0.05, 100, BB = 3L)
  A <- out$A
  expect_false(any(diag(A)))
  # only one triangle of the matrix can contain links
  expect_true(!any(A[upper.tri(A)]) || !any(A[lower.tri(A)]))
})

test_that("cascadeNetAssembly with no migration and no extinction keeps BB basal species", {
  set.seed(1)
  out <- cascadeNetAssembly(20, 0, 0.2, 0, 30, BB = 5L)
  expect_true(all(out$S == 5))
  expect_equal(sum(out$A), 0)
})

test_that("neutral and cascade models are reproducible with the same seed", {
  set.seed(5); a <- neutralNetAssembly(20, 0.1, 0.2, 0.1, 40, BB = 2L)
  set.seed(5); b <- neutralNetAssembly(20, 0.1, 0.2, 0.1, 40, BB = 2L)
  expect_identical(a, b)

  set.seed(5); a <- cascadeNetAssembly(20, 0.1, 0.2, 0.1, 40, BB = 2L)
  set.seed(5); b <- cascadeNetAssembly(20, 0.1, 0.2, 0.1, 40, BB = 2L)
  expect_identical(a, b)
})
