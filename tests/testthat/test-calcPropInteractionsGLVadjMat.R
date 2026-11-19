# 10 species: one pair of each of the 5 interaction types
adj_all_types <- function() {
  M <- matrix(0, 10, 10)
  M[1, 2] <- M[2, 1] <- 1      # mutualism    (+,+)
  M[3, 4] <- -1; M[4, 3] <- 1  # predator-prey (-,+)
  M[5, 6] <- M[6, 5] <- -1     # competition  (-,-)
  M[7, 8] <- 1                 # commensalism (+,0)
  M[9, 10] <- -1               # amensalism   (-,0)
  M
}

# The function counts directed matrix entries: symmetric interactions
# (mutualism, competition) contribute 2, the others contribute 1.
test_that("interaction types are counted as directed entries", {
  p <- calcPropInteractionsGLVadjMat(adj_all_types(), rep(1, 10))
  expect_length(p, 5)
  expect_equal(p, c(2, 1, 2, 1, 1) / 7)
})

test_that("proportions sum to one", {
  set.seed(1)
  M <- generateRandomGLVadjMat(50, 0.15, c(0.2, 0.3, 0.3, 0.1, 0.1))
  p <- calcPropInteractionsGLVadjMat(M, rep(1, 50))
  expect_equal(sum(p), 1)
})

test_that("absent species are not counted", {
  spc <- rep(1, 10)
  spc[1:2] <- 0   # removes the mutualistic pair
  p <- calcPropInteractionsGLVadjMat(adj_all_types(), spc)
  expect_equal(p, c(0, 1, 2, 1, 1) / 5)
})
