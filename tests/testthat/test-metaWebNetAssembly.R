make_metaweb <- function(n = 10, p = 0.3, seed = 1) {
  set.seed(seed)
  mw <- matrix(0, n, n)
  mw[upper.tri(mw)] <- rbinom(n * (n - 1) / 2, 1, p)
  mw
}

models <- list(discrete = metaWebNetAssembly, continuous = metaWebNetAssemblyCT)

for (nm in names(models)) {
  model <- models[[nm]]
  mw <- make_metaweb()
  n <- nrow(mw)

  test_that(paste(nm, "model returns the documented structure"), {
    set.seed(1)
    out <- model(mw, rep(0.2, n), rep(0.1, n), rep(0.1, n), 30)

    expect_named(out, c("S", "L", "STime", "A"))
    expect_length(out$S, 30)
    expect_length(out$L, 30)
    expect_equal(dim(out$STime), c(n, 30))
    expect_equal(dim(out$A), c(n, n))
    expect_true(all(out$S >= 0 & out$S <= n))
  })

  test_that(paste(nm, "model: no migration means no species"), {
    set.seed(1)
    out <- model(mw, rep(0, n), rep(0.1, n), rep(0.1, n), 10)
    expect_true(all(out$S == 0))
    expect_true(all(out$L == 0))
    expect_false(any(out$A))
  })

  test_that(paste(nm, "model: full migration and no extinction recovers the metaweb"), {
    set.seed(1)
    out <- model(mw, rep(1, n), rep(0, n), rep(0, n), 10)
    expect_equal(sum(out$A), sum(mw))
    expect_equal(out$L[10], sum(mw))
    expect_false(any(out$A & mw == 0))   # no links outside the metaweb
  })

  test_that(paste(nm, "model: local web is a subset of the metaweb"), {
    set.seed(2)
    out <- model(mw, rep(0.3, n), rep(0.2, n), rep(0.2, n), 50)
    expect_false(any(out$A & mw == 0))
    expect_true(all(out$L >= 0))
  })

  test_that(paste(nm, "model is reproducible with the same seed"), {
    set.seed(8); a <- model(mw, rep(0.3, n), rep(0.2, n), rep(0.2, n), 20)
    set.seed(8); b <- model(mw, rep(0.3, n), rep(0.2, n), rep(0.2, n), 20)
    expect_identical(a, b)
  })
}
