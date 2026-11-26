ed <- data.frame(from  = c("alga", "grazer", "alga",    "grazer2", "grazer",      "alga",         "grazer"),
                 to    = c("grazer", "fish", "grazer2", "fish",    "grazer2",     "worm",         "worm"),
                 layer = c("Trophic", "Trophic", "Trophic", "Trophic", "Competition", "Commensalism", "Amensalism"))

test_that("layers are translated with the sign convention", {
  b <- build_multilayer_assembly(ed, beta = c(pos = 1, neg = 2), obligate_on = "alga", normalize = FALSE)
  expect_true(b$metaW["alga", "grazer"])
  expect_false(b$metaW["grazer", "alga"])
  expect_equal(b$M["grazer", "grazer2"], 2)          # competition raises extinction, both ways
  expect_equal(b$M["grazer2", "grazer"], 2)
  expect_equal(b$M["worm", "alga"], -1)              # commensalism 0 -> +: alga helps worm
  expect_equal(b$M["alga", "worm"], 0)
  expect_equal(b$M["grazer", "worm"], 2)             # amensalism - -> 0: worm harms grazer
  expect_equal(b$M["worm", "grazer"], 0)
  expect_true(b$Dep["worm", "alga"])                 # obligate on the species that helps it
  expect_equal(sum(b$Dep), 1)
})

test_that("symmetric links recorded in both directions count once, and layer_beta scales", {
  ed2 <- rbind(ed, data.frame(from = "grazer2", to = "grazer", layer = "Competition"))
  b <- build_multilayer_assembly(ed2, beta = c(pos = 1, neg = 1), layer_beta = c(Competition = 0.5), normalize = FALSE)
  expect_equal(b$M["grazer", "grazer2"], 0.5)
})

test_that("unknown layers stop", {
  ed3 <- rbind(ed, data.frame(from = "a", to = "b", layer = "Parasitism"))
  expect_error(build_multilayer_assembly(ed3), "Parasitism")
})

test_that("igraph list input gives the same matrices as the edge list", {
  skip_if_not_installed("igraph")
  gl <- lapply(split(ed, ed$layer), function(x) igraph::graph_from_data_frame(x[, 1:2]))
  b1 <- build_multilayer_assembly(gl, beta = c(pos = 1, neg = 1))
  b2 <- build_multilayer_assembly(ed, beta = c(pos = 1, neg = 1))
  sp <- b2$species
  expect_equal(b1$M[sp, sp], b2$M)
  expect_equal(b1$metaW[sp, sp], b2$metaW)
})

test_that("with beta = 0 and no obligate rule it matches metaWebNetAssemblyCT in distribution", {
  set.seed(1)
  S <- 30
  A <- matrix(runif(S * S) < 0.1, S, S); diag(A) <- FALSE
  sp <- paste0("s", 1:S); dimnames(A) <- list(sp, sp)
  idx <- which(A, arr.ind = TRUE)
  el <- data.frame(from = sp[idx[, 1]], to = sp[idx[, 2]], layer = "Trophic")
  nsim <- 300; tf <- 50
  s_ml <- replicate(nsim, {
    r <- metaWebNetAssemblyML(el, m = .3, e = .1, se = .5, time = tf)
    sum(r$STime[colnames(A), tf - 1])
  })
  s_ct2 <- replicate(nsim, sum(metaWebNetAssemblyCT(A, rep(.3, S), rep(.1, S), rep(.5, S), tf)$STime[, tf - 1]))
  expect_gt(suppressWarnings(ks.test(s_ml, s_ct2)$p.value), 0.01)
})

test_that("mutualism lowers and competition raises the number of species", {
  set.seed(2)
  S <- 20; sp <- paste0("s", 1:S)
  basal <- data.frame(from = character(0), to = character(0), layer = character(0))
  pairs <- t(combn(sp, 2))[sample(choose(S, 2), 60), ]
  mk <- function(type) rbind(basal, data.frame(from = pairs[, 1], to = pairs[, 2], layer = type))
  run <- function(type, b) mean(replicate(100, {
    r <- metaWebNetAssemblyML(mk(type), beta = c(pos = b, neg = b), m = .1, e = .1, se = .5, time = 60)
    sum(r$STime[, 59])
  }))
  s0 <- run("Mutualism", 0)
  expect_gt(run("Mutualism", 0.5), s0)
  expect_lt(run("Competition", 0.5), s0)
})

test_that("obligate dependents go extinct when the foundation species is never present", {
  set.seed(3)
  ed4 <- data.frame(from = c("found", "found", "other"), to = c("dep", "dep2", "dep"),
                    layer = c("Commensalism", "Mutualism", "Trophic"))
  sp <- build_multilayer_assembly(ed4)$species
  mm <- ifelse(sp == "found", 0, 1)
  r <- metaWebNetAssemblyML(ed4, m = mm, e = 0.01, se = 1000, time = 50, obligate_on = "found")
  expect_true(mean(r$STime["dep2", 10:49]) < 0.05)
  r0 <- metaWebNetAssemblyML(ed4, m = mm, e = 0.01, se = 1000, time = 50)
  expect_true(mean(r0$STime["dep2", 10:49]) > 0.9)
})

test_that("a sign attribute per link overrides the layer signs", {
  ed5 <- data.frame(from = c("a", "b", "a", "c"), to = c("b", "a", "c", "a"),
                    layer = c("Mutualism", "Mutualism", "Amensalism", "Trophic"),
                    effect = c("+", "+", "-", NA))
  b <- build_multilayer_assembly(ed5, beta = c(pos = 1, neg = 2), normalize = FALSE)
  expect_equal(b$M["b", "a"], -1); expect_equal(b$M["a", "b"], -1)
  expect_equal(b$M["c", "a"], 2)                    # a harms c (effect of from on to)
  expect_equal(b$M["a", "c"], 0)
  b2 <- build_multilayer_assembly(ed5, beta = c(pos = 1, neg = 2), obligate_on = "a", obligate_layers = "Commensalism")
  expect_equal(sum(b2$Dep), 0)
  b3 <- build_multilayer_assembly(ed5, beta = c(pos = 1, neg = 2), obligate_on = "a")
  expect_true(b3$Dep["b", "a"])
})

test_that("normalize divides by the number of partners of each sign", {
  ed6 <- data.frame(from = c("a", "b", "c"), to = c("x", "x", "x"), layer = "Commensalism", effect = "+")
  b <- build_multilayer_assembly(ed6, beta = c(pos = 1, neg = 1))
  expect_equal(sum(b$M["x", ]), -1)
  expect_equal(b$M["x", "a"], -1/3)
})
