# meweasmo 0.3.0

## New features

* `metaWebNetAssemblyML()`: multilayer version of the continuous-time metaweb
  assembly model. Non-trophic layers change the extinction rates: a species
  present with a positive effect (mutualism, commensalism) lowers the
  extinction rate of its partner and one with a negative effect (competition,
  amensalism) raises it, `e_i(t) = e_i exp(beta_neg f-_i - beta_pos f+_i)`,
  with `f` the fraction of partners present (`normalize = FALSE` adds the
  effects per partner instead). Optional obligate dependence
  (`obligate_on`, `obligate_layers`): dependent species go extinct with rate
  `se` when none of the species they depend on (e.g. foundation species) is
  present. Takes an `mgraph` or list of igraph layers, an igraph with a layer
  attribute, or an edge data.frame; uses a per-link sign attribute
  (`effect_attr`) when present, otherwise the layer signs of
  `default_layer_signs()`. With `beta = 0` and no obligate rule it is
  `metaWebNetAssemblyCT()`.
* `build_multilayer_assembly()`: the matrices of the multilayer model, to check
  how the layers were translated.
* `metaWebNetAssemblyMLCT()`: the C++ core of the multilayer model
  (extinction rates simulated exactly by thinning).

# meweasmo 0.2.0

## New features

* `metaWebNetAssemblyCT()`: continuous-time version of the metaweb assembly
  model, simulated with the Gillespie algorithm (migration, local extinction
  and secondary extinction events with per-species rates).
* `generateGLVparmsFromAdjR()`: pure R implementation of
  `generateGLVparmsFromAdj()`.
* `generateGLVparmsFromAdj()` supports gamma-distributed interaction strengths
  (`rndType = 1`, with `predIntSd`) in addition to uniform ones.
* `metaWebNetAssemblyGLV()`: `ini` can now be a single value, recycled to all
  species (e.g. `0` for an empty initial community).
* Added a `testthat` test suite covering all exported functions.

## Bug fixes

* `generateGLVparmsFromAdjR()`: with `rndType = 1` the gamma scale was passed
  as `rate` to `rgamma()`, giving interaction strengths about 8 orders of
  magnitude too large. It now matches `generateGLVparmsFromAdj()`.
* `neutralNetAssembly()`: double interactions (A -> B and B -> A) were
  documented as forbidden but could be created; they are now prevented.
* `metaWebNetAssemblyGLV()`: a scalar `ini` (documented as allowed) shrank the
  population vector and read memory out of bounds. Arguments `ini`, `m` and `r`
  are now validated and give an error if their length does not match
  `nrow(metaW)`.
* `generateGLVparmsFromAdj()`: no longer modifies the input matrix in place
  when it has double storage.

## Documentation

* `metaWebNetAssembly()` and `metaWebNetAssemblyCT()`: documented the
  differences between the discrete and continuous time versions (probabilities
  vs. rates, synchronous vs. one event at a time, cascading secondary
  extinctions), cross-referenced them, and cite Saravia et al. (2022) as the
  main reference.
* `metaWebNetAssemblyGLV()`: rewrote the documentation. It now gives the model
  equation, explains the tau-leap integration, the meaning of the interaction
  matrix signs, that `m` is an immigration rate in individuals per unit of time
  (not a probability) and that `tau` is the step size (not a number of steps),
  describes each element of the returned list, and adds a reference and examples.
* `generateGLVparmsFromAdj()` and `generateGLVparmsFromAdjR()`: documented how
  interaction types, species types (mutualistic, basal, consumer),
  self-limitation, growth and migration rates are generated. Removed the
  statement that the Rcpp version "presumably has memory issues": the in-place
  modification of the input was fixed and both versions give the same results.
* roxygen2 now uses markdown (`Roxygen: list(markdown = TRUE)`).
* `calcPropInteractionsGLVadjMat()`: documented that proportions are computed
  over directed matrix entries, so mutualism and competition (two entries each)
  weigh twice as much as the other interaction types.
* `tests/testthat.R` referred to the old package name; fixed.

# meweasmo 0.1.8

* Initial versions: trophic theory of island biogeography based metaweb
  assembly (`metaWebNetAssembly()`), neutral and cascade assembly models,
  stochastic generalized Lotka-Volterra assembly (`metaWebNetAssemblyGLV()`)
  and helpers to generate GLV interaction matrices and parameters.
