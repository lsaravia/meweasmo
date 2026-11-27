# Simulation of an assembly process from a multilayer meta-web

Continuous time colonization-extinction model of
[`metaWebNetAssemblyCT()`](https://lsaravia.github.io/meweasmo/reference/metaWebNetAssemblyCT.md)
in which the non-trophic layers change the extinction rates and,
optionally, some species depend obligately on others (e.g. on foundation
species). Losing a partner with a positive effect (a mutualist, a
facilitator) raises the extinction rate of a species; losing a partner
with a negative effect (a competitor) lowers it.

## Usage

``` r
metaWebNetAssemblyML(
  mg,
  layers = default_layer_signs(),
  beta = c(pos = 0, neg = 0),
  m,
  e,
  se,
  time,
  layer_beta = NULL,
  obligate_on = NULL,
  type_attr = "type",
  effect_attr = "effect",
  obligate_layers = NULL,
  normalize = TRUE
)
```

## Arguments

- mg:

  the multilayer network: an `mgraph` (the list of igraph layers
  returned by
  [`multiweb::readMultiplex()`](https://lsaravia.github.io/multiweb/reference/readMultiplex.html))
  or any list of igraph layers, where the layer id is the list name or,
  if the list has no names, the `type_attr` edge attribute; a single
  igraph with the layer id in the edge attribute `type_attr`; or a
  data.frame with columns `from`, `to`, `layer` and, optionally, the
  sign column `effect_attr`. Vertices need a `name` attribute.

- layers:

  data.frame with one row per layer: `layer` (the id), `sign_from`,
  `sign_to` and `rule` (`"trophic"` or `"modifier"`); see
  [`default_layer_signs()`](https://lsaravia.github.io/meweasmo/reference/default_layer_signs.md).
  Every layer of `mg` must be in it, and exactly one must have
  `rule == "trophic"`.

- beta:

  named numeric vector `c(pos = , neg = )`: strength of the positive and
  of the negative effects on the log extinction rate. With
  `normalize = TRUE`, `exp(beta)` is the largest factor by which the
  partners can change the extinction rate.

- m, e, se:

  migration, baseline extinction and secondary extinction rates, one
  value or one per species, in the order of the `species` element of the
  result (the order in which species first appear in `mg`).

- time:

  simulated time, the state is recorded at integer times 1 to time-1.

- layer_beta:

  optional named numeric vector, a multiplier of `beta` for each layer
  (default 1).

- obligate_on:

  character vector of species on which others can depend obligately
  (e.g. foundation species), or `NULL` (default) for no obligate
  dependence.

- type_attr:

  name of the edge attribute with the layer id.

- effect_attr:

  name of the edge attribute (or data.frame column) with the sign of the
  effect of `from` on `to`. If it is `NULL` or missing, the layer signs
  are used.

- obligate_layers:

  layers whose positive effects create obligate dependence; `NULL`
  (default) uses all modifier layers.

- normalize:

  if `TRUE` (default) the positive and negative effects on each species
  are divided by its number of positive and negative partners, so that
  the effect depends on the fraction of partners present. With `FALSE`
  the effects add up per partner present, and a species with many
  mutualists can become practically immortal.

## Value

the list of
[`metaWebNetAssemblyMLCT()`](https://lsaravia.github.io/meweasmo/reference/metaWebNetAssemblyMLCT.md)
(`S`, `L`, `STime`, `A`, `rejected`) with species names in `STime` and
`A`, plus `metaW`, `M`, `Dep` and `species` from
[`build_multilayer_assembly()`](https://lsaravia.github.io/meweasmo/reference/build_multilayer_assembly.md).

## Details

The extinction rate of species \\i\\ while species are present is
\$\$e_i(t) = e_i \exp(\beta\_- f^-\_i(t) - \beta\_+ f^+\_i(t))\$\$ where
\\f^+\_i(t)\\ and \\f^-\_i(t)\\ are the fractions of the partners of
\\i\\ with a positive and a negative effect on it that are present at
time \\t\\ (`normalize = TRUE`, the default). Thus \\\exp(\beta\_+)\\
and \\\exp(\beta\_-)\\ are the largest factors by which the partners can
lower or raise the extinction rate. Migration and the trophic rule (a
predator without prey goes extinct with rate *se*) are as in
[`metaWebNetAssemblyCT()`](https://lsaravia.github.io/meweasmo/reference/metaWebNetAssemblyCT.md).
With `obligate_on`, a dependent species goes extinct with rate *se* when
none of the species it depends on is present.

The parameters `beta` (and `layer_beta`) are translated by
[`build_multilayer_assembly()`](https://lsaravia.github.io/meweasmo/reference/build_multilayer_assembly.md)
into the matrix `M` of effects on the log extinction rate used by the
C++ core
[`metaWebNetAssemblyMLCT()`](https://lsaravia.github.io/meweasmo/reference/metaWebNetAssemblyMLCT.md).
Call
[`build_multilayer_assembly()`](https://lsaravia.github.io/meweasmo/reference/build_multilayer_assembly.md)
to check how a network was translated.

With `beta = c(pos = 0, neg = 0)` and `obligate_on = NULL` the model is
[`metaWebNetAssemblyCT()`](https://lsaravia.github.io/meweasmo/reference/metaWebNetAssemblyCT.md)
on the trophic layer.

## See also

[`build_multilayer_assembly()`](https://lsaravia.github.io/meweasmo/reference/build_multilayer_assembly.md),
[`default_layer_signs()`](https://lsaravia.github.io/meweasmo/reference/default_layer_signs.md),
[`vignette("multilayer-assembly")`](https://lsaravia.github.io/meweasmo/articles/multilayer-assembly.md)

## Examples

``` r
ed <- data.frame(from  = c("alga", "grazer", "alga",       "grazer2", "grazer"),
                 to    = c("grazer", "fish",  "grazer2",   "fish",    "grazer2"),
                 layer = c("Trophic", "Trophic", "Trophic", "Trophic", "Competition"))
r <- metaWebNetAssemblyML(ed, beta = c(pos = 1, neg = 1), m = 0.1, e = 0.05, se = 0.5, time = 100)
r$STime[, 99]
#>    alga  grazer grazer2    fish 
#>       1       1       0       0 
```
