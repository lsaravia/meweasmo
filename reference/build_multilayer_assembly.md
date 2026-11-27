# Matrices of the multilayer assembly model

Translates a multilayer network and the parameters `beta` into the three
matrices used by
[`metaWebNetAssemblyMLCT()`](https://lsaravia.github.io/meweasmo/reference/metaWebNetAssemblyMLCT.md):
the trophic meta-web `metaW`, the matrix `M` of effects on the log
extinction rate, and the obligate dependence matrix `Dep`. It is called
by
[`metaWebNetAssemblyML()`](https://lsaravia.github.io/meweasmo/reference/metaWebNetAssemblyML.md);
use it directly to check how the layers were translated.

## Usage

``` r
build_multilayer_assembly(
  mg,
  layers = default_layer_signs(),
  beta = c(pos = 0, neg = 0),
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

a list with `metaW` (logical, `metaW[j, i]` is `TRUE` if `j` is a prey
of `i`), `M` (numeric), `Dep` (logical), all with the species as
dimnames, and `species`.

## Details

**Signs of the links.** Each modifier link gives the sign (+1, 0 or -1)
of the effect that each of the two species receives from the other. The
signs are taken from one of two sources:

- a sign attribute per link (`effect_attr`, values `"+"`, `"-"` or
  numbers): the link `a -> b` is the effect of `a` on `b` only; an
  effect of `b` on `a` must be another link `b -> a`. A layer must have
  the attribute in all its links or in none.

- otherwise, the signs of the layer in `layers`: in a link `a -> b`, `a`
  receives `sign_from` and `b` receives `sign_to`.

**From signs to `M`.** `M[i, j]` is the effect of the presence of `j` on
the log extinction rate of `i`. A positive effect of `j` on `i` adds
`-beta["pos"]` to `M[i, j]` and a negative one adds `+beta["neg"]`, both
multiplied by `layer_beta[layer]` if given. Pairs linked in several
layers add their effects. In symmetric layers (`sign_from == sign_to`)
each pair counts once, however many directions are recorded. With
`normalize = TRUE` the negative and positive entries of each row are
then divided by the number of partners of that sign.

**Obligate dependence.** With `obligate_on`, `Dep[i, j]` is `TRUE` if
`j` is in `obligate_on`, `i` is not, and `j` has a positive effect on
`i` in one of `obligate_layers`.

Self-loops are kept in the trophic layer (cannibals are their own prey)
and dropped in the modifier layers.

## See also

[`metaWebNetAssemblyML()`](https://lsaravia.github.io/meweasmo/reference/metaWebNetAssemblyML.md),
[`default_layer_signs()`](https://lsaravia.github.io/meweasmo/reference/default_layer_signs.md)
