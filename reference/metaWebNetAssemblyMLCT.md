# Simulation of an assembly process from a multilayer meta-web: C++ core

Low level simulation used by
[`metaWebNetAssemblyML()`](https://lsaravia.github.io/meweasmo/reference/metaWebNetAssemblyML.md).
Most users should call
[`metaWebNetAssemblyML()`](https://lsaravia.github.io/meweasmo/reference/metaWebNetAssemblyML.md),
which takes the network as layers and the effects as the parameters
`beta`; this function takes instead the three matrices built from them
by
[`build_multilayer_assembly()`](https://lsaravia.github.io/meweasmo/reference/build_multilayer_assembly.md).
Use it directly to run many simulations of the same network without
rebuilding the matrices.

## Usage

``` r
metaWebNetAssemblyMLCT(metaW, M, Dep, m, e, se, time)
```

## Arguments

- metaW:

  Logical trophic adjacency matrix of the meta-web, `metaW[j, i]` is
  `TRUE` if *j* is a prey of *i*.

- M:

  Numeric matrix, `M[i, j]` is the effect of the presence of *j* on the
  log extinction rate of *i* (the diagonal is ignored).

- Dep:

  Logical matrix, `Dep[i, j]` is `TRUE` if *i* depends obligately on
  *j*.

- m:

  A numeric vector of species' migration rates from the meta-web

- e:

  A numeric vector of species' baseline extinction rates

- se:

  A numeric vector of species' secondary extinction rates (no prey, or
  no obligate partner)

- time:

  Simulated time, the state is recorded at integer times 1 to time-1

## Value

A list as
[`metaWebNetAssemblyCT()`](https://lsaravia.github.io/meweasmo/reference/metaWebNetAssemblyCT.md):
`S` and `L` (number of species with trophic links and of trophic links
at each time), `STime` (species presence, species x time) and `A`
(trophic adjacency matrix of the final local web), plus `rejected`, the
number of rejected extinction proposals.

## Details

It is the continuous time model of
[`metaWebNetAssemblyCT()`](https://lsaravia.github.io/meweasmo/reference/metaWebNetAssemblyCT.md)
with two additions:

- **Non-trophic modifiers of extinction**: the extinction rate of
  species *i* depends on which species are present, \\e_i(t) = e_i
  \exp(\sum_j M\_{ij} x_j(t))\\, where \\x_j(t)\\ is 1 if species *j* is
  present. `M` is in units of log extinction rate:
  [`build_multilayer_assembly()`](https://lsaravia.github.io/meweasmo/reference/build_multilayer_assembly.md)
  puts \\-\beta\_+\\ (divided by the number of positive partners if
  `normalize = TRUE`) where *j* has a positive effect on *i*, and
  \\+\beta\_-\\ where it has a negative one. So a negative \\M\_{ij}\\
  (a mutualist, a facilitator) lowers the extinction rate of *i* while
  *j* is present and a positive one (a competitor) raises it. With
  `M = 0` the model is
  [`metaWebNetAssemblyCT()`](https://lsaravia.github.io/meweasmo/reference/metaWebNetAssemblyCT.md).

- **Obligate dependence**: if species *i* depends on at least one
  species (`Dep[i, ]` has a `TRUE`) and none of them is present, *i*
  goes extinct with the secondary extinction rate *se*, as a predator
  without prey does.

The state dependent extinction rate is simulated exactly by thinning:
extinction events are proposed with the upper bound \\e_i \exp(\sum_j
\max(M\_{ij},0))\\ and accepted with probability \\e_i(t)\\ divided by
that bound. Large positive values in `M` make the bound large and the
simulation slower (more rejected proposals).

## See also

[`metaWebNetAssemblyML()`](https://lsaravia.github.io/meweasmo/reference/metaWebNetAssemblyML.md),
[`build_multilayer_assembly()`](https://lsaravia.github.io/meweasmo/reference/build_multilayer_assembly.md),
[`metaWebNetAssemblyCT()`](https://lsaravia.github.io/meweasmo/reference/metaWebNetAssemblyCT.md)
