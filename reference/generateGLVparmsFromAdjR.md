# R implementation of generateGLVparmsFromAdj

Generates random parameters of a generalized Lotka-Volterra model from
an adjacency matrix. It is equivalent to
[`generateGLVparmsFromAdj()`](https://lsaravia.github.io/meweasmo/reference/generateGLVparmsFromAdj.md):
it has the same arguments and, for the same random seed, it gives the
same results (this is checked by the package tests). It is written in R,
so it is slower but easier to read and modify, and it can be used as a
reference implementation. See
[`generateGLVparmsFromAdj()`](https://lsaravia.github.io/meweasmo/reference/generateGLVparmsFromAdj.md)
for a detailed description of how the parameters are generated.

## Usage

``` r
generateGLVparmsFromAdjR(
  adjM,
  ef,
  predIntAvg = 0.01,
  selfLimAvg = c(0.01, 0.01, 0.01),
  migrAvg = 0,
  preserveInt = 0,
  predIntSd = 0.01,
  rndType = 0
)
```

## Arguments

- adjM:

  Numeric square matrix with the structure of the interactions. Only the
  signs and the zeros are used when `preserveInt = 0`; with
  `preserveInt` 1 or 2 the values of the non-zero off-diagonal entries
  are used as well.

- ef:

  Efficiency of predator-prey interactions: effect on the predator /
  effect on the prey.

- predIntAvg:

  Average strength of the interactions (predator-prey, competition and
  mutualism) when `preserveInt = 0`.

- selfLimAvg:

  Numeric vector of length 3 with the average self-limitation of 1 =
  mutualistic, 2 = basal and 3 = consumer species (diagonal entries of
  the interaction matrix).

- migrAvg:

  Average migration rate, the migration of each species is uniform in
  \[0, 2 \* `migrAvg`\].

- preserveInt:

  How to define the strength of the interactions:

  - 0: random values defined by `rndType` with mean `predIntAvg`.

  - 1: the values of `adjM` are preserved (off-diagonal entries only).

  - 2: random uniform values with mean given by the value in `adjM` for
    each interaction, that is, uniform in \[0, 2 \* abs(`adjM[i, j]`)\].

- predIntSd:

  Standard deviation of the strengths when `preserveInt = 0` and
  `rndType = 1` (gamma distribution). For the uniform distribution the
  sd is fixed at `predIntAvg * 2 / sqrt(12)`.

- rndType:

  Distribution of the strengths when `preserveInt = 0`: 0 = uniform, 1 =
  gamma. Other values of `preserveInt` always use the uniform
  distribution.

## Value

A list with the interaction matrix `interM` (with the self-limitation in
the diagonal), the intrinsic growth rates `r` and the migration rates
`m`. They can be passed to
[`metaWebNetAssemblyGLV()`](https://lsaravia.github.io/meweasmo/reference/metaWebNetAssemblyGLV.md)
as `metaW`, `r` and `m`.

## See also

[`generateGLVparmsFromAdj()`](https://lsaravia.github.io/meweasmo/reference/generateGLVparmsFromAdj.md)
for the Rcpp implementation,
[`generateRandomGLVadjMat()`](https://lsaravia.github.io/meweasmo/reference/generateRandomGLVadjMat.md)
to generate `adjM`, and
[`metaWebNetAssemblyGLV()`](https://lsaravia.github.io/meweasmo/reference/metaWebNetAssemblyGLV.md).

## Examples

``` r
set.seed(1)
adj <- generateRandomGLVadjMat(20, 0.2, c(0.1, 0.3, 0.3, 0.15, 0.15))
set.seed(2)
a <- generateGLVparmsFromAdj(adj, ef = 0.1)
set.seed(2)
b <- generateGLVparmsFromAdjR(adj, ef = 0.1)
all.equal(a, b)
#> [1] TRUE
```
