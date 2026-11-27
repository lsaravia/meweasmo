# Generate random parameters of a generalized Lotka-Volterra model from an adjacency matrix

Builds the parameters needed by
[`metaWebNetAssemblyGLV()`](https://lsaravia.github.io/meweasmo/reference/metaWebNetAssemblyGLV.md)
(interaction matrix, intrinsic growth rates and migration rates) from a
matrix that only defines the structure of the interactions, for example
one generated with
[`generateRandomGLVadjMat()`](https://lsaravia.github.io/meweasmo/reference/generateRandomGLVadjMat.md).
The input matrix is not modified.

## Usage

``` r
generateGLVparmsFromAdj(
  adjM,
  ef,
  predIntAvg = 0.01,
  selfLimAvg = as.numeric(c(0.01, 0.01, 0.01)),
  migrAvg = 0,
  preserveInt = 0L,
  predIntSd = 0.01,
  rndType = 0L
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

## Details

The type of each interaction is read from the signs of `adjM`
(`adjM[i, j]` is the effect of species `j` on species `i`):

- `adjM[i, j] < 0` and `adjM[j, i] > 0`: predator-prey, `j` is the
  predator and `i` the prey. The strength of the effect on the prey is
  random and the effect on the predator is `ef` times the effect on the
  prey, with opposite sign.

- both negative: competition; both positive: mutualism. Each direction
  gets its own random strength.

- `adjM[i, j] > 0` and `adjM[j, i] == 0`: commensalism (`i` benefits);
  `adjM[i, j] < 0` and `adjM[j, i] == 0`: amensalism (`i` is harmed).

The strengths of the interactions (off-diagonal values) are drawn
according to `preserveInt` and `rndType`.

The remaining parameters depend on the type of each species, which is
also deduced from `adjM`: *mutualistic* if it takes part in any
mutualistic interaction, otherwise *basal* if it has no resources (it is
never the predator in a predator-prey interaction, nor the species `j`
in a commensalism), otherwise *consumer*. For each type:

- Self-limitation (diagonal of the interaction matrix):
  `-runif(0, 2 * selfLimAvg[k])`, where `k` is 1 for mutualistic, 2 for
  basal and 3 for consumer species. The diagonal of `adjM` is ignored.

- Intrinsic growth rate `r`: positive, `runif(0, 1)`, for basal species;
  negative, `-runif(0, 1)`, for consumers (they need resources to
  persist); and of random sign, with the same magnitude, for mutualistic
  species (random obligate mutualism).

- Migration `m`: `runif(0, 2 * migrAvg)` for all species.

There is an equivalent implementation in R,
[`generateGLVparmsFromAdjR()`](https://lsaravia.github.io/meweasmo/reference/generateGLVparmsFromAdjR.md),
that gives the same results for the same random seed.

## See also

[`generateRandomGLVadjMat()`](https://lsaravia.github.io/meweasmo/reference/generateRandomGLVadjMat.md)
to generate `adjM`,
[`generateGLVparmsFromAdjR()`](https://lsaravia.github.io/meweasmo/reference/generateGLVparmsFromAdjR.md)
for the R implementation, and
[`metaWebNetAssemblyGLV()`](https://lsaravia.github.io/meweasmo/reference/metaWebNetAssemblyGLV.md).

## Examples

``` r
set.seed(1)
adj <- generateRandomGLVadjMat(20, 0.2, c(0.1, 0.3, 0.3, 0.15, 0.15))
parms <- generateGLVparmsFromAdj(adj, ef = 0.1, predIntAvg = 0.01, migrAvg = 0.01)
str(parms)
#> List of 3
#>  $ interM: num [1:20, 1:20] -0.01213 0 0 0 -0.00846 ...
#>  $ r     : num [1:20] 0.0294 -0.0432 -0.1737 -0.9557 -0.1971 ...
#>  $ m     : num [1:20] 0.00673 0.00741 0.01244 0.01307 0.00231 ...

# Gamma distributed interaction strengths
parms <- generateGLVparmsFromAdj(adj, ef = 0.1, predIntAvg = 0.01, predIntSd = 0.001, rndType = 1)
```
