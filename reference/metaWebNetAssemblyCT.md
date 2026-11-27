# Simulation of an Assembly process from a Meta-web (continuous time)

This is a dynamical model of colonization and extinction process \[1\],
similar to the TTIB model described in \[2\]. The interactions of the
meta-web are conserved in the local web. With an additional secondary
extinction rate that controls when a predator goes extinct if it has no
prey.

## Usage

``` r
metaWebNetAssemblyCT(metaW, m, e, se, time)
```

## Arguments

- metaW:

  metacommunity adyacency matrix

- m:

  A numeric vector of species' migration rates (events per unit of time)
  from the meta-web

- e:

  A numeric vector of species' extinction rates

- se:

  A numeric vector of species' secondary extinction rates

- time:

  Simulated time, the state is recorded at integer times 1 to time-1

## Value

A list with the number of species with at least one link by time S, the
number of links by time L, the time series of species presence STime
(one row per species, one column per recorded time) and the adjacency
matrix A with effective links.

## Details

This is the continuous time version of
[`metaWebNetAssembly()`](https://lsaravia.github.io/meweasmo/reference/metaWebNetAssembly.md),
simulated with the Gillespie algorithm \[3\]. There are 3 possible
events: migration with rate *m*, local extinction with rate *e*, and
secondary extinction, when a predator has no prey, with rate *se*. The
parameters are vectors so each species could have their own parameter
values. At each iteration only one event happens: the species and the
event type are chosen with probability proportional to the rates, and
the waiting time until the event is drawn from an exponential
distribution with parameter equal to the sum of all rates. The state of
the system is recorded at integer times 1, 2, ..., *time*-1.

The main differences with the discrete time version
[`metaWebNetAssembly()`](https://lsaravia.github.io/meweasmo/reference/metaWebNetAssembly.md)
are:

- *m*, *e* and *se* are rates (events per unit of time) and can be
  larger than 1, instead of probabilities per time step in \[0,1\].

- Events are asynchronous, one at a time, instead of all species being
  updated in each time step.

- Secondary extinctions are not resolved in cascade within a step: a
  predator without prey goes extinct only when its own secondary
  extinction event happens, so cascades unfold over time.

- The time series are samples of the state at integer times, and the
  number of events between samples is variable.

## References

1.  Saravia, L. A., Marina, T. I., Kristensen, N. P., De Troch, M., &
    Momo, F. R. (2022). Ecological network assembly: How the regional
    metaweb influences local food webs. Journal of Animal Ecology,
    91(3), 630–642. https://doi.org/10.1111/1365-2656.13652

2.  Galiana, N., Lurgi, M., Claramunt-López, B., Fortin, M.-J., Leroux,
    S., Cazelles, K., et al. (2018). The spatial scaling of species
    interaction networks. Nat. Ecol. Evol., 2, 782–790

3.  Gillespie, D. T. (1976). A general method for numerically simulating
    the stochastic time evolution of coupled chemical reactions. Journal
    of Computational Physics, 22(4), 403–434. doi:
    10.1016/0021-9991(76)90041-3

## See also

[`metaWebNetAssembly()`](https://lsaravia.github.io/meweasmo/reference/metaWebNetAssembly.md)
for the discrete time version.
