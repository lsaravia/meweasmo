# Simulation of a Cascade fod web assembly procesess

This model simulates an assembly process starting from basal species,
that could be fixed with the parameter BB or allowed to migrate with
rate m. The interactions are established with probability q if species
i\>j then i becomes predator and j prey (i\<–j). The model simulates
secondary extinctions, when a predator has all his preys extinct it will
become locally exctinct.

## Usage

``` r
cascadeNetAssembly(rho, m, q, a, time, BB = 0L)
```

## Arguments

- rho:

  metacommunity number of species

- m:

  migration rate (probability) from the meta-web

- q:

  probability of interaction

- a:

  extinctions of species (not proportional to C)

- time:

  time of simulation

- BB:

  number of basal species if 0 will be ~poison(rho\*m)

## Value

A list with the final the number of species by time S, the number of
links by time L, the number of basal species and the adjacency matrix A.
