# Simulation of a neutral network assembly procesess

This model simulates an assembly process starting from basal species,
that could be fixed with the parameter BB or allowed to migrate with
rate m. The interactions are established with probability q with no
restrictions, execpt that double interactions A-\>B, B-\>A and self
interactions are not allowed. The model simulates secondary extinctions,
when a predator has all his preys extinct it will become locally
exctinct.

## Usage

``` r
neutralNetAssembly(rho, m, q, a, time, BB = 0L)
```

## Arguments

- rho:

  metacommunity number of species

- m:

  migration rate (probability) from the meta-web

- q:

  probability of interaction

- a:

  probability of extinction

- time:

  time of simulation

- BB:

  number of initial basal species present
