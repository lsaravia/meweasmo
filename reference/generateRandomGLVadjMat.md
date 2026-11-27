# Generate random Lotka-Volterra adjacency matrix with a given proportion of interactions (approximately)

Generate random Lotka-Volterra adjacency matrix with a given proportion
of interactions (approximately)

## Usage

``` r
generateRandomGLVadjMat(numSp, C, propInt)
```

## Arguments

- numSp:

  Integer, number of species

- C:

  Double, Probability of interaction = Connectivity = number of
  interacions/(numSp^2)

- propInt:

  Numeric vector, Proportion of interactions, there are 5 types of
  interactions the elements of the vector represent 1= mutualism (+,+),
  2=Predator-prey (+,-), 3= competition (-,-), 4=Comensalims (+,0),
  5=Amensalism (-,0).

## Value

A NumericMatrix with the GLV type of adjacency matrix
