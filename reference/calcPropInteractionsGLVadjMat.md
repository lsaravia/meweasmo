# Calculate the proportion of interactions types in a Lotka-Volterra adjacency matrix

Calculate the proportion of interactions types in a Lotka-Volterra
adjacency matrix

## Usage

``` r
calcPropInteractionsGLVadjMat(adjM, spc)
```

## Arguments

- adjM:

  Numeric Matrix Lotka-Volterra adjcency matrix

- sp:

  Numeric vector with the present species, if species=0 the interaction
  is not counted

## Value

Numeric vector with 5 elements Matrix 1=Mutualistic (+,+),
2=Predator-prey (-,+), 3=Competitive (-,-) 4=Comensalism (+,0),
5=Amensalism (-,0). Proportions are computed over directed matrix
entries, so symmetric interactions (mutualism and competition, two
entries each) weigh twice as much as the other types (one entry each).
