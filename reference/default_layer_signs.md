# Default signs of the interaction layers

A table with the layers expected by
[`build_multilayer_assembly()`](https://lsaravia.github.io/meweasmo/reference/build_multilayer_assembly.md)
when the links carry no sign attribute. `sign_from` and `sign_to` are
the signs of the effect received by the species at the origin and at the
end of each link. The defaults assume that links go from the species
with the lower sign to the species with the higher sign: trophic links
from the prey (-) to the predator (+), commensalism from the unaffected
species (0) to the beneficiary (+), amensalism from the harmed species
(-) to the unaffected one (0). Mutualism (+,+) and competition (-,-) are
symmetric, and non-trophic antagonism is coded (-,+) like a trophic
link.

## Usage

``` r
default_layer_signs()
```

## Value

a data.frame with columns `layer`, `sign_from`, `sign_to` and `rule`.

## Details

Datasets use different conventions, so check the direction of the links
and edit the table (or use a sign attribute per link, see `effect_attr`
in
[`build_multilayer_assembly()`](https://lsaravia.github.io/meweasmo/reference/build_multilayer_assembly.md))
before simulating.

`rule` is `"trophic"` for the layer that defines the meta-web and the
prey rule, and `"modifier"` for the layers whose links change the
extinction rates.

## See also

[`build_multilayer_assembly()`](https://lsaravia.github.io/meweasmo/reference/build_multilayer_assembly.md),
[`metaWebNetAssemblyML()`](https://lsaravia.github.io/meweasmo/reference/metaWebNetAssemblyML.md)

## Examples

``` r
default_layer_signs()
#>          layer sign_from sign_to     rule
#> 1      Trophic        -1       1  trophic
#> 2   Antagonism        -1       1 modifier
#> 3    Mutualism         1       1 modifier
#> 4  Competition        -1      -1 modifier
#> 5 Commensalism         0       1 modifier
#> 6   Amensalism        -1       0 modifier
```
