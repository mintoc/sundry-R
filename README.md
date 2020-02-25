# sundry-R
Various R functions

## sp_uncertainty()
Plots uncertainty of a [spict](https://github.com/DTUAqua/spict)-estimated production curve by sampling from the joint posterior distribution of the parameters. Written by Cóilín Minto and Paul Bouch.

Plotting of empirical surplus-production inspired by [Hilborn (2001)](https://www.nrcresearchpress.com/doi/10.1139/f01-018#.XlVLP3X7RhE).
 
```R
library(mvtnorm)
library(spict)
source("sp_uncertainty.R")
data(pol)
rep <- fit.spict(pol$hake)
sp_uncertainty(rep)
```
<img src = "./figures/sp_uncertainty.png" width = 600, class="center">
