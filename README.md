# sundry-R
Various R functions

## sp_uncertainty()
Plots uncertainty of a spict-estimated production curve by sampling from the posterior joint distribution of the parameters. 

 
```R
library(mvtnorm)
library(spict)
source("sp_uncertainty.R")
data(pol)
rep <- fit.spict(pol$hake)
plot_production(rep)
```
