# l1ideal: L1 Norm Multidimensional Ideal Point Estimation

## Authors
- [Sooahn Shin](http://sooahnshin.com/), sooahnshin@g.harvard.edu (Creator)
- Johan Lim
- Jong Hee Park

## Paper
Sooahn Shin, Johan Lim, and Jong Hee Park. "l1-based Bayesian Ideal Point Model for Multidimensional Politics" Working Paper.

## Installation

``` r
library(devtools)
install_github("sooahnshin/l1ideal", dependencies = TRUE, ref = "master")
```

## Examples
``` r
library(l1ideal)
dat <- generate.data(seed = 1)
rc <- pscl::rollcall(dat$votes,
                     yea = 1, nay = 0, missing = 2, notInLegis = 3,
                     legis.names = dat$legis_data$name, legis.data = dat$legis_data,
                     vote.names = dat$votes_data$name, vote.data = dat$votes_data)
res <- l1ideal(rc, dimensions = 2, mcmc  = 100, thin = 10, burnin = 100, minvotes = 20,
               lop = 0, polarity = c(86,97), verbose = 100, seed = 123)
p <- plot.simulation(dat, res)
multiplot(plotlist = p, cols = 4)
```
