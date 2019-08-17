# l1ideal: L1 Norm Multidimensional Ideal Point Estimation
[![Build Status](https://travis-ci.com/sooahn/l1ideal.svg?token=yMJkyPKby9ixFy1a4BPp&branch=master)](https://travis-ci.com/sooahn/l1ideal)

## Authors
[Sooahn Shin](http://sooahnshin.com/), Yohan Lim, [Jong Hee Park](http://jhp.snu.ac.kr/)

## Paper
Sooahn Shin, Yohan Lim, and Jong Hee Park 2019. "L1 norm Based Multidimensional Ideal Point Estimation: With Application to Roll Call Voting Data" Working Paper.

## Installation

``` r
library(devtools)
install_github("sooahn/l1ideal", dependencies = TRUE, ref = "master")
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
