# NBA ordinal graphical mixture

This directory contains the NBA application of the mixture of finite mixtures
probit graphical model (MFM--PGM). The archived posterior sample contains 20,000 MCMC iterations, and the reported summaries use iterations 4,001--20,000.

## Start here

1. Download this directory and open `NBA_chain1.Rproj` in RStudio.
2. In the R console, run `source("run.R")`.
3. The default `action <- "summary"` recomputes the Dahl partition and graph
   summaries from the supplied posterior draws. **It does not start MCMC** and
   does not require Rtools or a C++ compilation. This calculation can take several
   minutes. New output is written under `runs/summary_<timestamp>_<pid>/`.

The already computed results are available immediately in
`results/reported_chain1/`. Nothing in that directory is overwritten by `run.R`.

## What the files do

| File | Purpose |
| --- | --- |
| `run.R` | Main entry point for summary, sampling, continuation and preprocessing. |
| `model.R`, `model.cpp` | Model and MCMC sampler implementation. |
| `chain.R` | MCMC driver and checkpoint support. |
| `preprocess.R` | Season-record selection and ordinal coding. |
| `analysis.R` | Player-profile, position, graph and iteration-window summaries. |
| `raw/` | Spreadsheet used in the NBA analysis. |
| `processed/` | Ordinal input, numerical measurements, player metadata and processing audit. |
| `results/reported_chain1/` | Posterior draws, reported summaries and diagnostic tables. |
| `provenance/` | Model settings, RNG information, initial state, source hashes and R session information. |


## Run the chain from the beginning

In the R console, install missing packages once:

```r
install.packages(c("Rcpp", "RcppArmadillo", "MASS", "mclust", "readxl"))
```

On Windows, install Rtools compatible with your R version to compile `model.cpp`.
In `run.R`, set:

```r
action <- "sample"
target_iterations <- 20000L
checkpoint_every <- 4000L
```

Then run `source("run.R")`. 
The automatic initializer used normal scores and `mclust` EEE models
with 1--6 components. 

The original sampling environment was R 4.3.1 on Windows, with MASS 7.3-60,
mclust 6.1.3, Rcpp 1.1.2 and RcppArmadillo 15.6.0-1, as recorded in the supplied
session logs.  The supplied draws reproduce
the **reported** summaries without requiring a new run or identical environments.
