# NBA ordinal graphical mixture: reported chain 1

This directory contains the NBA application of the mixture of finite mixtures
probit graphical model (MFM--PGM). It preserves the original **chain 01**, including
its continuation to **20,000 iterations**. The reported summaries use iterations
**4,001--20,000**. It contains no other chain and launches no parallel workers.

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
| `run.R` | One entry point; choose summary, sampling, continuation or preprocessing. |
| `model.R`, `model.cpp` | Original archived model and sampler, unchanged byte for byte. |
| `chain.R` | Single-chain driver and checkpoint support around the original R loop. |
| `preprocess.R` | Original season-record selection and ordinal coding. |
| `analysis.R` | Player, profile, position, graph and iteration-window tables. |
| `raw/` | Spreadsheet supplied for the original analysis. |
| `processed/` | Ordinal input, numeric measurements, player metadata and processing audit. |
| `results/reported_chain1/` | Actual archived-chain draws, summaries and diagnostic tables. |
| `provenance/` | Original settings, RNG stream, initial state, source hashes and R sessions. |
| `checkpoints/` | Destination for the optional full checkpoint downloaded from a release. |

The historical generic defaults inside `model_settings()` are unchanged.
The NBA entry point uses the archived NBA settings in
`provenance/chain1_initial_manifest.rds`, including **gamma = 1**, not the generic
default of 0.1. Use `run.R` rather than calling the generic defaults directly.

## Reproduce preprocessing

Set `action <- "preprocess"` in `run.R`, then source the file. Install `readxl` if
needed. Processing retains a player's full-season TOT row where present, removes
team-stint duplicates, excludes missing/nonfinite retained measurements, and
removes WS and BPM from the model inputs. There is no minimum-playing-time filter.
The 664 original rows represent 540 players; 124 duplicate stint rows and three
players with missing retained values are excluded, leaving 537 players and 18
indicators. Each indicator is categorised using type-7 sample quintiles, with
cutpoint equality assigned to the lower category and ties kept together.
These observed-data cutpoints differ from the latent model thresholds.

The original external data provider has not been independently established in
this package; the spreadsheet is the one supplied with the analysis. The numeric
input and row ordering are checked against the archived chain.

## Run only chain 1 from the beginning

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

Then run `source("run.R")`. This is a long calculation. The original archived run
recorded approximately 54.6 hours of cumulative sampling time; that is historical
timing, not an estimate for every computer. Checkpoints are saved every 4,000
iterations and the driver **continues automatically** to the target; checkpoint
saving is not a pause. Progress is in `runs/.../chain_01/progress.txt`.

The driver restores the **actual saved pre-fit RNG stream** from the original
doRNG task. It does not replace that stream with `set.seed()` and does not run
chains 2--4. The automatic initializer used normal scores and `mclust` EEE models
with 1--6 components. The original initial partition had six groups of sizes
68, 86, 17, 79, 71 and 216. The archived initial labels/state are included for
verification; they are not posterior conclusions. The driver checks initialization
against this state and records any discrepancy.

The original sampling environment was R 4.3.1 on Windows, with MASS 7.3-60,
mclust 6.1.3, Rcpp 1.1.2 and RcppArmadillo 15.6.0-1, as recorded in the supplied
session logs. Those are archived version records, not a claim that currently
installed packages are identical. Different versions, compilers and numerical
libraries can produce different rerun trajectories. The supplied draws reproduce
the **reported** summaries without requiring a new run or identical environments.

## Resume after interruption or extend the archived chain

For a new run interrupted on your computer, point `source_checkpoint` in `run.R`
to its latest fully written `checkpoint_XXXXX.rds` file. Do not use a `.partial`
file. For the original reported run, download the separate full-checkpoint release
and extract `checkpoint_20000.rds` into `checkpoints/`. Set:

```r
action <- "resume"
source_checkpoint <- "checkpoints/checkpoint_20000.rds"
target_iterations <- 24000L
checkpoint_every <- 4000L
```

Then source `run.R`. The driver preserves the original burn-in of 4,000,
restores the saved state and RNG, and continues at the next **absolute** iteration.
It does not repeat the initial warm-up. A new run directory preserves previous
files. No continuation is started by the default summary action.

## Reported results and interpretation

- N = 537; p = 18; five levels for every model variable.
- 20,000 iterations; 4,000 discarded; 16,000 retained partitions and graphs.
- Dahl reference partition: four groups, sizes **117, 25, 339, 56**; selected
  iteration **12,476**.
- Occupied-count frequencies in retained draws: K+ = 4: **0.666375**;
  K+ = 5: **0.333625**. These refer to K+, not the total latent mixture count K.
- Graph summaries use the original greedy one-to-one Jaccard matching with a
  minimum similarity of 0.3. Frequencies divide by matched draws; edges with
  frequency **at least 0.5** are selected. The four graphs contain **28, 11, 58,
  25** edges. These are latent conditional-dependence graphs, not causal graphs.

`posterior_draws.rds` keeps every partition and every retained graph needed by the
unchanged summarizer. It omits the large latent/parameter snapshots and final RNG
state, so it is **not a continuation checkpoint**. The separate full checkpoint
preserves the complete original fit and supports continuation.

Window-comparison tables and selected parameter traces from the previous analysis
are included alongside the summaries. They document a finite single-chain run;
the package does not claim that all parameters or graph summaries have converged.
`single_chain_split_diagnostics.csv` uses splits of one chain, not independent
chains. The archived graph/threshold window diagnostics are preserved.

## Scope and provenance

This is the **NBA MFM analysis only**. It does not contain the simulation study,
other datasets, or the separately fitted continuous-data comparator. A data/code
availability statement should describe that scope accurately.

The posterior kernel, preprocessing and matching rules are unchanged. New code
only replaces the four-chain launcher with a single-chain launcher, adds
checkpoint orchestration and exports tables. Source hashes and verification
results are supplied in `provenance/`. No new formal NBA experiment was run to
create the release. The full checkpoint is distributed separately to keep normal
repository files small.
