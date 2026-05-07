# poco

<!-- badges: start -->
[![Lifecycle:experimental](https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

<p align="left">
  <img src="inst/poco_logo.png" alt="poco logo" style="width: 200px; border-radius: 12px; border: 1px solid rgba(159, 122, 234, 0.35);" />
</p>

---

`poco` (POsterior COmpression, "little" in Italian) provides a small, unified API for compressing MCMC
posterior draws into compact density approximations (Gaussian mixture
models, kernel density estimates), regenerating samples on demand, and
evaluating the density at arbitrary points.

It is designed for situations where you need to store, ship, or
post-process *many* Bayesian models with limited disk space - for example
when caching hundreds of `brms` or `sccomp` fits, or sharing pre-trained
posteriors with collaborators. With Gaussian-mixture compression, the
on-disk footprint can drop by 100-1000x while preserving correlation
structure, multimodality and most downstream summaries.

## Why?

A Stan model with 1000 parameters x 8000 draws x 8 bytes = ~64 MB per
chain. With four chains that is 256 MB, and 100 such models is 25 GB.

A 10-component Gaussian mixture in 1000 dimensions stores only
`weights + means + covariances`, typically a few hundred KB.

## Installation

```r
# install.packages("remotes")
remotes::install_github("MangiolaLaboratory/poco")
```

`brms` / `cmdstanr` / `sccomp` are *Suggested* dependencies; you only
need them if you plan to compress fits from those packages.

## Quick tour

`poco` exposes three entry points that all share the same
`method` argument:

| Function              | Input                                  |
| --------------------- | -------------------------------------- |
| `compress_posterior()` | a draws matrix / data.frame / `posterior::draws_*` |
| `compress_fit()`       | a `cmdstanr::CmdStanMCMC` or vector of CSV paths |
| `compress_brmsfit()`   | a `brms::brmsfit` (cmdstanr backend) |
| `compress_sccomp()`    | an `sccomp` fit (cmdstanr backend) |

All return an S3 object with class `posterior_compressed`, which is
consumed by:

- `sample_posterior()` - regenerate draws (any number),
- `density_posterior()` - evaluate p.d.f. (or log p.d.f.) at given points,
- `evaluate_compression()` - quantify how much information was lost
  (a "98.3% reproduction"-style score, see below),
- `reconstruct_brmsfit()` / `reconstruct_sccomp()` - rebuild a usable
  fit object so the standard post-processing keeps working.

### From a draws matrix

```r
library(poco)

set.seed(1)
draws <- matrix(rnorm(2000 * 3), ncol = 3,
                dimnames = list(NULL, c("alpha", "beta", "sigma")))

comp <- compress_posterior(
  draws,
  method = "mclust",        # or "mvdens_gmm", "mvdens_kde"
  n_components = 5
)

new_draws <- sample_posterior(comp, n_draws = 4000)
density_posterior(comp, x = matrix(c(0, 0, 1), nrow = 1))

saveRDS(comp, "model.rds", compress = "xz")
```

### From a brms fit

```r
library(brms)
fit <- brm(y ~ x + (1 | g), data = dat, backend = "cmdstanr")

result <- compress_brmsfit(fit, method = "mclust", n_components = 5)

saveRDS(result$compressed, "model_compressed.rds", compress = "xz")
saveRDS(result$structure,  "model_structure.rds")

# Later, reload + reconstruct + use brms machinery as normal:
result_loaded <- list(
  compressed = readRDS("model_compressed.rds"),
  structure  = readRDS("model_structure.rds")
)
fit_recon <- reconstruct_brmsfit(result_loaded)

predict(fit_recon, newdata = new_data)
posterior_predict(fit_recon)
```

### From a cmdstanr fit (or CSV files)

```r
library(cmdstanr)
fit <- mod$sample(data = stan_data, chains = 4)

comp <- compress_fit(fit, method = "mclust", n_components = 5)
new_draws <- sample_posterior(comp, n_draws = 4000)
```

## Methods

| `method`        | Backend                       | Notes |
| --------------- | ----------------------------- | ----- |
| `"mclust"`      | `mclust::Mclust()`            | parametric GMM, BIC-based selection, fast & compact |
| `"mvdens_gmm"`  | `mvdens::fit.gmm()`           | EM-based GMM, requires the suggested `mvdens` package |
| `"mvdens_kde"`  | `mvdens::fit.kde()`           | kernel density estimate, larger but non-parametric |

`compression_methods()` returns the current set.

## How much information was lost?

Like JPEG quality scores, `evaluate_compression()` reports a
`reproduction %` for the compressed posterior. It uses two
*backend-independent* two-sample diagnostics on samples regenerated
from the compressed object versus the reference draws — so the score
is not circular (it doesn't reuse the same density model that did the
compression):

- **Energy distance** (Szekely & Rizzo) anchored against a bootstrap
  Monte Carlo noise envelope from the reference itself.
- **Classifier two-sample test (C2ST)** — a random forest via
  `ranger::ranger()` under cross-validation tries to tell reference
  draws from reconstructed draws; AUC of 0.5 means indistinguishable.
  Pass `classifier = "knn"` for a k-NN alternative instead of the forest.

```r
fidelity <- evaluate_compression(comp, reference_draws = draws)
fidelity
#> <compression_fidelity>
#>   method        : mclust
#>   parameters    : 3
#>   reference n   : 1500
#>   eval n        : 1500
#>   ----------------------------------------
#>   energy        : 100.0% reproduction
#>     distance    : 0.054   noise envelope (q90): 0.062   ratio: 0.87x
#>   C2ST          :  97.7% reproduction
#>     AUC         : 0.511   classifier: ranger   cv_folds: 5
#>   ----------------------------------------
#>   reproduction  : 98.9%
```

For a strict out-of-sample evaluation, hold the reference draws out
*before* fitting:

```r
idx  <- sample.int(nrow(draws), 0.8 * nrow(draws))
comp <- compress_posterior(draws[idx, ], method = "mclust")
evaluate_compression(comp, reference_draws = draws[-idx, ])
```

## Vignettes

- `vignette("introduction", package = "poco")` - hello world,
  matrix in / samples out.
- `vignette("compress-brms", package = "poco")` - end-to-end
  `brms` workflow, including reconstruction.
- `vignette("compress-funnel", package = "poco")` - `cmdstanr`
  workflow on Neal's funnel.
- `vignette("compress-sccomp", package = "poco")` - end-to-end
  `sccomp` workflow.
- `vignette("methods-benchmark", package = "poco")` - side-by-side
  comparison of the supported compression methods.
- `vignette("evaluate-compression", package = "poco")` - measure
  reproduction quality with `evaluate_compression()` (sanity checks
  on lost correlation, mode collapse, and a sliding-difficulty sweep).

## License

GPL-3
