# poco

<!-- badges: start -->
[![Lifecycle:experimental](https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

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

## Vignettes

- `vignette("introduction", package = "poco")` - Hello world,
  matrix in / samples out.
- `vignette("compress-brms", package = "poco")` - end-to-end
  brms workflow, including reconstruction.
- `vignette("compress-funnel", package = "poco")` - cmdstanr
  workflow on Neal's funnel.

## License

GPL-3
