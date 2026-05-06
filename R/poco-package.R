#' poco (POsterior COmpression): Compress Bayesian Posteriors
#'
#' @description
#' `poco` (POsterior COmpression, "little" in Italian) provides a small, unified API for compressing MCMC
#' posterior draws into compact parametric or non-parametric density
#' approximations, and for regenerating samples or evaluating density on
#' demand.
#'
#' The package exposes three entry points that all share the same
#' `method` argument:
#'
#' * [compress_posterior()] for a posterior draws matrix or any object
#'   convertible via the `posterior` package.
#' * [compress_fit()] for a `cmdstanr` `CmdStanMCMC` fit.
#' * [compress_brmsfit()] for a `brms::brmsfit` fit (cmdstanr backend).
#'
#' A convenience wrapper [compress_sccomp()] is also provided for `sccomp`
#' fits, which store the underlying cmdstanr fit in `attr(x, "fit")`.
#'
#' Currently supported methods (see `?compression_methods`):
#'
#' * `"mclust"`: Gaussian mixture model via [mclust::Mclust()].
#' * `"mvdens_gmm"`: Gaussian mixture model via the `mvdens` package
#'   (suggested).
#' * `"mvdens_kde"`: kernel density estimate via the `mvdens` package
#'   (suggested).
#'
#' All compress wrappers return an S3 object of class
#' `posterior_compressed` (with a method-specific subclass) that can be
#' fed to [sample_posterior()] and [density_posterior()].
#'
#' @keywords internal
#' @importFrom utils head object.size
#' @importFrom stats cov sd
#' @importFrom mclust Mclust mclustBIC
#' @importFrom mvtnorm rmvnorm dmvnorm
"_PACKAGE"

## quiets concerns of R CMD check re: the .'s that appear in pipelines
utils::globalVariables(c("."))
