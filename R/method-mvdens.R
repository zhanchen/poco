#' mvdens GMM backend
#'
#' Internal worker that fits a Gaussian mixture model with the `mvdens`
#' package (`mvdens::fit.gmm`) and stores the mixture parameters.
#'
#' Most users should call [compress_posterior()] with
#' `method = "mvdens_gmm"`.
#'
#' @keywords internal
#' @noRd
.compress_mvdens_gmm <- function(
    draws_mat,
    n_components = 3L,
    verbose = FALSE,
    ...) {
  if (!requireNamespace("mvdens", quietly = TRUE)) {
    stop(
      "Method 'mvdens_gmm' requires the 'mvdens' package.\n",
      "  remotes::install_github('NKI-CCB/mvdens')",
      call. = FALSE
    )
  }

  param_names <- colnames(draws_mat)
  fit <- mvdens::fit.gmm(draws_mat, K = n_components, verbose = verbose, ...)

  out <- list(
    method       = "mvdens_gmm",
    param_names  = param_names,
    n_params     = length(param_names),
    n_components = n_components,
    mvdens_fit   = fit,
    n_draws      = nrow(draws_mat)
  )
  class(out) <- c("posterior_compressed_mvdens_gmm",
                  "posterior_compressed_mvdens",
                  "posterior_compressed", "list")
  out
}


#' mvdens KDE backend
#'
#' @keywords internal
#' @noRd
.compress_mvdens_kde <- function(
    draws_mat,
    verbose = FALSE,
    ...) {
  if (!requireNamespace("mvdens", quietly = TRUE)) {
    stop(
      "Method 'mvdens_kde' requires the 'mvdens' package.\n",
      "  remotes::install_github('NKI-CCB/mvdens')",
      call. = FALSE
    )
  }

  param_names <- colnames(draws_mat)
  fit <- mvdens::fit.kde(draws_mat, verbose = verbose, ...)

  out <- list(
    method       = "mvdens_kde",
    param_names  = param_names,
    n_params     = length(param_names),
    mvdens_fit   = fit,
    n_draws      = nrow(draws_mat)
  )
  class(out) <- c("posterior_compressed_mvdens_kde",
                  "posterior_compressed_mvdens",
                  "posterior_compressed", "list")
  out
}


#' @keywords internal
#' @noRd
.sample_mvdens_gmm <- function(comp, n_draws = NULL) {
  if (is.null(n_draws)) n_draws <- comp$n_draws
  n_draws <- as.integer(n_draws)

  fit <- comp$mvdens_fit
  weights <- fit$proportions
  means   <- fit$centers       # K x d
  sigmas  <- fit$covariances   # list

  K <- length(weights)
  components <- sample.int(K, size = n_draws, replace = TRUE, prob = weights)

  n_params <- length(comp$param_names)
  samples <- matrix(NA_real_, nrow = n_draws, ncol = n_params)
  for (k in seq_len(K)) {
    idx <- which(components == k)
    nk <- length(idx)
    if (nk == 0L) next
    samples[idx, ] <- mvtnorm::rmvnorm(
      n = nk, mean = means[k, ], sigma = sigmas[[k]]
    )
  }
  colnames(samples) <- comp$param_names
  samples
}


#' @keywords internal
#' @noRd
.sample_mvdens_kde <- function(comp, n_draws = NULL) {
  if (is.null(n_draws)) n_draws <- comp$n_draws
  n_draws <- as.integer(n_draws)

  fit <- comp$mvdens_fit
  data_pts <- fit$x
  H <- fit$H

  idx <- sample.int(nrow(data_pts), size = n_draws, replace = TRUE)
  noise <- mvtnorm::rmvnorm(
    n = n_draws,
    mean = rep(0, ncol(data_pts)),
    sigma = H
  )
  samples <- data_pts[idx, , drop = FALSE] + noise
  colnames(samples) <- comp$param_names
  samples
}


#' @keywords internal
#' @noRd
.density_mvdens <- function(comp, x, log = FALSE) {
  if (!is.matrix(x)) {
    if (is.vector(x)) x <- matrix(x, nrow = 1L) else x <- as.matrix(x)
  }
  d <- mvdens::mvd.pdf(comp$mvdens_fit, x)
  if (log) log(d) else d
}
