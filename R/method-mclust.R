#' Compress a draws matrix with an mclust Gaussian mixture model
#'
#' Internal worker that fits a finite Gaussian mixture model to a posterior
#' draws matrix using [mclust::Mclust()] and returns a compact list
#' containing only the mixture parameters.
#'
#' Most users should call [compress_posterior()] with `method = "mclust"`
#' instead of this function directly.
#'
#' @param draws_mat Numeric matrix of draws (rows = draws,
#'   cols = parameters).
#' @param n_components Integer number of mixture components (`G`).
#'   Can be a single integer or a vector to let mclust choose by BIC.
#' @param model_name Character string passed to `modelNames` in
#'   [mclust::Mclust()] (e.g. `"VVV"` for variable, `"EEE"` for equal).
#' @param verbose Logical; print mclust progress.
#' @param ... Additional arguments passed to [mclust::Mclust()].
#'
#' @return A list with class `c("posterior_compressed_mclust",
#'   "posterior_compressed")`.
#' @keywords internal
#' @noRd
.compress_mclust <- function(
    draws_mat,
    n_components = 3L,
    model_name = "VVV",
    verbose = FALSE,
    ...) {
  if (!requireNamespace("mclust", quietly = TRUE)) {
    stop(
      "Method 'mclust' requires the 'mclust' package.\n",
      "  install.packages('mclust')",
      call. = FALSE
    )
  }

  param_names <- colnames(draws_mat)

  if (length(param_names) == 1L &&
      !model_name %in% c("E", "V")) {
    model_name <- "V"
  }

  fit <- mclust::Mclust(
    data = draws_mat,
    G = n_components,
    modelNames = model_name,
    verbose = verbose,
    ...
  )

  d <- length(param_names)
  G <- fit$G

  means_raw <- fit$parameters$mean
  if (is.matrix(means_raw)) {
    means_mat <- means_raw
  } else {
    means_mat <- matrix(means_raw, nrow = d, ncol = G)
  }

  cov_raw <- fit$parameters$variance$sigma
  if (d == 1L) {
    sigsq <- fit$parameters$variance$sigmasq
    if (is.null(sigsq)) sigsq <- as.numeric(cov_raw)
    if (length(sigsq) == 1L) sigsq <- rep(sigsq, G)
    cov_arr <- array(0, dim = c(1L, 1L, G))
    for (k in seq_len(G)) cov_arr[, , k] <- sigsq[k]
  } else if (is.array(cov_raw) && length(dim(cov_raw)) == 3L) {
    cov_arr <- cov_raw
  } else if (is.matrix(cov_raw)) {
    cov_arr <- array(cov_raw, dim = c(d, d, G))
  } else {
    stop("Unrecognised mclust covariance shape.", call. = FALSE)
  }

  out <- list(
    method       = "mclust",
    param_names  = param_names,
    n_params     = d,
    n_components = G,
    weights      = fit$parameters$pro,
    means        = means_mat,
    covariances  = cov_arr,
    n_draws      = nrow(draws_mat),
    model_name   = fit$modelName,
    loglik       = fit$loglik,
    bic          = fit$bic
  )
  class(out) <- c("posterior_compressed_mclust", "posterior_compressed", "list")
  out
}


#' @keywords internal
#' @noRd
.sample_mclust <- function(comp, n_draws = NULL) {
  if (is.null(n_draws)) n_draws <- comp$n_draws
  n_draws <- as.integer(n_draws)

  components <- sample.int(
    n = comp$n_components,
    size = n_draws,
    replace = TRUE,
    prob = comp$weights
  )

  n_params <- length(comp$param_names)
  samples <- matrix(NA_real_, nrow = n_draws, ncol = n_params)

  for (k in seq_len(comp$n_components)) {
    idx <- which(components == k)
    nk <- length(idx)
    if (nk == 0L) next

    mean_k <- comp$means[, k]
    sigma_k <- .mclust_sigma(comp$covariances, k)
    samples[idx, ] <- mvtnorm::rmvnorm(n = nk, mean = mean_k, sigma = sigma_k)
  }

  colnames(samples) <- comp$param_names
  samples
}


#' @keywords internal
#' @noRd
.density_mclust <- function(comp, x, log = FALSE) {
  if (!is.matrix(x)) {
    if (is.vector(x)) {
      x <- matrix(x, nrow = 1L)
    } else {
      x <- as.matrix(x)
    }
  }
  if (ncol(x) != length(comp$param_names)) {
    stop(
      "x has ", ncol(x), " columns but the compressed posterior has ",
      length(comp$param_names), " parameters.",
      call. = FALSE
    )
  }

  n_pts <- nrow(x)
  densities <- numeric(n_pts)

  for (k in seq_len(comp$n_components)) {
    weight_k <- comp$weights[k]
    mean_k <- comp$means[, k]
    sigma_k <- .mclust_sigma(comp$covariances, k)
    densities <- densities +
      weight_k * mvtnorm::dmvnorm(x, mean = mean_k, sigma = sigma_k)
  }

  if (log) log(densities) else densities
}


#' Extract the kth covariance matrix regardless of mclust shape
#' @keywords internal
#' @noRd
.mclust_sigma <- function(cov_obj, k) {
  if (is.array(cov_obj) && length(dim(cov_obj)) == 3L) {
    cov_obj[, , k]
  } else if (is.matrix(cov_obj)) {
    cov_obj
  } else {
    stop("Unknown mclust covariance structure.", call. = FALSE)
  }
}
