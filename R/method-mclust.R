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
#' @param model_name Character vector passed to `modelNames` in
#'   [mclust::Mclust()] (e.g. `"VVV"` for variable, `"EEE"` for equal).
#'   When `NULL` (default) `poco` auto-selects a suitable set of model
#'   names based on the shape of `draws_mat`: for high-dimensional cases
#'   (`nrow(draws_mat) <= ncol(draws_mat)`) the spherical and diagonal
#'   models `c("EII", "VII", "EEI", "EVI", "VEI", "VVI")` are used so
#'   covariances remain identifiable; otherwise mclust's full default
#'   set is used and BIC picks the best.
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
    model_name = NULL,
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
  d <- length(param_names)
  n <- nrow(draws_mat)

  user_picked <- !is.null(model_name)
  model_name <- .resolve_mclust_model_names(model_name, n, d)
  effective_models <- .effective_mclust_model_names(model_name)

  fit <- .call_mclust(
    data = draws_mat,
    G = n_components,
    modelNames = model_name,
    verbose = verbose,
    ...
  )

  if (is.null(fit) || is.null(fit$parameters$mean)) {
    stop(
      "mclust::Mclust() failed to fit a Gaussian mixture with G = ",
      paste(n_components, collapse = ","),
      " and modelNames = c(",
      paste(shQuote(effective_models), collapse = ", "), ") on a ",
      n, " x ", d, " draws matrix.\n",
      "  Try a smaller n_components, an even simpler covariance ",
      "(e.g. model_name = 'EII' or 'VII'), or method = 'mvdens_kde'.",
      call. = FALSE
    )
  }

  G <- fit$G

  if (!user_picked) {
    .message(
      "mclust: selected model '", fit$modelName, "' with G = ", G,
      " (BIC = ", format(round(fit$bic, 2), big.mark = ","), ") ",
      "out of ", length(effective_models), " candidate models: ",
      paste(effective_models, collapse = ", "), "."
    )
  }

  means_raw <- fit$parameters$mean
  if (is.matrix(means_raw)) {
    means_mat <- means_raw
  } else {
    means_mat <- matrix(as.numeric(means_raw), nrow = d, ncol = G)
  }
  # mclust sometimes returns a G x d matrix (e.g. univariate G means as a column).
  if (nrow(means_mat) == G && ncol(means_mat) == d) {
    means_mat <- t(means_mat)
  }
  if (nrow(means_mat) != d || ncol(means_mat) != G) {
    stop(
      "Unrecognised mclust mean dimensions: got ", nrow(means_mat), " x ",
      ncol(means_mat), ", expected ", d, " x ", G, ".",
      call. = FALSE
    )
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

  cov_storage <- .compact_covariances(cov_arr)

  out <- list(
    method          = "mclust",
    param_names     = param_names,
    n_params        = d,
    n_components    = G,
    weights         = fit$parameters$pro,
    means           = means_mat,
    covariances     = cov_storage$values,
    covariance_type = cov_storage$type,
    n_draws         = nrow(draws_mat),
    model_name      = fit$modelName,
    loglik          = fit$loglik,
    bic             = fit$bic
  )
  class(out) <- c("posterior_compressed_mclust", "posterior_compressed", "list")
  out
}


#' @keywords internal
#' @noRd
.sample_mclust <- function(comp, n_draws = NULL) {
  if (is.null(n_draws)) n_draws <- comp$n_draws
  n_draws <- as.integer(n_draws)

  cov_type <- .resolve_covariance_type(comp)

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

    if (cov_type == "diagonal") {
      # Each parameter independent: avoid building the d x d sigma entirely.
      sd_k <- sqrt(comp$covariances[, k])
      z <- matrix(stats::rnorm(nk * n_params), nrow = nk, ncol = n_params)
      samples[idx, ] <- z *
        matrix(sd_k,   nrow = nk, ncol = n_params, byrow = TRUE) +
        matrix(mean_k, nrow = nk, ncol = n_params, byrow = TRUE)
    } else {
      sigma_k <- .mclust_sigma(comp, k, cov_type)
      samples[idx, ] <- mvtnorm::rmvnorm(n = nk, mean = mean_k, sigma = sigma_k)
    }
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

  cov_type <- .resolve_covariance_type(comp)

  n_pts <- nrow(x)
  d <- ncol(x)
  densities <- numeric(n_pts)

  for (k in seq_len(comp$n_components)) {
    weight_k <- comp$weights[k]
    mean_k <- comp$means[, k]

    if (cov_type == "diagonal") {
      # Independent normals: log p_k(x) = sum_j log dnorm(x[, j], mu_j, sd_j)
      sd_k <- sqrt(comp$covariances[, k])
      log_pk <- rowSums(stats::dnorm(
        x,
        mean = matrix(mean_k, nrow = n_pts, ncol = d, byrow = TRUE),
        sd   = matrix(sd_k,   nrow = n_pts, ncol = d, byrow = TRUE),
        log  = TRUE
      ))
      densities <- densities + weight_k * exp(log_pk)
    } else {
      sigma_k <- .mclust_sigma(comp, k, cov_type)
      densities <- densities +
        weight_k * mvtnorm::dmvnorm(x, mean = mean_k, sigma = sigma_k)
    }
  }

  if (log) log(densities) else densities
}


#' Detect whether mclust's d x d x G sigma is effectively diagonal across
#' every component, and if so collapse to a d x G matrix of diagonals.
#'
#' This is what makes axis-aligned mclust models (`EII`/`VII`/`EEI`/
#' `VEI`/`EVI`/`VVI`) usable as a high-dimensional compression: storage
#' goes from `d^2 G` to `d G`.
#'
#' @return `list(values, type)` where `type` is one of `"full"`,
#'   `"diagonal"`, or `"shared"`.
#' @keywords internal
#' @noRd
.compact_covariances <- function(sigma) {
  if (is.null(sigma)) {
    return(list(values = sigma, type = "full"))
  }
  if (is.matrix(sigma)) {
    # Single shared d x d covariance reused across components (legacy).
    return(list(values = sigma, type = "shared"))
  }
  if (!is.array(sigma) || length(dim(sigma)) != 3L) {
    return(list(values = sigma, type = "full"))
  }
  d <- dim(sigma)[1]
  G <- dim(sigma)[3]
  if (d <= 1L) {
    # Univariate: store as diagonal (1 x G) so sampling uses the diagonal
    # path. If we kept a 1 x 1 x G array as "full", sigma[, , k] drops to a
    # scalar and mvtnorm::rmvnorm / dmvnorm call isSymmetric() on a vector.
    diag_mat <- matrix(NA_real_, nrow = 1L, ncol = G)
    for (kk in seq_len(G)) {
      diag_mat[1L, kk] <- sigma[1L, 1L, kk]
    }
    return(list(values = diag_mat, type = "diagonal"))
  }
  diag_mat <- matrix(NA_real_, nrow = d, ncol = G)
  for (k in seq_len(G)) {
    s_k <- sigma[, , k]
    diag_k <- diag(s_k)
    tol <- 1e-10 * max(1, max(abs(diag_k)))
    off_max <- max(abs(s_k - diag(diag_k, nrow = d)))
    if (off_max > tol) {
      return(list(values = sigma, type = "full"))
    }
    diag_mat[, k] <- diag_k
  }
  list(values = diag_mat, type = "diagonal")
}


#' Resolve the covariance storage type of a compressed object,
#' inferring from shape for legacy objects that lack `covariance_type`.
#' @keywords internal
#' @noRd
.resolve_covariance_type <- function(comp) {
  ct <- comp$covariance_type
  if (!is.null(ct)) return(ct)
  if (is.array(comp$covariances) && length(dim(comp$covariances)) == 3L) {
    return("full")
  }
  if (is.matrix(comp$covariances)) {
    d <- length(comp$param_names)
    G <- comp$n_components
    # Per-component diagonal storage is d x G (including univariate d = 1).
    # Legacy objects saved without `covariance_type` used to treat any matrix
    # as "shared", which breaks univariate 1 x G blocks.
    if (
      !is.null(G) &&
      nrow(comp$covariances) == d &&
      ncol(comp$covariances) == G
    ) {
      if (d == 1L || !identical(d, G) || !isSymmetric(comp$covariances)) {
        return("diagonal")
      }
    }
    if (
      nrow(comp$covariances) == ncol(comp$covariances) &&
      nrow(comp$covariances) == d
    ) {
      return("shared")
    }
    stop(
      "Legacy mclust object: cannot infer covariance layout from a ",
      nrow(comp$covariances), " x ", ncol(comp$covariances),
      " matrix (n_params = ", d, ", n_components = ", G, ").",
      call. = FALSE
    )
  }
  stop("Unknown mclust covariance structure.", call. = FALSE)
}


#' Extract the kth covariance matrix regardless of mclust storage shape.
#' @keywords internal
#' @noRd
.mclust_sigma <- function(comp, k, cov_type = .resolve_covariance_type(comp)) {
  switch(
    cov_type,
    full     = {
      sig <- comp$covariances[, , k]
      if (is.matrix(sig)) {
        sig
      } else {
        d0 <- nrow(comp$means)
        matrix(as.numeric(sig), nrow = d0, ncol = d0)
      }
    },
    shared   = comp$covariances,
    diagonal = diag(comp$covariances[, k], nrow = nrow(comp$covariances)),
    stop("Unknown mclust covariance_type: ", cov_type, call. = FALSE)
  )
}


#' Call mclust::Mclust() in a way that survives mclust's internal
#' eval(parse()) lookups when poco itself has not attached mclust.
#' @keywords internal
#' @noRd
.call_mclust <- function(...) {
  fn <- get("Mclust", envir = asNamespace("mclust"))
  do.call(fn, list(...), envir = asNamespace("mclust"))
}


#' Resolve an mclust modelNames argument given the data shape.
#'
#' Defaults follow mclust conventions but auto-restrict to spherical and
#' diagonal models when the per-component sample size is too small for
#' full covariance estimation (`n <= d`). Whenever poco picks a set on
#' the user's behalf, the choice is reported via `.message()` so users
#' know the auto-selection happened.
#' @keywords internal
#' @noRd
.resolve_mclust_model_names <- function(model_name, n, d) {
  diag_set <- c("EII", "VII", "EEI", "EVI", "VEI", "VVI")

  if (d == 1L) {
    if (is.null(model_name)) {
      .message(
        "mclust: trying univariate models c('E', 'V') ",
        "and picking the best by BIC (n = ", n, ", d = ", d, ")."
      )
      return(c("E", "V"))
    }
    uni <- intersect(model_name, c("E", "V"))
    if (length(uni)) {
      return(uni)
    }
    # Multivariate-only names (e.g. diagonal EII/VII/... from blockwise remainder)
    # do not apply when d = 1; mclust univariate mixtures use E and V.
    .message(
      "mclust: univariate data (d = 1); ",
      "using c('E', 'V') instead of multivariate model names ",
      "(n = ", n, ")."
    )
    return(c("E", "V"))
  }

  if (is.null(model_name)) {
    if (n <= d) {
      .message(
        "mclust: high-dim case (n = ", n, " <= d = ", d, "); ",
        "trying spherical/diagonal covariance models c(",
        paste(shQuote(diag_set), collapse = ", "), ") ",
        "and picking the best by BIC. ",
        "Pass `model_name = ` to override."
      )
      return(diag_set)
    }
    full_set <- .effective_mclust_model_names(NULL)
    .message(
      "mclust: trying all ", length(full_set),
      " covariance models c(", paste(full_set, collapse = ", "),
      ") and picking the best by BIC (n = ", n, ", d = ", d, ")."
    )
    return(NULL)
  }

  model_name
}


#' Resolve the actual modelNames mclust will try, given our argument.
#'
#' When we pass `NULL` to mclust it falls back to its package-level
#' default. We use `mclust::mclust.options("emModelNames")` so the
#' message can list the concrete names rather than "NULL".
#' @keywords internal
#' @noRd
.effective_mclust_model_names <- function(model_name) {
  if (!is.null(model_name)) {
    return(model_name)
  }
  defaults <- tryCatch(
    mclust::mclust.options("emModelNames"),
    error = function(e) NULL
  )
  if (is.null(defaults) || !length(defaults)) {
    defaults <- c(
      "EII", "VII", "EEI", "VEI", "EVI", "VVI",
      "EEE", "VEE", "EVE", "VVE",
      "EEV", "VEV", "EVV", "VVV"
    )
  }
  defaults
}
