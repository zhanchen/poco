#' Sample from a compressed posterior
#'
#' Draws fresh posterior samples from an object returned by
#' [compress_posterior()] or its wrappers.
#'
#' @param object A `posterior_compressed` object (single-block or
#'   `posterior_compressed_blockwise`), or a path to an `.rds` file
#'   containing one.
#' @param n_draws Integer number of samples to draw. If `NULL` (default),
#'   uses the original number of MCMC draws.
#' @param ... Backend-specific options. Currently the only honoured option
#'   is `verbose = TRUE` for blockwise objects, which prints a per-block
#'   sampling progress message.
#'
#' @return A numeric matrix of draws (rows = draws, cols = parameters).
#' @export
#'
#' @examples
#' set.seed(1)
#' draws <- matrix(rnorm(2000 * 2), ncol = 2,
#'                 dimnames = list(NULL, c("a", "b")))
#' comp <- compress_posterior(draws, method = "mclust", n_components = 2)
#' s <- sample_posterior(comp, n_draws = 500)
#' dim(s)
sample_posterior <- function(object, n_draws = NULL, ...) {
  object <- .resolve_compressed(object)
  if (inherits(object, "posterior_compressed_blockwise")) {
    dots <- list(...)
    return(.sample_blockwise(
      object,
      n_draws = n_draws,
      verbose = isTRUE(dots$verbose)
    ))
  }
  if (inherits(object, "posterior_compressed_mclust")) {
    return(.sample_mclust(object, n_draws = n_draws))
  }
  if (inherits(object, "posterior_compressed_mvdens_gmm")) {
    return(.sample_mvdens_gmm(object, n_draws = n_draws))
  }
  if (inherits(object, "posterior_compressed_mvdens_kde")) {
    return(.sample_mvdens_kde(object, n_draws = n_draws))
  }
  stop(
    "No sample_posterior() method for object of class '",
    paste(class(object), collapse = "/"), "'.",
    call. = FALSE
  )
}


#' Evaluate the density of a compressed posterior
#'
#' Computes the density (or log-density) at one or more points using a
#' compressed posterior.
#'
#' @param object A `posterior_compressed` object, or a path to an `.rds`
#'   file containing one.
#' @param x A matrix or data.frame of points to evaluate (rows = points,
#'   cols = parameters in the same order as `object$param_names`). A
#'   plain numeric vector is treated as a single point.
#' @param log Logical; if `TRUE`, return the log-density. Default `FALSE`.
#' @param ... Currently unused.
#'
#' @return A numeric vector of densities (or log-densities), one entry
#'   per row of `x`.
#'
#' @export
density_posterior <- function(object, x, log = FALSE, ...) {
  object <- .resolve_compressed(object)
  if (inherits(object, "posterior_compressed_blockwise")) {
    return(.density_blockwise(object, x, log = log))
  }
  if (inherits(object, "posterior_compressed_mclust")) {
    return(.density_mclust(object, x, log = log))
  }
  if (inherits(object, "posterior_compressed_mvdens_gmm")) {
    return(.density_mvdens(object, x, log = log))
  }
  if (inherits(object, "posterior_compressed_mvdens_kde")) {
    return(.density_mvdens(object, x, log = log))
  }
  stop(
    "No density_posterior() method for object of class '",
    paste(class(object), collapse = "/"), "'.",
    call. = FALSE
  )
}


#' @keywords internal
#' @noRd
.resolve_compressed <- function(object) {
  if (is.character(object) && length(object) == 1L && file.exists(object)) {
    object <- readRDS(object)
  }
  if (inherits(object, "compressed_fit") && is.list(object) &&
      !is.null(object$compressed)) {
    object <- object$compressed
  }
  object
}


#' Print method for compressed posteriors
#'
#' @param x A `posterior_compressed` object.
#' @param ... Currently unused.
#' @return Invisibly returns `x`.
#' @export
print.posterior_compressed <- function(x, ...) {
  cat("<posterior_compressed:", x$method, ">\n")
  cat(" parameters: ", x$n_params, "\n", sep = "")
  if (!is.null(x$n_components)) {
    cat(" components: ", x$n_components, "\n", sep = "")
  }
  cat(" original draws: ", x$n_draws, "\n", sep = "")
  if (!is.null(x$model_name)) {
    cat(" mclust model: ", x$model_name, "\n", sep = "")
  }
  if (!is.null(x$covariance_type)) {
    cat(" covariance type: ", x$covariance_type, "\n", sep = "")
  }
  if (!is.null(x$bic)) {
    cat(" BIC: ", round(x$bic, 2), "\n", sep = "")
  }
  invisible(x)
}


#' Summary method for compressed posteriors
#'
#' @param object A `posterior_compressed` object.
#' @param ... Currently unused.
#' @return A list with class `summary.posterior_compressed`.
#' @export
summary.posterior_compressed <- function(object, ...) {
  out <- list(
    method = object$method,
    n_params = object$n_params,
    n_components = object$n_components,
    n_draws = object$n_draws,
    object_size = utils::object.size(object),
    parameters = utils::head(object$param_names, 20)
  )
  class(out) <- "summary.posterior_compressed"
  out
}

#' @export
print.summary.posterior_compressed <- function(x, ...) {
  cat("Compressed posterior\n")
  cat("  method        : ", x$method, "\n", sep = "")
  cat("  parameters    : ", x$n_params, "\n", sep = "")
  cat("  components    : ", if (is.null(x$n_components)) "(NA)"
      else x$n_components, "\n", sep = "")
  cat("  original draws: ", x$n_draws, "\n", sep = "")
  cat("  in-memory size: ", format(x$object_size, units = "auto"), "\n",
      sep = "")
  cat("  first params  : ", paste(x$parameters, collapse = ", "), "\n",
      sep = "")
  invisible(x)
}
