#' Reconstruct a brms fit object from a compressed posterior
#'
#' Replaces the posterior draws of a `brmsfit` with samples regenerated
#' from a compressed posterior, so the standard `brms` post-processing
#' methods (`predict()`, `posterior_predict()`, `pp_check()`,
#' `conditional_effects()`, ...) keep working with the lighter object.
#'
#' MCMC diagnostics (`Rhat`, `ESS`, ...) are not preserved.
#'
#' @param x The output of [compress_brmsfit()] (a list with `compressed`
#'   and `structure`), or a similarly-shaped list constructed from
#'   separately-saved objects.
#' @param n_draws Integer number of samples to regenerate. If `NULL`,
#'   uses the original number of draws.
#'
#' @return A `brmsfit` object with regenerated draws.
#' @export
reconstruct_brmsfit <- function(x, n_draws = NULL) {
  .check_reconstruct_input(x)
  comp <- x$compressed
  fit  <- x$structure

  if (!inherits(fit, "brmsfit")) {
    stop("`x$structure` must be a brmsfit object.", call. = FALSE)
  }
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop(
      "Package 'posterior' is required to reconstruct brms fits.",
      call. = FALSE
    )
  }

  original_draws <- posterior::as_draws_array(fit)
  n_chains <- posterior::nchains(original_draws)
  n_iter   <- posterior::niterations(original_draws)
  if (is.null(n_draws)) n_draws <- n_chains * n_iter

  param_names <- posterior::variables(fit)
  samples <- sample_posterior(comp, n_draws = n_draws)

  common <- intersect(param_names, colnames(samples))
  if (length(common) < length(param_names)) {
    .message(
      "Some parameters were not in the compressed posterior and will be ",
      "missing from the reconstructed fit (",
      length(param_names) - length(common), " of ", length(param_names), ")."
    )
  }
  samples <- samples[, common, drop = FALSE]

  per_chain <- floor(n_draws / n_chains)
  used <- per_chain * n_chains
  samples <- samples[seq_len(used), , drop = FALSE]

  arr <- array(
    dim = c(per_chain, n_chains, length(common)),
    dimnames = list(
      iteration = seq_len(per_chain),
      chain     = seq_len(n_chains),
      variable  = common
    )
  )
  for (chain in seq_len(n_chains)) {
    rows <- ((chain - 1L) * per_chain + 1L):(chain * per_chain)
    arr[, chain, ] <- samples[rows, , drop = FALSE]
  }
  regenerated_draws <- posterior::as_draws_array(arr)

  fit_recon <- fit
  if (!is.null(fit_recon$fit) && inherits(fit_recon$fit, "CmdStanMCMC")) {
    fit_recon$fit <- NULL
  }
  attr(fit_recon, "regenerated_draws")  <- regenerated_draws
  attr(fit_recon, "compression_method") <- comp$method
  attr(fit_recon, "reconstructed")      <- TRUE

  fit_recon
}


#' Reconstruct an sccomp object from a compressed posterior
#'
#' The compressed posterior and regenerated samples are stored in
#' `attr(x, "fit_compressed")`. The original `attr(x, "fit")` is removed.
#'
#' @param x The output of [compress_sccomp()].
#' @param n_draws Integer number of samples to regenerate. If `NULL`,
#'   uses the original number of draws.
#'
#' @return The sccomp object with `attr(x, "fit_compressed")` populated.
#' @export
reconstruct_sccomp <- function(x, n_draws = NULL) {
  .check_reconstruct_input(x)
  comp <- x$compressed
  obj  <- x$structure

  if (is.null(n_draws)) n_draws <- comp$n_draws
  samples <- sample_posterior(comp, n_draws = n_draws)

  attr(obj, "fit") <- NULL
  attr(obj, "fit_compressed") <- list(
    compressed          = comp,
    regenerated_samples = samples,
    n_draws             = nrow(samples),
    param_names         = comp$param_names
  )
  attr(obj, "compression_method") <- comp$method
  attr(obj, "reconstructed")      <- TRUE
  obj
}


#' @keywords internal
#' @noRd
.check_reconstruct_input <- function(x) {
  if (!is.list(x) || !all(c("compressed", "structure") %in% names(x))) {
    stop(
      "`x` must be a list with 'compressed' and 'structure' components, ",
      "as returned by compress_brmsfit() or compress_sccomp().",
      call. = FALSE
    )
  }
  if (!inherits(x$compressed, "posterior_compressed")) {
    stop(
      "`x$compressed` must be a posterior_compressed object.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}
