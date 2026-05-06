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
#' @param n_draws Integer number of samples to regenerate. Defaults to
#'   the original number of draws stored in `x$structure`.
#'
#' @return A `brmsfit` object with regenerated draws.
#' @export
reconstruct_brmsfit <- function(x, n_draws = get_n_draws_from_fit(x$structure)) {
  .check_reconstruct_input(x)
  comp <- x$compressed
  fit_structure <- x$structure

  if (!inherits(fit_structure, "brmsfit")) {
    stop("`x$structure` must be a brmsfit object.", call. = FALSE)
  }
  sf <- fit_structure$fit
  if (is.null(sf) || !methods::is(sf, "stanfit") || is.null(sf@sim)) {
    stop(
      "`x$structure$fit` is not a stanfit; cannot reconstruct.",
      call. = FALSE
    )
  }

  fnames   <- sf@sim$fnames_oi
  n_chains <- length(sf@sim$samples)

  samples <- sample_posterior(comp, n_draws = n_draws)
  missing_params <- setdiff(fnames, colnames(samples))
  if (length(missing_params)) {
    stop(
      "The reconstructed fit needs draws for parameters that are not in ",
      "the compressed posterior:\n  ",
      paste(missing_params, collapse = ", "),
      "\n  Re-run compress_brmsfit() without restricting `variables`.",
      call. = FALSE
    )
  }

  per_chain <- floor(n_draws / n_chains)
  used <- per_chain * n_chains
  samples <- samples[seq_len(used), fnames, drop = FALSE]

  for (chain in seq_len(n_chains)) {
    rows <- ((chain - 1L) * per_chain + 1L):(chain * per_chain)
    sf@sim$samples[[chain]] <- as.list(
      as.data.frame(samples[rows, , drop = FALSE])
    )
  }

  sf@sim$iter    <- per_chain
  sf@sim$warmup  <- 0L
  sf@sim$thin    <- 1L
  sf@sim$n_save  <- rep(per_chain, n_chains)
  sf@sim$warmup2 <- rep(0L, n_chains)

  fit_recon <- fit_structure
  fit_recon$fit <- sf
  attr(fit_recon, "compression_method") <- comp$method
  attr(fit_recon, "reconstructed")      <- TRUE

  fit_recon
}


#' Default number of regenerated draws inferred from a fit's stanfit shell.
#'
#' Used as the default for `n_draws` in [reconstruct_brmsfit()] and
#' [reconstruct_sccomp()] so callers don't have to supply a magic number.
#' Falls back to `4000` when no MCMC metadata is present.
#'
#' @param fit_structure A `brmsfit` (or any object with `attr(., "fit")`
#'   pointing to a `stanfit` / `CmdStanMCMC`).
#' @return An integer number of post-warmup draws across all chains.
#' @export
get_n_draws_from_fit <- function(fit_structure) {
  sf <- if (inherits(fit_structure, "brmsfit")) {
    fit_structure$fit
  } else {
    attr(fit_structure, "fit")
  }
  if (!is.null(sf) && methods::is(sf, "stanfit") && !is.null(sf@sim)) {
    n_chains <- length(sf@sim$samples)
    per_chain <- as.integer(sf@sim$iter - sf@sim$warmup)[[1L]]
    if (!is.na(per_chain) && per_chain > 0L) {
      return(n_chains * per_chain)
    }
  }
  if (!is.null(sf) && inherits(sf, "CmdStanMCMC")) {
    return(as.integer(sf$num_chains() * sf$metadata()$iter_sampling))
  }
  4000L
}


#' Reconstruct an sccomp object from a compressed posterior
#'
#' The compressed posterior and regenerated samples are stored in
#' `attr(x, "fit_compressed")`. The original `attr(x, "fit")` is removed.
#'
#' @param x The output of [compress_sccomp()].
#' @param n_draws Integer number of samples to regenerate. Defaults to
#'   the number of draws stored in the compressed posterior.
#'
#' @return The sccomp object with `attr(x, "fit_compressed")` populated.
#' @export
reconstruct_sccomp <- function(x, n_draws = x$compressed$n_draws) {
  .check_reconstruct_input(x)
  comp <- x$compressed
  obj  <- x$structure

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
