#' Reconstruct a brms fit object from a compressed posterior
#'
#' Replaces the posterior draws of a `brmsfit` with samples regenerated
#' from a compressed posterior, so the standard `brms` post-processing
#' methods (`predict()`, `posterior_predict()`, `pp_check()`,
#' `conditional_effects()`, ...) keep working with the lighter object.
#'
#' Works with both single-block compressions
#' (`posterior_compressed_mclust`, `posterior_compressed_mvdens_*`) and
#' the hybrid blockwise variant (`posterior_compressed_blockwise`).
#' For blockwise objects, [sample_posterior()] regenerates each block
#' independently and re-aligns columns to the original parameter order
#' before the per-chain refill below — no special handling is required.
#'
#' MCMC diagnostics (`Rhat`, `ESS`, ...) are not preserved.
#'
#' @param x The output of [compress_brmsfit()] (a list with `compressed`
#'   and `structure`), or a similarly-shaped list constructed from
#'   separately-saved objects.
#' @param n_draws Integer number of samples to regenerate. Defaults to
#'   the original number of draws stored in `x$structure`.
#'
#' @return A `brmsfit` object with regenerated draws and the following
#'   informational attributes:
#'   \describe{
#'     \item{`compression_method`}{`comp$method` (e.g. `"mclust"` or
#'       `"blockwise"`);}
#'     \item{`compression_base_method`}{for blockwise compressions, the
#'       per-block backend (e.g. `"mclust"`); `NA_character_` otherwise;}
#'     \item{`compression_blocks`}{for blockwise compressions, an integer
#'       vector of per-block parameter counts (named after the block);}
#'     \item{`reconstructed`}{`TRUE`.}
#'   }
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
      paste(utils::head(missing_params, 5L), collapse = ", "),
      if (length(missing_params) > 5L) ", ..." else "",
      "\n  Re-run compress_brmsfit() without restricting `variables`",
      if (inherits(comp, "posterior_compressed_blockwise")) {
        " (or check the partition covers every parameter)."
      } else {
        "."
      },
      call. = FALSE
    )
  }

  per_chain <- floor(n_draws / n_chains)
  used <- per_chain * n_chains
  samples <- samples[seq_len(used), fnames, drop = FALSE]

  # Build placeholder sampler-param *lists* once; rstan / brms read this
  # both as a stanfit-level attribute (list of matrices, one per chain)
  # AND per-chain as `attr(sim$samples[[i]], "sampler_params")`. Without
  # the stanfit-level entry, auto-printing the reconstructed brmsfit dies
  # with: do.call(cbind, attr(x, "sampler_params")) :
  #         second argument must be a list.
  placeholder_sp <- replicate(
    n_chains,
    .placeholder_sampler_params(per_chain),
    simplify = FALSE
  )

  for (chain in seq_len(n_chains)) {
    rows <- ((chain - 1L) * per_chain + 1L):(chain * per_chain)
    chain_samples <- as.list(
      as.data.frame(samples[rows, , drop = FALSE])
    )
    attr(chain_samples, "sampler_params") <- placeholder_sp[[chain]]
    sf@sim$samples[[chain]] <- chain_samples
  }

  # Stanfit-level attribute (list of per-chain matrices). Some rstan code
  # paths read this directly via `attr(stanfit, "sampler_params")`, e.g.
  # `do.call(cbind, attr(stanfit, "sampler_params"))`.
  attr(sf, "sampler_params") <- lapply(
    placeholder_sp,
    function(sp) do.call(cbind, sp)
  )

  # Stanfit-level `elapsed_time` is also commonly read by print/summary;
  # set a benign zero matrix per chain so `get_elapsed_time()` does not
  # crash on a stripped fit. (`do.call(rbind, NULL)` errors similarly.)
  if (is.null(attr(sf, "elapsed_time"))) {
    attr(sf, "elapsed_time") <- replicate(
      n_chains,
      stats::setNames(c(0, 0), c("warmup", "sample")),
      simplify = FALSE
    )
  }

  sf@sim$iter    <- per_chain
  sf@sim$warmup  <- 0L
  sf@sim$thin    <- 1L
  sf@sim$n_save  <- rep(per_chain, n_chains)
  sf@sim$warmup2 <- rep(0L, n_chains)

  fit_recon <- fit_structure
  fit_recon$fit <- sf
  cls <- class(fit_recon)
  class(fit_recon) <- cls[cls != "brmsfit_stripped"]

  attr(fit_recon, "compression_method") <- comp$method
  attr(fit_recon, "reconstructed")      <- TRUE
  if (inherits(comp, "posterior_compressed_blockwise")) {
    attr(fit_recon, "compression_base_method") <- comp$base_method
    attr(fit_recon, "compression_blocks") <- vapply(
      comp$blocks,
      function(b) as.integer(b$n_params %||% length(b$param_names)),
      integer(1L)
    )
  } else {
    attr(fit_recon, "compression_base_method") <- NA_character_
  }

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


#' Standard rstan sampler-param placeholders used by HMC/NUTS chains.
#' `reconstruct_brmsfit()` attaches a finite-valued *named list* of vectors
#' (one vector per sampler-param column), because rstan internally calls
#' `do.call(cbind, attr(chain, "sampler_params"))` and therefore expects
#' this attribute to be list-like. Values are benign defaults (e.g.
#' `divergent__ = 0`) so brms/rstan summary code that uses logical checks
#' like `if (div_trans > 0)` does not fail on `NA`.
#' @keywords internal
#' @noRd
.placeholder_sampler_params <- function(n_iter) {
  list(
    accept_stat__ = rep(1.0, n_iter),
    stepsize__    = rep(1.0, n_iter),
    treedepth__   = rep(1.0, n_iter),
    n_leapfrog__  = rep(1.0, n_iter),
    divergent__   = rep(0.0, n_iter),
    energy__      = rep(0.0, n_iter)
  )
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
