#' Compress a posterior draws matrix
#'
#' Fits a compact density approximation to a matrix of MCMC draws and
#' returns an object of class `posterior_compressed`. The compressed
#' object can be sampled from via [sample_posterior()] and evaluated via
#' [density_posterior()].
#'
#' @param draws An object containing posterior draws. Accepted:
#'   \itemize{
#'     \item a numeric matrix (rows = draws, cols = parameters);
#'     \item a `data.frame` of draws (one column per parameter);
#'     \item any object the `posterior` package can convert with
#'       `posterior::as_draws_matrix()` (`draws_matrix`, `draws_array`, ...);
#'     \item a `brms::brmsfit` or `rstan::stanfit` (via
#'       `posterior::as_draws_matrix()`; requires `posterior`);
#'     \item a `cmdstanr::CmdStanMCMC` object (requires `cmdstanr`);
#'     \item a 3-D array of shape `iter x chain x var`.
#'   }
#' @param method One of `r paste(shQuote(compression_methods()), collapse = ", ")`.
#'   See [compression_methods()].
#' @param variables Optional character vector of parameter names to keep.
#'   If `NULL` (default) all columns are used.
#' @param n_components Integer number of mixture components (used by
#'   `"mclust"` and `"mvdens_gmm"`). Default `3`.
#' @param model_name `mclust` covariance structure (e.g. `"VVV"`, `"EEE"`,
#'   or a vector of allowed model names). Ignored by other methods.
#'   When `NULL` (default) `poco` auto-selects a sensible set: the
#'   spherical and diagonal models `c("EII", "VII", "EEI", "EVI", "VEI",
#'   "VVI")` are used when `nrow(draws) <= ncol(draws)` so covariances
#'   remain identifiable, otherwise mclust's full default set is used
#'   and BIC picks the best.
#' @param verbose Logical; print blockwise progress messages. Default `FALSE`
#'   (silent). When `TRUE` and blockwise `partition` is used: prints whether
#'   cluster blocks run in parallel or sequentially, worker count
#'   (BiocParallel), enables the BiocParallel **progress bar** when supported
#'   (`progressbar` on `MulticoreParam` / `SnowParam`), prints cluster-fit
#'   completion counts in sequential mode, and reports the remainder block.
#'   Backend (`mclust` / `mvdens`) fitting is kept silent inside blockwise mode.
#' @param partition Optional [partition_parameters_clusters()] return value
#'   (default **`simple_output = TRUE`**: plain `list()` of cluster vectors),
#'   its full **`poco_partition_clusters`** form when
#'   **`simple_output = FALSE`**, or a **manual** cluster list from
#'   [partition_blocks()] / a plain `list()` of character vectors (see
#'   [partition_blocks()]; **`list()`** means compress the full posterior as
#'   a single remainder block). When supplied, the posterior
#'   is compressed **blockwise**: each cluster is fitted on its own
#'   (typically with a richer covariance family,
#'   `cluster_model_name`), the **remainder** is fitted with the cheaper
#'   diagonal `remainder_model_name`, and the per-block fits are wrapped in
#'   a `posterior_compressed_blockwise` object. Blocks are then assumed to
#'   be independent (within-block correlations preserved, between-block
#'   correlations dropped).
#'
#'   Names are aligned to the draws **before** fitting: any partition entry
#'   absent from the draws is **dropped** with [warning()], and any draw
#'   column not appearing in a cluster or in `remainder` is **added to
#'   remainder** with [warning()]. This helps when the partition was built on
#'   a different variable set than `variables` / the `brmsfit` posterior.
#' @param cluster_model_name `mclust` covariance structure to use for the
#'   *cluster* blocks when `partition` is provided. If `NULL` (default) and
#'   `method = "mclust"`, `poco` restricts BIC selection to the
#'   **non-diagonal** ellipsoidal families (`EEE`/`VEE`/`EVE`/`VVE`/`EEV`/
#'   `VEV`/`EVV`/`VVV`) so within-cluster correlations can be represented
#'   when identifiable. If that fit errors for a block (e.g. very many
#'   parameters), the same block is retried once with the diagonal/spherical
#'   family set and an informative message is shown.
#' @param remainder_model_name `mclust` covariance structure to use for the
#'   *remainder* block when `partition` is provided. If `NULL` (default) and
#'   `method = "mclust"`, `poco` restricts BIC selection to the
#'   **diagonal/spherical** families (`EII`/`VII`/`EEI`/`EVI`/`VEI`/`VVI`).
#' @param cluster_BPPARAM How to run **cluster** block fits when `partition` is
#'   set (the **remainder** block is always fitted afterward in the main
#'   process). One of:
#'   \describe{
#'     \item{`NULL` (default)}{**Automatic.** If **BiocParallel** is installed
#'       and there are at least two cluster blocks, cluster blocks are
#'       compressed in parallel via `BiocParallel::bplapply()` with worker count
#'       `min(n_cluster_blocks, max(1, parallel::detectCores() - 1))` (on
#'       Windows `SnowParam`, else `MulticoreParam`). If **BiocParallel** is not
#'       installed, a [message()] encourages installation and blocks are fitted
#'       **sequentially** (no error). With fewer than two cluster blocks,
#'       parallel is skipped.}
#'     \item{`FALSE`}{**Sequential:** always compress cluster blocks one by one,
#'       even if **BiocParallel** is installed.}
#'     \item{Single number}{**Worker count:** a single finite numeric value
#'       (e.g. `2` or `2L`). Values below `2` (including `1`, `0.5`, negatives,
#'       and non-finite values) are treated like `FALSE` (sequential), with an
#'       informative message. Non-integers use [ceiling()] (e.g. \code{2.5}
#'       becomes three workers). The count is capped by the number of cluster blocks and
#'       `parallel::detectCores() - 1`; a message is shown when capping or
#'       rounding applies. Builds `BiocParallel::MulticoreParam` (Unix-like) or
#'       `BiocParallel::SnowParam` (Windows).}
#'     \item{`BiocParallelParam`}{Use that object as `BPPARAM` in
#'       `BiocParallel::bplapply()` (requires **BiocParallel** to be installed).}
#'   }
#' @param ... Additional arguments forwarded to the backend (e.g.
#'   [mclust::Mclust()]).
#'
#' @return An S3 list of class `c("posterior_compressed_<method>",
#'   "posterior_compressed", "list")` containing the parameters of the
#'   fitted approximation. When `partition` is provided, the class is
#'   `c("posterior_compressed_blockwise", "posterior_compressed", "list")`
#'   and `$blocks` is a named list of per-block compressed objects.
#'
#' @seealso [compress_fit()], [compress_brmsfit()], [compress_sccomp()],
#'   [sample_posterior()], [density_posterior()],
#'   [partition_parameters_clusters()], [partition_blocks()].
#'
#' @examples
#' set.seed(1)
#' draws <- matrix(rnorm(2000 * 3), ncol = 3,
#'                 dimnames = list(NULL, c("alpha", "beta", "sigma")))
#' comp <- compress_posterior(draws, method = "mclust", n_components = 2)
#' comp
#' new_samples <- sample_posterior(comp, n_draws = 500)
#' dim(new_samples)
#'
#' @export
compress_posterior <- function(
    draws,
    method = c("mclust", "mvdens_gmm", "mvdens_kde"),
    variables = NULL,
    n_components = 3L,
    model_name = NULL,
    verbose = FALSE,
    partition            = NULL,
    cluster_model_name   = NULL,
    remainder_model_name = NULL,
    cluster_BPPARAM      = NULL,
    ...) {
  method <- .match_method(method)
  draws_mat <- .as_draws_matrix(draws, variables = variables)

  if (!is.null(partition)) {
    return(.compress_blockwise(
      draws_mat,
      partition            = partition,
      method               = method,
      n_components         = n_components,
      cluster_model_name   = cluster_model_name,
      remainder_model_name = remainder_model_name,
      verbose              = verbose,
      cluster_BPPARAM      = cluster_BPPARAM,
      ...
    ))
  }

  comp <- switch(
    method,
    mclust     = .compress_mclust(
      draws_mat,
      n_components = n_components,
      model_name   = model_name,
      verbose      = verbose,
      ...
    ),
    mvdens_gmm = .compress_mvdens_gmm(
      draws_mat,
      n_components = n_components,
      verbose = verbose,
      ...
    ),
    mvdens_kde = .compress_mvdens_kde(
      draws_mat,
      verbose = verbose,
      ...
    )
  )

  comp
}


#' Compress the posterior of a cmdstanr fit
#'
#' Convenience wrapper around [compress_posterior()] that accepts either a
#' `cmdstanr::CmdStanMCMC` fit or a vector of cmdstan CSV file paths.
#'
#' @param fit Either a `CmdStanMCMC` object (from `cmdstanr::sample()`) or
#'   a character vector of cmdstan CSV file paths.
#' @inheritParams compress_posterior
#' @param remove_csvs Logical; if `TRUE` and `fit` is a CmdStanMCMC,
#'   delete the CSV output files after compression. Default `FALSE`.
#'
#' @return A `posterior_compressed` object (see [compress_posterior()]).
#'
#' @seealso [compress_posterior()], [compress_brmsfit()],
#'   [reconstruct_brmsfit()].
#'
#' @examples
#' \dontrun{
#' fit <- mod$sample(data = stan_data, chains = 4)
#' comp <- compress_fit(fit, method = "mclust", n_components = 5)
#' samples <- sample_posterior(comp, n_draws = 1000)
#' }
#'
#' @export
compress_fit <- function(
    fit,
    method = c("mclust", "mvdens_gmm", "mvdens_kde"),
    variables = NULL,
    n_components = 3L,
    model_name = NULL,
    verbose = FALSE,
    remove_csvs = FALSE,
    ...) {
  method <- .match_method(method)

  if (is.character(fit)) {
    csv_files <- fit
    draws_mat <- .draws_matrix_from_csvs(csv_files, variables)
    original_size <- sum(file.size(csv_files))
  } else if (inherits(fit, "CmdStanMCMC")) {
    csv_files <- fit$output_files()
    if (length(csv_files) > 0L && all(file.exists(csv_files))) {
      draws_mat <- .draws_matrix_from_csvs(csv_files, variables)
      original_size <- sum(file.size(csv_files))
    } else {
      draws_mat <- fit$draws(variables = variables, format = "matrix")
      storage.mode(draws_mat) <- "double"
      original_size <- NA_real_
    }
  } else {
    stop(
      "`fit` must be a CmdStanMCMC object or a character vector of CSV ",
      "file paths.",
      call. = FALSE
    )
  }

  comp <- compress_posterior(
    draws_mat,
    method = method,
    n_components = n_components,
    model_name = model_name,
    verbose = verbose,
    ...
  )

  if (!is.na(original_size)) {
    comp$original_size <- original_size
    .report_compression(comp, original_size)
  }

  if (remove_csvs && length(csv_files) > 0L && all(file.exists(csv_files))) {
    unlink(csv_files)
    .message("Removed ", length(csv_files), " cmdstan CSV files.")
  }

  comp
}


#' Compress the posterior of a brms fit
#'
#' Convenience wrapper around [compress_posterior()] for `brms::brmsfit`
#' objects. Returns a list with both the compressed posterior and the
#' original fit `structure` (with `fit$fit` cleared) so the model can be
#' reconstructed later via [reconstruct_brmsfit()].
#'
#' Requires the `brms` model to have been fit with `backend = "cmdstanr"`.
#' Draws are extracted via `posterior::as_draws_matrix()`, so brms's
#' user-facing parameter names (e.g. `b_x`, `sd_group__Intercept`) are
#' preserved.
#'
#' @param brmsfit A `brms::brmsfit` object using the cmdstanr backend.
#' @param method Passed to [compress_posterior()]. See [compression_methods()].
#' @param variables Passed to [compress_posterior()]: optional subset of
#'   parameter columns after `posterior::as_draws_matrix(brmsfit)`.
#' @param n_components Passed to [compress_posterior()].
#' @param model_name Passed to [compress_posterior()] (ignored unless
#'   `method = "mclust"`).
#' @param verbose Passed to [compress_posterior()].
#' @param partition Passed to [compress_posterior()] for blockwise compression;
#'   see that function and [partition_parameters_clusters()], [partition_blocks()].
#' @param cluster_model_name Passed to [compress_posterior()] when `partition`
#'   is non-`NULL` (typically `method = "mclust"`).
#' @param remainder_model_name Passed to [compress_posterior()] when `partition`
#'   is non-`NULL` (typically `method = "mclust"`).
#' @param cluster_BPPARAM Passed to [compress_posterior()] when `partition` is
#'   non-`NULL` (parallel cluster blocks). See [compress_posterior()] for
#'   allowed values (`NULL`, `FALSE`, a worker count, or a `BiocParallelParam`).
#' @param ... Passed to [compress_posterior()] and then to the compression
#'   backend (e.g. [mclust::Mclust()]).
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{compressed}{a `posterior_compressed` object (or
#'       `posterior_compressed_blockwise` when `partition` is supplied);}
#'     \item{structure}{the brmsfit (so [reconstruct_brmsfit()] can rebuild
#'       a usable model).}
#'   }
#'
#' @seealso [compress_posterior()], [reconstruct_brmsfit()], [sample_posterior()],
#'   [partition_parameters_clusters()], [partition_blocks()].
#'
#' @examples
#' \dontrun{
#' fit <- brms::brm(y ~ x, data = dat, backend = "cmdstanr")
#' # Single-block compression (existing default).
#' result <- compress_brmsfit(fit, method = "mclust", n_components = 5)
#'
#' # Hybrid blockwise compression: rich covariance per cluster, diagonal
#' # ("VVI") for the remainder block.
#' draws <- posterior::as_draws_matrix(fit)
#' cm    <- posterior_correlation(draws)
#' part  <- partition_parameters_clusters(cm, threshold = 0.4, min_size = 0.02)
#' result_block <- compress_brmsfit(
#'   fit,
#'   method               = "mclust",
#'   n_components         = 3,
#'   partition            = part,
#'   cluster_model_name   = NULL,    # auto: ellipsoidal set; diagonal retry on error
#'   remainder_model_name = "VVI"
#' )
#'
#' saveRDS(result$compressed, "model_compressed.rds", compress = "xz")
#' saveRDS(result$structure,  "model_structure.rds")
#'
#' fit_recon <- reconstruct_brmsfit(result)
#' }
#'
#' @export
compress_brmsfit <- function(
    brmsfit,
    method = c("mclust", "mvdens_gmm", "mvdens_kde"),
    variables = NULL,
    n_components = 3L,
    model_name = NULL,
    verbose = FALSE,
    partition            = NULL,
    cluster_model_name   = NULL,
    remainder_model_name = NULL,
    cluster_BPPARAM      = NULL,
    ...) {
  method <- .match_method(method)

  if (!inherits(brmsfit, "brmsfit")) {
    stop("Input must be a brmsfit object.", call. = FALSE)
  }
  if (!identical(brmsfit$backend, "cmdstanr")) {
    stop(
      "brms fit must use the cmdstanr backend.\n",
      "  Refit with: brms::brm(..., backend = 'cmdstanr')",
      call. = FALSE
    )
  }
  if (is.null(brmsfit$fit)) {
    stop("brmsfit object has no `$fit` slot.", call. = FALSE)
  }
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop(
      "Package 'posterior' is required to compress brms fits.",
      call. = FALSE
    )
  }

  draws_mat <- as.matrix(posterior::as_draws_matrix(brmsfit))
  if (!is.null(variables)) {
    missing_vars <- setdiff(variables, colnames(draws_mat))
    if (length(missing_vars)) {
      stop(
        "The following variables are not in the brms posterior: ",
        paste(missing_vars, collapse = ", "),
        call. = FALSE
      )
    }
    draws_mat <- draws_mat[, variables, drop = FALSE]
  }
  storage.mode(draws_mat) <- "double"

  comp <- compress_posterior(
    draws_mat,
    method               = method,
    n_components         = n_components,
    model_name           = model_name,
    verbose              = verbose,
    partition            = partition,
    cluster_model_name   = cluster_model_name,
    remainder_model_name = remainder_model_name,
    cluster_BPPARAM      = cluster_BPPARAM,
    ...
  )

  structure_obj <- .strip_brmsfit_draws(brmsfit)

  .new_compressed_fit(comp, structure_obj, kind = "brmsfit")
}


#' Strip the heavy stanfit draws from a brmsfit, keeping the shell so
#' [reconstruct_brmsfit()] can refill `@sim$samples` later.
#'
#' We also prepend the class `"brmsfit_stripped"` so that auto-printing
#' the structure (e.g. `result$structure` at the console) dispatches to
#' our small [print.brmsfit_stripped()] method instead of
#' `brms::print.brmsfit`, which would call into the rstan summary
#' machinery and fail with `do.call(cbind, attr(x, "sampler_params"))`
#' on the empty stanfit. `inherits(., "brmsfit")` still holds.
#' @keywords internal
#' @noRd
.strip_brmsfit_draws <- function(brmsfit) {
  sf <- brmsfit$fit
  if (is.null(sf) || !methods::is(sf, "stanfit") || is.null(sf@sim)) {
    return(.tag_brmsfit_stripped(brmsfit))
  }
  for (chain in seq_along(sf@sim$samples)) {
    sf@sim$samples[[chain]] <- lapply(sf@sim$samples[[chain]], function(x) {
      numeric(0L)
    })
  }
  sf@sim$n_save  <- rep(0L, length(sf@sim$samples))
  sf@sim$warmup2 <- rep(0L, length(sf@sim$samples))
  brmsfit$fit <- sf
  .tag_brmsfit_stripped(brmsfit)
}


#' Add/remove the `"brmsfit_stripped"` class tag.
#' @keywords internal
#' @noRd
.tag_brmsfit_stripped <- function(brmsfit) {
  cls <- class(brmsfit)
  if (!"brmsfit_stripped" %in% cls) {
    class(brmsfit) <- c("brmsfit_stripped", cls)
  }
  brmsfit
}

#' @keywords internal
#' @noRd
.untag_brmsfit_stripped <- function(brmsfit) {
  cls <- class(brmsfit)
  cls <- cls[cls != "brmsfit_stripped"]
  class(brmsfit) <- cls
  brmsfit
}


#' Print method for a draws-stripped brmsfit shell
#'
#' Avoids dispatching to `brms::print.brmsfit`, whose summary machinery
#' fails on a stanfit whose `@sim$samples` have been zeroed out by
#' [compress_brmsfit()] (typical error:
#' `do.call(cbind, attr(x, "sampler_params")) : second argument must be a list`).
#'
#' @param x A brmsfit returned inside `compress_brmsfit()$structure`.
#' @param ... Currently unused.
#' @return Invisibly returns `x`.
#' @export
print.brmsfit_stripped <- function(x, ...) {
  cls <- setdiff(class(x), "brmsfit_stripped")
  cat("<brmsfit, draws stripped by poco::compress_brmsfit()>\n")
  cat("  classes : ", paste(cls, collapse = "/"), "\n", sep = "")
  if (!is.null(x$formula)) {
    f <- tryCatch(format(x$formula), error = function(e) NULL)
    if (!is.null(f)) {
      cat("  formula : ", paste(f, collapse = " "), "\n", sep = "")
    }
  }
  if (!is.null(x$family)) {
    fam <- tryCatch(x$family$family, error = function(e) NULL)
    if (!is.null(fam)) cat("  family  : ", fam, "\n", sep = "")
  }
  sf <- x$fit
  if (!is.null(sf) && methods::is(sf, "stanfit") && !is.null(sf@sim)) {
    cat(
      "  chains  : ", length(sf@sim$samples),
      "  (per-chain draws stored: 0)\n",
      sep = ""
    )
  }
  cat("  Use reconstruct_brmsfit() to refill draws from the compressed posterior.\n")
  invisible(x)
}


#' Compress the posterior of an sccomp fit
#'
#' Convenience wrapper around [compress_posterior()] for `sccomp` fits.
#' `sccomp` stores the underlying `cmdstanr` fit in `attr(x, "fit")` rather
#' than `x$fit`.
#'
#' @param sccomp_obj An sccomp fit object (with `attr(x, "fit")` a
#'   `CmdStanMCMC`).
#' @inheritParams compress_posterior
#' @param remove_csvs Logical; if `TRUE`, delete the cmdstan CSV files
#'   after compression. Default `FALSE`.
#'
#' @return A list with `compressed` and `structure` (the sccomp object
#'   with `attr(x, "fit")` cleared). Use [reconstruct_sccomp()] to rebuild.
#'
#' @seealso [reconstruct_sccomp()].
#'
#' @export
compress_sccomp <- function(
    sccomp_obj,
    method = c("mclust", "mvdens_gmm", "mvdens_kde"),
    variables = NULL,
    n_components = 3L,
    model_name = NULL,
    verbose = FALSE,
    remove_csvs = FALSE,
    ...) {
  method <- .match_method(method)

  fit <- attr(sccomp_obj, "fit")
  if (is.null(fit)) {
    stop(
      "sccomp object must have a cmdstanr fit in attr(x, 'fit').",
      call. = FALSE
    )
  }
  if (!inherits(fit, "CmdStanMCMC")) {
    stop(
      "sccomp fit must use the cmdstanr backend (attr(x, 'fit') is a ",
      "CmdStanMCMC).",
      call. = FALSE
    )
  }

  comp <- compress_fit(
    fit,
    method = method,
    variables = variables,
    n_components = n_components,
    model_name = model_name,
    verbose = verbose,
    remove_csvs = remove_csvs,
    ...
  )

  structure_obj <- sccomp_obj
  attr(structure_obj, "fit") <- NULL

  .new_compressed_fit(comp, structure_obj, kind = "sccomp")
}


#' Build a classed wrapper around a compressed posterior + fit structure
#'
#' Returns a list with `compressed` and `structure`, classed as
#' `c("compressed_<kind>", "compressed_fit", "list")` so [print()] does
#' not dispatch to e.g. `print.brmsfit` on the draws-stripped structure
#' (which can fail).
#' @keywords internal
#' @noRd
.new_compressed_fit <- function(compressed, structure_obj,
                                kind = c("brmsfit", "sccomp")) {
  kind <- match.arg(kind)
  out <- list(compressed = compressed, structure = structure_obj)
  class(out) <- c(paste0("compressed_", kind), "compressed_fit", "list")
  out
}


#' Print method for `compressed_fit` wrappers
#'
#' Bypasses the auto-print of the (draws-stripped) `structure` element to
#' keep `brms::print.brmsfit` from running on an empty stanfit.
#'
#' @param x A `compressed_fit` object from [compress_brmsfit()] or
#'   [compress_sccomp()].
#' @param ... Currently unused.
#' @return Invisibly returns `x`.
#' @export
print.compressed_fit <- function(x, ...) {
  kind <- if (inherits(x, "compressed_brmsfit")) {
    "brmsfit"
  } else if (inherits(x, "compressed_sccomp")) {
    "sccomp"
  } else {
    "fit"
  }
  cat("<compressed_", kind, ">\n", sep = "")
  cat("$compressed\n")
  print(x$compressed)
  cat("\n$structure\n")
  cls <- class(x$structure)
  recon <- if (kind == "brmsfit") "reconstruct_brmsfit()"
           else if (kind == "sccomp") "reconstruct_sccomp()"
           else NULL
  cat(
    "  <", paste(cls, collapse = "/"), "> ",
    "(draws stripped",
    if (!is.null(recon)) paste0(", use ", recon, " to rebuild") else "",
    ")\n",
    sep = ""
  )
  invisible(x)
}


#' @describeIn print.compressed_fit S3 method for [compress_brmsfit()] results.
#' @export
print.compressed_brmsfit <- function(x, ...) {
  print.compressed_fit(x, ...)
}


#' @describeIn print.compressed_fit S3 method for [compress_sccomp()] results.
#' @export
print.compressed_sccomp <- function(x, ...) {
  print.compressed_fit(x, ...)
}


#' @keywords internal
#' @noRd
.report_compression <- function(comp, original_size) {
  comp_size <- as.numeric(utils::object.size(comp))
  ratio <- original_size / comp_size
  .message(
    "Compression: ", .size_pretty(original_size),
    " B raw -> ~", .size_pretty(round(comp_size)),
    " B in memory (", round(ratio), "x)."
  )
  invisible(NULL)
}
