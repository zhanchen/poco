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
#'     \item a 3-D array of shape `iter x chain x var`.
#'   }
#' @param method One of `r paste(shQuote(compression_methods()), collapse = ", ")`.
#'   See [compression_methods()].
#' @param variables Optional character vector of parameter names to keep.
#'   If `NULL` (default) all columns are used.
#' @param n_components Integer number of mixture components (used by
#'   `"mclust"` and `"mvdens_gmm"`). Default `3`.
#' @param model_name `mclust` covariance structure (e.g. `"VVV"`, `"EEE"`).
#'   Ignored by other methods. Default `"VVV"`.
#' @param verbose Logical; print backend progress. Default `FALSE`.
#' @param ... Additional arguments forwarded to the backend (e.g.
#'   [mclust::Mclust()]).
#'
#' @return An S3 list of class `c("posterior_compressed_<method>",
#'   "posterior_compressed", "list")` containing the parameters of the
#'   fitted approximation.
#'
#' @seealso [compress_fit()], [compress_brmsfit()], [compress_sccomp()],
#'   [sample_posterior()], [density_posterior()].
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
    model_name = "VVV",
    verbose = FALSE,
    ...) {
  method <- .match_method(method)
  draws_mat <- .as_draws_matrix(draws, variables = variables)

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
    model_name = "VVV",
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
#'
#' @param brmsfit A `brms::brmsfit` object using the cmdstanr backend.
#' @inheritParams compress_posterior
#' @param remove_csvs Logical; if `TRUE`, delete the cmdstan CSV files
#'   after compression. Default `FALSE`.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{compressed}{a `posterior_compressed` object;}
#'     \item{structure}{the brmsfit (so [reconstruct_brmsfit()] can rebuild
#'       a usable model).}
#'   }
#'
#' @seealso [reconstruct_brmsfit()], [sample_posterior()].
#'
#' @examples
#' \dontrun{
#' fit <- brms::brm(y ~ x, data = dat, backend = "cmdstanr")
#' result <- compress_brmsfit(fit, method = "mclust", n_components = 5)
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
    model_name = "VVV",
    verbose = FALSE,
    remove_csvs = FALSE,
    ...) {
  method <- .match_method(method)

  if (!inherits(brmsfit, "brmsfit")) {
    stop("Input must be a brmsfit object.", call. = FALSE)
  }
  if (is.null(brmsfit$fit) || !inherits(brmsfit$fit, "CmdStanMCMC")) {
    stop(
      "brms fit must use the cmdstanr backend.\n",
      "  Refit with: brms::brm(..., backend = 'cmdstanr')",
      call. = FALSE
    )
  }

  if (is.null(variables)) {
    if (requireNamespace("posterior", quietly = TRUE)) {
      variables <- posterior::variables(brmsfit)
    } else {
      variables <- colnames(brmsfit$fit$draws(format = "matrix"))
    }
  }

  comp <- compress_fit(
    brmsfit$fit,
    method = method,
    variables = variables,
    n_components = n_components,
    model_name = model_name,
    verbose = verbose,
    remove_csvs = remove_csvs,
    ...
  )

  list(compressed = comp, structure = brmsfit)
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
    model_name = "VVV",
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

  list(compressed = comp, structure = structure_obj)
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
