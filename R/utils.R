#' Coerce an input into a numeric draws matrix
#'
#' Internal helper used by all `compress_*()` functions to normalise input
#' into a matrix of draws (rows = draws, cols = parameters).
#'
#' Supported inputs:
#' * matrix: returned as-is (with a sanity check)
#' * data.frame: coerced via [as.matrix()]
#' * `posterior::draws_*` objects: converted via
#'   `posterior::as_draws_matrix()` (requires the `posterior` package)
#' * `brms::brmsfit` and `rstan::stanfit`: draws via
#'   `posterior::as_draws_matrix()` (requires `posterior`; `brms` for
#'   `brmsfit` conversion)
#' * `cmdstanr::CmdStanMCMC`: `fit$draws(..., format = "matrix")` (requires
#'   `cmdstanr`)
#' * any other object `posterior::as_draws_matrix()` accepts, if `posterior`
#'   is installed (attempted last)
#' * 3-D array (iter x chain x var): collapsed across iter/chain
#'
#' @param draws An object containing posterior draws.
#' @param variables Optional character vector of parameters to keep.
#'
#' @return A numeric matrix with column names.
#' @keywords internal
#' @noRd
.as_draws_matrix <- function(draws, variables = NULL) {
  mat <- NULL

  if (is.matrix(draws) && is.numeric(draws)) {
    mat <- draws
  } else if (is.data.frame(draws)) {
    mat <- as.matrix(draws)
  } else if (inherits(draws, c("draws_matrix", "draws_array",
                               "draws_df", "draws_list", "draws_rvars"))) {
    if (!requireNamespace("posterior", quietly = TRUE)) {
      stop(
        "Package 'posterior' is required to convert objects of class '",
        paste(class(draws), collapse = "/"), "'.",
        " Install with install.packages('posterior').",
        call. = FALSE
      )
    }
    mat <- as.matrix(posterior::as_draws_matrix(draws))
  } else if (inherits(draws, "brmsfit") || inherits(draws, "stanfit")) {
    if (!requireNamespace("posterior", quietly = TRUE)) {
      stop(
        "Package 'posterior' is required to extract draws from objects of class '",
        paste(class(draws), collapse = "/"), "'.\n",
        "  install.packages('posterior')",
        call. = FALSE
      )
    }
    mat <- as.matrix(posterior::as_draws_matrix(draws))
  } else if (inherits(draws, "CmdStanMCMC")) {
    if (!requireNamespace("cmdstanr", quietly = TRUE)) {
      stop(
        "Package 'cmdstanr' is required for CmdStanMCMC objects.\n",
        "  https://mc-stan.org/cmdstanr",
        call. = FALSE
      )
    }
    mat <- draws$draws(variables = variables, format = "matrix")
    storage.mode(mat) <- "double"
  } else if (is.array(draws) && length(dim(draws)) == 3L) {
    d <- dim(draws)
    nm <- dimnames(draws)[[3]]
    mat <- matrix(
      aperm(draws, c(1, 2, 3)),
      nrow = d[1] * d[2],
      ncol = d[3]
    )
    colnames(mat) <- nm
  } else if (requireNamespace("posterior", quietly = TRUE)) {
    mat <- tryCatch(
      as.matrix(posterior::as_draws_matrix(draws)),
      error = function(e) NULL
    )
  }

  if (is.null(mat)) {
    stop(
      "Don't know how to coerce object of class '",
      paste(class(draws), collapse = "/"),
      "' to a posterior draws matrix.",
      call. = FALSE
    )
  }

  if (is.null(colnames(mat))) {
    colnames(mat) <- paste0("V", seq_len(ncol(mat)))
  }

  if (!is.null(variables)) {
    missing_vars <- setdiff(variables, colnames(mat))
    if (length(missing_vars) > 0) {
      stop(
        "The following requested variables were not found in draws: ",
        paste(missing_vars, collapse = ", "),
        call. = FALSE
      )
    }
    mat <- mat[, variables, drop = FALSE]
  }

  storage.mode(mat) <- "double"
  mat
}


#' Read a draws matrix from cmdstan CSV files
#'
#' @param csv_files Character vector of cmdstan CSV file paths.
#' @param variables Optional character vector of parameter names.
#'
#' @return A numeric draws matrix.
#' @keywords internal
#' @noRd
.draws_matrix_from_csvs <- function(csv_files, variables = NULL) {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop(
      "Package 'cmdstanr' is required to read Stan CSV files.\n",
      "  Install from https://mc-stan.org/cmdstanr",
      call. = FALSE
    )
  }
  if (length(csv_files) == 0L || !all(file.exists(csv_files))) {
    stop(
      "All `csv_files` must exist; ",
      sum(!file.exists(csv_files)),
      " missing.",
      call. = FALSE
    )
  }
  fit <- cmdstanr::as_cmdstan_fit(csv_files)
  draws <- fit$draws(variables = variables, format = "matrix")
  storage.mode(draws) <- "double"
  draws
}


#' List of all supported compression methods
#'
#' Returns the methods accepted by the `method` argument of
#' [compress_posterior()] and friends.
#'
#' @return A character vector of method identifiers.
#' @export
#' @examples
#' compression_methods()
compression_methods <- function() {
  c("mclust", "mvdens_gmm", "mvdens_kde")
}


#' @keywords internal
#' @noRd
.match_method <- function(method) {
  method <- match.arg(method, compression_methods())
  if (method %in% c("mvdens_gmm", "mvdens_kde") &&
      !requireNamespace("mvdens", quietly = TRUE)) {
    stop(
      "Method '", method, "' requires the 'mvdens' package, which is not ",
      "installed.\n",
      "  remotes::install_github('NKI-CCB/mvdens')",
      call. = FALSE
    )
  }
  method
}


#' @keywords internal
#' @noRd
.message <- function(...) {
  if (requireNamespace("cli", quietly = TRUE)) {
    cli::cli_alert_info(paste0(...))
  } else {
    message(...)
  }
}


#' @keywords internal
#' @noRd
.size_pretty <- function(x) {
  format(x, big.mark = ",", scientific = FALSE)
}
