#' CmdStan helpers
#'
#' `poco` itself does not ship Stan models, but several of its
#' wrappers and vignettes interact with `cmdstanr` fits (and so with the
#' underlying CmdStan install). The helpers below mirror the pattern used by
#' the `sccomp` package, which routes all CmdStan availability checks through
#' the lightweight `instantiate` package.
#'
#' @name cmdstan_utils
#' @keywords internal
NULL


#' Default cache directory for compiled Stan models
#'
#' Mirrors `sccomp::sccomp_stan_models_cache_dir`. Used as the default cache
#' location when `poco` (or downstream packages built on top of it) compile or
#' load Stan models via `cmdstanr` / `instantiate`.
#'
#' @return A character string with the path to the cache directory.
#' @export
#' @examples
#' poco_stan_models_cache_dir
poco_stan_models_cache_dir <- file.path(path.expand("~"), ".poco_models")


#' Check (and optionally install) `cmdstanr` and CmdStan
#'
#' Light wrapper around the same logic used by `sccomp`. Verifies that
#' `cmdstanr` (>= 0.9.0) and a working CmdStan install are available; if not,
#' it offers an actionable message (and, in non-interactive sessions, attempts
#' to install `cmdstanr` from the Stan-dev R-Universe).
#'
#' @return Invisibly `TRUE` if everything is available, `FALSE` otherwise.
#' @export
#' @importFrom instantiate stan_cmdstan_exists
#' @examples
#' \dontrun{
#' check_and_install_cmdstanr()
#' }
check_and_install_cmdstanr <- function() {
  ok_pkg <- requireNamespace("cmdstanr", quietly = TRUE) &&
    utils::packageVersion("cmdstanr") >= "0.9.0"

  if (!ok_pkg) {
    if (requireNamespace("rlang", quietly = TRUE)) {
      rlang::check_installed(
        pkg     = "cmdstanr",
        version = "0.9.0",
        reason  = paste(
          "{cmdstanr} (>= 0.9.0) is required to fit / read Stan models",
          "used by poco's compress_fit() and compress_brmsfit()."
        ),
        action  = function(...) utils::install.packages(
          ...,
          repos = c(
            "https://stan-dev.r-universe.dev",
            "https://cloud.r-project.org"
          )
        )
      )
    } else {
      stop(
        "Package 'cmdstanr' (>= 0.9.0) is required.\n",
        "  install.packages('cmdstanr',\n",
        "    repos = c('https://stan-dev.r-universe.dev/', getOption('repos')))",
        call. = FALSE
      )
    }
  }

  if (!instantiate::stan_cmdstan_exists()) {
    message(
      "CmdStan is not installed.\n",
      "  Run the following to install it:\n",
      "    cmdstanr::check_cmdstan_toolchain(fix = TRUE)\n",
      "    cmdstanr::install_cmdstan()\n",
      "  See https://mc-stan.org/users/interfaces/cmdstan for details."
    )
    return(invisible(FALSE))
  }

  invisible(TRUE)
}


#' Has CmdStan been installed?
#'
#' Re-exports `instantiate::stan_cmdstan_exists()` so users (and vignettes)
#' have a single, package-friendly entry point to ask "is the CmdStan
#' toolchain available on this machine?". Useful for gating
#' `cmdstanr` / `brms` examples in CI environments where CmdStan is not
#' available.
#'
#' @return Logical scalar.
#' @export
#' @importFrom instantiate stan_cmdstan_exists
#' @examples
#' has_cmdstan()
has_cmdstan <- function() {
  isTRUE(instantiate::stan_cmdstan_exists())
}


#' Clear the Stan model cache
#'
#' Removes the cache directory used to store compiled Stan models (see
#' [poco_stan_models_cache_dir]). Useful when a cached model is stale or
#' fails to load. Mirrors `sccomp::clear_stan_model_cache()`.
#'
#' @param cache_dir Path to the cache directory to remove. Defaults to
#'   [poco_stan_models_cache_dir].
#'
#' @return Invisibly `NULL`.
#' @export
#' @examples
#' \dontrun{
#' clear_stan_model_cache()
#' }
clear_stan_model_cache <- function(cache_dir = poco_stan_models_cache_dir) {
  if (dir.exists(cache_dir)) {
    unlink(cache_dir, recursive = TRUE)
    message("Cache deleted: ", cache_dir)
  } else {
    message("Cache does not exist: ", cache_dir)
  }
  invisible(NULL)
}
