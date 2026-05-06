make_two_blob_draws <- function(n_per = 1000L, seed = 1L) {
  set.seed(seed)
  m1 <- mvtnorm::rmvnorm(n_per, mean = c(-2, -1), sigma = diag(c(0.5, 0.3)))
  m2 <- mvtnorm::rmvnorm(n_per, mean = c( 2,  1), sigma = diag(c(0.4, 0.6)))
  draws <- rbind(m1, m2)
  colnames(draws) <- c("alpha", "beta")
  draws
}

ensure_brms_cmdstanr <- function(fit) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Package 'brms' must be installed for this test.", call. = FALSE)
  }
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("Package 'cmdstanr' must be installed for this test.", call. = FALSE)
  }
  if (!identical(fit$backend, "cmdstanr")) {
    stop("brms fit did not use cmdstanr backend", call. = FALSE)
  }
  invisible(TRUE)
}
