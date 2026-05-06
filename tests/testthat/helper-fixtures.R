make_two_blob_draws <- function(n_per = 1000L, seed = 1L) {
  set.seed(seed)
  m1 <- mvtnorm::rmvnorm(n_per, mean = c(-2, -1), sigma = diag(c(0.5, 0.3)))
  m2 <- mvtnorm::rmvnorm(n_per, mean = c( 2,  1), sigma = diag(c(0.4, 0.6)))
  draws <- rbind(m1, m2)
  colnames(draws) <- c("alpha", "beta")
  draws
}

skip_if_no_brms_cmdstanr <- function(fit) {
  testthat::skip_if_not_installed("brms")
  testthat::skip_if_not_installed("cmdstanr")
  if (is.null(fit$fit) || !inherits(fit$fit, "CmdStanMCMC")) {
    testthat::skip("brms fit did not use cmdstanr backend")
  }
  invisible(TRUE)
}
