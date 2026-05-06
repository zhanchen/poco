test_that("compress_brmsfit input validation", {
  expect_error(compress_brmsfit("not a brmsfit"), "brmsfit")
  expect_error(compress_brmsfit(NULL), "brmsfit")

  mock <- structure(list(), class = "brmsfit")
  expect_error(compress_brmsfit(mock), "cmdstanr backend")
})

test_that("reconstruct_brmsfit input validation", {
  expect_error(reconstruct_brmsfit("not a list"), "compressed.*structure")
  expect_error(
    reconstruct_brmsfit(list(compressed = NULL, structure = NULL)),
    "posterior_compressed"
  )
})

test_that("compress_brmsfit + reconstruct_brmsfit (integration)", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("posterior")

  tryCatch(cmdstanr::cmdstan_path(), error = function(e) skip("cmdstan not installed"))

  set.seed(1)
  dat <- data.frame(x = rnorm(40), y = NA_real_)
  dat$y <- 1 + 0.5 * dat$x + rnorm(40)

  fit <- suppressMessages(suppressWarnings(brms::brm(
    y ~ x,
    data = dat,
    chains = 2,
    iter = 600,
    warmup = 300,
    backend = "cmdstanr",
    refresh = 0,
    silent = 2
  )))
  skip_if_no_brms_cmdstanr(fit)

  res <- compress_brmsfit(fit, method = "mclust", n_components = 2,
                          remove_csvs = FALSE)
  expect_named(res, c("compressed", "structure"))
  expect_s3_class(res$compressed, "posterior_compressed_mclust")
  expect_s3_class(res$structure, "brmsfit")

  recon <- reconstruct_brmsfit(res, n_draws = 400)
  expect_s3_class(recon, "brmsfit")
  expect_true(!is.null(attr(recon, "regenerated_draws")))
  expect_equal(attr(recon, "compression_method"), "mclust")
})
