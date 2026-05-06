test_that("compression_methods() returns the expected set", {
  expect_setequal(
    compression_methods(),
    c("mclust", "mvdens_gmm", "mvdens_kde")
  )
})

test_that("compress_posterior fits an mclust GMM to a draws matrix", {
  draws <- make_two_blob_draws()
  comp <- compress_posterior(draws, method = "mclust", n_components = 2)

  expect_s3_class(comp, "posterior_compressed")
  expect_s3_class(comp, "posterior_compressed_mclust")
  expect_equal(comp$method, "mclust")
  expect_equal(comp$n_params, 2L)
  expect_equal(comp$n_components, 2L)
  expect_equal(comp$param_names, c("alpha", "beta"))
  expect_true(is.numeric(comp$weights))
  expect_equal(length(comp$weights), 2L)
  expect_equal(sum(comp$weights), 1, tolerance = 1e-6)
})

test_that("compress_posterior accepts a data.frame and selects variables", {
  draws <- as.data.frame(make_two_blob_draws())
  comp <- compress_posterior(
    draws,
    method = "mclust",
    variables = "alpha",
    n_components = 2
  )
  expect_equal(comp$param_names, "alpha")
  expect_equal(comp$n_params, 1L)
})

test_that("compress_posterior errors on unknown variables", {
  draws <- make_two_blob_draws()
  expect_error(
    compress_posterior(draws, method = "mclust", variables = "missing"),
    "missing"
  )
})

test_that("compress_posterior errors on unknown methods", {
  draws <- make_two_blob_draws()
  expect_error(
    compress_posterior(draws, method = "not_a_method"),
    "should be one of"
  )
})

test_that("sample_posterior returns the right shape and names", {
  draws <- make_two_blob_draws()
  comp <- compress_posterior(draws, method = "mclust", n_components = 2)

  s_default <- sample_posterior(comp)
  expect_equal(nrow(s_default), comp$n_draws)
  expect_equal(colnames(s_default), comp$param_names)

  s_custom <- sample_posterior(comp, n_draws = 500)
  expect_equal(nrow(s_custom), 500L)
  expect_equal(ncol(s_custom), 2L)
  expect_true(all(is.finite(s_custom)))
})

test_that("density_posterior returns non-negative values", {
  draws <- make_two_blob_draws()
  comp <- compress_posterior(draws, method = "mclust", n_components = 2)

  pts <- rbind(c(-2, -1), c(2, 1), c(0, 0))
  d <- density_posterior(comp, pts)
  expect_length(d, 3L)
  expect_true(all(d >= 0))
  expect_true(all(is.finite(d)))

  ld <- density_posterior(comp, pts, log = TRUE)
  expect_equal(d, exp(ld), tolerance = 1e-8)
})

test_that("density_posterior errors on wrong column count", {
  draws <- make_two_blob_draws()
  comp <- compress_posterior(draws, method = "mclust", n_components = 2)
  expect_error(
    density_posterior(comp, matrix(0, ncol = 3, nrow = 1)),
    "3 columns"
  )
})

test_that("sample_posterior accepts an .rds path", {
  draws <- make_two_blob_draws()
  comp <- compress_posterior(draws, method = "mclust", n_components = 2)
  tf <- tempfile(fileext = ".rds")
  on.exit(unlink(tf), add = TRUE)
  saveRDS(comp, tf)

  s <- sample_posterior(tf, n_draws = 100)
  expect_equal(nrow(s), 100L)
  expect_equal(colnames(s), comp$param_names)
})

test_that("print and summary produce sensible output", {
  draws <- make_two_blob_draws()
  comp <- compress_posterior(draws, method = "mclust", n_components = 2)

  expect_output(print(comp), "posterior_compressed")
  s <- summary(comp)
  expect_s3_class(s, "summary.posterior_compressed")
  expect_output(print(s), "Compressed posterior")
})

test_that("default S3 methods error helpfully on wrong class", {
  expect_error(sample_posterior(list()), "No sample_posterior")
  expect_error(density_posterior(list(), matrix(1)), "No density_posterior")
})
