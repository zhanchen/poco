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

test_that("univariate mclust compression samples without isSymmetric error", {
  set.seed(1L)
  n <- 300L
  draws <- matrix(rnorm(n), ncol = 1L, dimnames = list(NULL, "x"))
  comp <- compress_posterior(draws, method = "mclust", n_components = 2L)
  expect_equal(comp$covariance_type, "diagonal")
  expect_equal(dim(comp$covariances), c(1L, comp$n_components))
  s <- sample_posterior(comp, n_draws = 100L)
  expect_equal(dim(s), c(100L, 1L))
  expect_true(all(is.finite(s)))
  d <- density_posterior(comp, matrix(c(0, 0.5, -1), ncol = 1L))
  expect_true(all(is.finite(d)) && all(d >= 0))
})

test_that("univariate mclust infers diagonal layout if covariance_type was stripped", {
  set.seed(3L)
  comp <- compress_posterior(
    matrix(rnorm(250L), ncol = 1L, dimnames = list(NULL, "x")),
    method = "mclust",
    n_components = 2L
  )
  comp$covariance_type <- NULL
  expect_no_error(s <- sample_posterior(comp, n_draws = 40L))
  expect_equal(ncol(s), 1L)
  expect_true(all(is.finite(s)))
  expect_no_error(density_posterior(comp, matrix(0.1, nrow = 1L, ncol = 1L)))
})

test_that("axis-aligned mclust models are stored as diagonal covariances", {
  draws <- make_two_blob_draws()
  comp <- compress_posterior(
    draws,
    method = "mclust",
    n_components = 2,
    model_name = "VVI"
  )

  expect_equal(comp$covariance_type, "diagonal")
  expect_true(is.matrix(comp$covariances))
  expect_equal(dim(comp$covariances), c(2L, comp$n_components))

  set.seed(42)
  s <- sample_posterior(comp, n_draws = 500)
  expect_equal(dim(s), c(500L, 2L))
  expect_equal(colnames(s), comp$param_names)
  expect_true(all(is.finite(s)))

  pts <- rbind(c(-2, -1), c(2, 1), c(0, 0))
  d <- density_posterior(comp, pts)
  expect_length(d, 3L)
  expect_true(all(is.finite(d)) && all(d >= 0))

  expect_equal(
    density_posterior(comp, pts),
    exp(density_posterior(comp, pts, log = TRUE)),
    tolerance = 1e-8
  )
})

test_that("diagonal storage matches a manually-built full covariance", {
  draws <- make_two_blob_draws()
  comp_diag <- compress_posterior(
    draws,
    method = "mclust",
    n_components = 2,
    model_name = "VVI"
  )

  d <- comp_diag$n_params
  G <- comp_diag$n_components
  full_arr <- array(0, dim = c(d, d, G))
  for (k in seq_len(G)) {
    full_arr[, , k] <- diag(comp_diag$covariances[, k], nrow = d)
  }
  comp_full <- comp_diag
  comp_full$covariances <- full_arr
  comp_full$covariance_type <- "full"

  pts <- rbind(c(-2, -1), c(2, 1), c(0, 0), c(1.2, 0.4))
  expect_equal(
    density_posterior(comp_diag, pts),
    density_posterior(comp_full, pts),
    tolerance = 1e-10
  )
})
