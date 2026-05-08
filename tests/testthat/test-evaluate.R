test_that("evaluate_compression accepts compressed_fit wrapper (compress_brmsfit-style)", {
  draws <- make_two_blob_draws()
  comp <- compress_posterior(draws, method = "mclust", n_components = 2)
  wrapped <- structure(
    list(compressed = comp, structure = list(shell = TRUE)),
    class = c("compressed_brmsfit", "compressed_fit", "list")
  )
  fid <- evaluate_compression(
    wrapped,
    reference_draws = draws,
    seed            = 1L,
    n_self_reps     = 5L,
    max_n           = 600L
  )
  expect_s3_class(fid, "compression_fidelity")
  expect_equal(fid$n_params, 2L)
})

test_that("evaluate_compression returns a compression_fidelity object", {
  draws <- make_two_blob_draws()
  comp <- compress_posterior(draws, method = "mclust", n_components = 2)

  fid <- evaluate_compression(
    comp,
    reference_draws = draws,
    seed            = 1L,
    n_self_reps     = 5L,
    max_n           = 600L
  )

  expect_s3_class(fid, "compression_fidelity")
  expect_named(fid$metrics, c("energy", "c2st"))
  expect_true(is.finite(fid$reproduction_pct))
  expect_gte(fid$reproduction_pct, 0)
  expect_lte(fid$reproduction_pct, 100)
  expect_equal(fid$n_params, 2L)
  expect_equal(fid$method, "mclust")
})

test_that("a faithful compression scores high on both metrics", {
  draws <- make_two_blob_draws(n_per = 1500L)
  comp <- compress_posterior(draws, method = "mclust", n_components = 2)

  fid <- evaluate_compression(
    comp,
    reference_draws = draws,
    seed            = 1L,
    n_self_reps     = 5L,
    max_n           = 600L
  )

  expect_gt(fid$metrics$energy$reproduction_pct, 80)
  expect_gt(fid$metrics$c2st$reproduction_pct,   60)
  expect_lt(abs(fid$metrics$c2st$auc - 0.5),     0.2)
})

test_that("a deliberately bad compression scores lower than a good one", {
  draws <- make_two_blob_draws(n_per = 1500L)

  good <- compress_posterior(draws, method = "mclust", n_components = 2)
  bad  <- compress_posterior(draws, method = "mclust", n_components = 1L,
                             model_name = "EII")

  fid_good <- evaluate_compression(
    good, reference_draws = draws,
    metric = "energy", seed = 1L, n_self_reps = 5L, max_n = 600L
  )
  fid_bad  <- evaluate_compression(
    bad,  reference_draws = draws,
    metric = "energy", seed = 1L, n_self_reps = 5L, max_n = 600L
  )

  expect_gt(fid_good$reproduction_pct, fid_bad$reproduction_pct)
})

test_that("diagonal-only mclust scores lower than full covariance on a correlated posterior", {
  set.seed(2)
  Sigma <- matrix(c(1.0, 0.85, 0.85, 1.0), 2, 2)
  draws <- mvtnorm::rmvnorm(2000, mean = c(0, 0), sigma = Sigma)
  colnames(draws) <- c("alpha", "beta")

  comp_diag <- compress_posterior(
    draws, method = "mclust", n_components = 1L,
    model_name = c("EII", "VII", "EEI", "VEI", "EVI", "VVI")
  )
  comp_full <- compress_posterior(
    draws, method = "mclust", n_components = 1L,
    model_name = c("EEE", "VVV")
  )

  fid_diag <- evaluate_compression(
    comp_diag, draws,
    seed = 1L, n_self_reps = 10L, max_n = 1000L
  )
  fid_full <- evaluate_compression(
    comp_full, draws,
    seed = 1L, n_self_reps = 10L, max_n = 1000L
  )

  expect_lt(fid_diag$reproduction_pct, fid_full$reproduction_pct - 20)
  expect_lt(
    fid_diag$metrics$energy$reproduction_pct,
    fid_full$metrics$energy$reproduction_pct
  )
  expect_lt(
    fid_diag$metrics$c2st$reproduction_pct,
    fid_full$metrics$c2st$reproduction_pct
  )

  expect_gt(fid_full$reproduction_pct, 90)
  expect_lt(fid_diag$reproduction_pct, 70)
})

test_that("evaluate_compression accepts a single metric", {
  draws <- make_two_blob_draws()
  comp <- compress_posterior(draws, method = "mclust", n_components = 2)

  fid_e <- evaluate_compression(
    comp, draws, metric = "energy",
    seed = 1L, n_self_reps = 5L, max_n = 400L
  )
  expect_named(fid_e$metrics, "energy")

  fid_c <- evaluate_compression(
    comp, draws, metric = "c2st",
    seed = 1L, max_n = 400L
  )
  expect_named(fid_c$metrics, "c2st")
})

test_that("evaluate_compression errors on mismatched parameters", {
  draws <- make_two_blob_draws()
  comp <- compress_posterior(draws, method = "mclust", n_components = 2)

  bad_ref <- matrix(rnorm(20), ncol = 2,
                    dimnames = list(NULL, c("foo", "bar")))
  expect_error(evaluate_compression(comp, bad_ref), regexp = "matching|found")
})

test_that("evaluate_compression supports a held-out reference set", {
  draws <- make_two_blob_draws(n_per = 1500L)
  set.seed(42)
  idx <- sample.int(nrow(draws), size = floor(0.8 * nrow(draws)))
  comp <- compress_posterior(draws[idx, ], method = "mclust", n_components = 2)

  fid <- evaluate_compression(
    comp,
    reference_draws = draws[-idx, ],
    metric          = "energy",
    seed            = 1L,
    n_self_reps     = 5L,
    max_n           = 400L
  )
  expect_s3_class(fid, "compression_fidelity")
  expect_gt(fid$reproduction_pct, 50)
})

test_that("print.compression_fidelity prints a reproduction line", {
  draws <- make_two_blob_draws()
  comp <- compress_posterior(draws, method = "mclust", n_components = 2)

  fid <- evaluate_compression(
    comp, draws,
    metric = "energy", seed = 1L, n_self_reps = 5L, max_n = 400L
  )
  expect_output(print(fid), "reproduction")
})

test_that("default C2ST classifier is ranger", {
  draws <- make_two_blob_draws()
  comp <- compress_posterior(draws, method = "mclust", n_components = 2)

  fid <- evaluate_compression(
    comp, draws,
    metric = "c2st",
    seed = 1L, max_n = 400L
  )
  expect_equal(fid$metrics$c2st$classifier, "ranger")
  expect_true(is.finite(fid$metrics$c2st$auc))
})

test_that("C2ST ranger works with brms-like non-formula-safe colnames", {
  draws <- make_two_blob_draws()
  colnames(draws) <- c("b[1]", "cor__(1,2):(Intercept__foo)")

  comp <- compress_posterior(draws, method = "mclust", n_components = 2)
  fid <- evaluate_compression(
    comp, draws,
    metric = "c2st",
    seed = 1L, max_n = 400L
  )
  expect_equal(fid$metrics$c2st$classifier, "ranger")
  expect_true(is.finite(fid$metrics$c2st$auc))
})

test_that("classifier = knn uses k-NN path", {
  draws <- make_two_blob_draws()
  comp <- compress_posterior(draws, method = "mclust", n_components = 2)

  fid <- evaluate_compression(
    comp, draws,
    metric = "c2st", classifier = "knn",
    seed = 1L, max_n = 300L
  )
  expect_equal(fid$metrics$c2st$classifier, "knn")
  expect_true(is.finite(fid$metrics$c2st$auc))
})
