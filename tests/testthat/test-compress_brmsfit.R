test_that("print compressed_brmsfit dispatches safe printer (regression)", {
  comp <- structure(
    list(
      method      = "mclust",
      n_params    = 1L,
      n_components = 1L,
      n_draws     = 1L,
      param_names = "x"
    ),
    class = c("posterior_compressed_mclust", "posterior_compressed", "list")
  )
  st <- structure(
    list(fit = NULL),
    class = c("brmsfit_stripped", "brmsfit", "list")
  )
  res <- structure(
    list(compressed = comp, structure = st),
    class = c("compressed_brmsfit", "compressed_fit", "list")
  )
  expect_no_error(out <- capture.output(print(res)))
  expect_true(any(grepl("compressed", out, fixed = TRUE)))
})

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
  # Same gate as vignettes: rworkflows has no CmdStan; workflow_with_cmdstanr
  # installs it and runs this test in full.
  skip_if_not(
    instantiate::stan_cmdstan_exists(),
    "CmdStan not installed (full Stan coverage runs in workflow_with_cmdstanr)"
  )
  expect_true(requireNamespace("brms", quietly = TRUE))
  expect_true(requireNamespace("cmdstanr", quietly = TRUE))
  expect_true(requireNamespace("posterior", quietly = TRUE))

  set.seed(1)
  dat <- data.frame(x = rnorm(40), y = NA_real_)
  dat$y <- 1 + 0.5 * dat$x + rnorm(40)

  fit <- brms::brm(
    y ~ x,
    data = dat,
    chains = 2,
    iter = 600,
    warmup = 300,
    backend = "cmdstanr",
    refresh = 0,
    silent = 2
  )
  expect_no_error(ensure_brms_cmdstanr(fit))

  res <- compress_brmsfit(fit, method = "mclust", n_components = 2)
  expect_named(res, c("compressed", "structure"))
  expect_s3_class(res$compressed, "posterior_compressed_mclust")
  expect_s3_class(res$structure, "brmsfit")
  # The structure must also carry the stripped-shell tag so that
  # auto-printing it doesn't dispatch to brms::print.brmsfit and crash on
  # the empty stanfit (regression for "do.call(cbind, attr(x,
  # 'sampler_params')) : second argument must be a list").
  expect_s3_class(res$structure, "brmsfit_stripped")
  expect_no_error(capture.output(print(res$structure)))
  expect_no_error(capture.output(print(res)))

  recon <- reconstruct_brmsfit(res, n_draws = 400)
  expect_s3_class(recon, "brmsfit")
  expect_false(inherits(recon, "brmsfit_stripped"))
  expect_true(methods::is(recon$fit, "stanfit"))
  expect_true(length(recon$fit@sim$samples) > 0L)
  expect_true(length(recon$fit@sim$samples[[1]][["b_Intercept"]]) > 0L)
  draws_mat <- posterior::as_draws_matrix(recon)
  expect_true(nrow(draws_mat) > 0L)
  expect_equal(attr(recon, "compression_method"), "mclust")
})
