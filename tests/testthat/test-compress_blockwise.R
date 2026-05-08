make_blockwise_draws <- function(n = 800L, seed = 1L) {
  set.seed(seed)
  # Cluster 1: 3 strongly correlated parameters (linear).
  c1a <- rnorm(n)
  c1b <- c1a + 0.05 * rnorm(n)
  c1c <- c1a - 0.05 * rnorm(n)
  # Cluster 2: 3 strongly correlated parameters around a different mean.
  c2a <- rnorm(n, mean = 5)
  c2b <- c2a + 0.05 * rnorm(n)
  c2c <- c2a - 0.05 * rnorm(n)
  # Remainder: 4 (nearly) independent parameters.
  rem <- matrix(rnorm(n * 4L), ncol = 4L)
  draws <- cbind(c1a, c1b, c1c, c2a, c2b, c2c, rem)
  colnames(draws) <- c(
    "g1_a", "g1_b", "g1_c",
    "g2_a", "g2_b", "g2_c",
    "r1", "r2", "r3", "r4"
  )
  draws
}

test_that("compress_posterior(partition = ...) builds a blockwise object", {
  draws <- make_blockwise_draws()
  cm    <- posterior_correlation(draws)
  part  <- partition_parameters_clusters(cm, threshold = 0.5, min_size = 2L)

  comp <- compress_posterior(
    draws,
    method               = "mclust",
    n_components         = 2L,
    partition            = part,
    remainder_model_name = "VVI"
  )

  expect_s3_class(comp, "posterior_compressed_blockwise")
  expect_s3_class(comp, "posterior_compressed")
  expect_equal(comp$method, "blockwise")
  expect_equal(comp$base_method, "mclust")
  expect_equal(comp$n_params, ncol(draws))
  expect_equal(comp$n_draws, nrow(draws))
  expect_setequal(comp$param_names, colnames(draws))

  expect_true("remainder" %in% names(comp$blocks))
  expect_true(comp$has_remainder)
  expect_gt(length(comp$cluster_block_names), 0L)

  # Each block compressed object should be a posterior_compressed.
  for (b in comp$blocks) {
    expect_s3_class(b, "posterior_compressed")
    expect_true(all(b$param_names %in% colnames(draws)))
  }

  # Remainder block should be restricted to diagonal/spherical families.
  rem_block <- comp$blocks$remainder
  expect_true(rem_block$model_name %in% c("EII", "VII", "EEI", "EVI", "VEI", "VVI"))
  expect_equal(rem_block$covariance_type, "diagonal")
})

test_that("sample_posterior on a blockwise object recovers all params and order", {
  draws <- make_blockwise_draws()
  cm    <- posterior_correlation(draws)
  part  <- partition_parameters_clusters(cm, threshold = 0.5, min_size = 2L)

  comp <- compress_posterior(
    draws,
    method       = "mclust",
    n_components = 2L,
    partition    = part
  )

  s <- sample_posterior(comp, n_draws = 200L)
  expect_true(is.matrix(s))
  expect_equal(dim(s), c(200L, ncol(draws)))
  expect_equal(colnames(s), comp$param_names)
  expect_false(anyNA(s))

  # Marginal sanity: means within a few SDs of original.
  orig_means <- colMeans(draws)
  new_means  <- colMeans(s)
  expect_true(all(abs(orig_means - new_means) < 0.5))
})

test_that("density_posterior on blockwise sums per-block log-densities", {
  draws <- make_blockwise_draws(n = 400L)
  cm    <- posterior_correlation(draws)
  part  <- partition_parameters_clusters(cm, threshold = 0.5, min_size = 2L)

  comp <- compress_posterior(
    draws,
    method       = "mclust",
    n_components = 2L,
    partition    = part
  )

  test_pts <- draws[1:5, , drop = FALSE]
  ld <- density_posterior(comp, test_pts, log = TRUE)
  expect_equal(length(ld), 5L)
  expect_true(all(is.finite(ld)))

  # Compare to manual sum across blocks.
  manual <- numeric(5L)
  for (b in comp$blocks) {
    manual <- manual + density_posterior(
      b, test_pts[, b$param_names, drop = FALSE], log = TRUE
    )
  }
  expect_equal(ld, manual)
})

test_that("blockwise errors when partition references missing parameters", {
  draws <- make_blockwise_draws()
  cm    <- posterior_correlation(draws)
  part  <- partition_parameters_clusters(cm, threshold = 0.5, min_size = 2L)

  # Inject a parameter name that is not in draws.
  part$blocks$cluster_1 <- c(part$blocks$cluster_1, "ghost_param")

  expect_error(
    compress_posterior(
      draws,
      method       = "mclust",
      n_components = 2L,
      partition    = part
    ),
    "ghost_param"
  )
})

test_that("summary / get_compressed_block work on blockwise objects", {
  draws <- make_blockwise_draws()
  cm    <- posterior_correlation(draws)
  part  <- partition_parameters_clusters(cm, threshold = 0.5, min_size = 2L)

  comp <- compress_posterior(
    draws,
    method       = "mclust",
    n_components = 2L,
    partition    = part
  )

  s <- summary(comp)
  expect_s3_class(s, "summary.posterior_compressed_blockwise")
  expect_equal(s$method, "blockwise")
  expect_equal(s$base_method, "mclust")
  expect_true(is.data.frame(s$blocks))
  expect_setequal(
    s$blocks$block,
    names(comp$blocks)
  )
  expect_true(all(c("cluster", "remainder") %in% s$blocks$role))

  rem <- get_compressed_block(comp, "remainder")
  expect_s3_class(rem, "posterior_compressed")
  expect_equal(rem$covariance_type, "diagonal")

  expect_error(get_compressed_block(comp, "ghost_block"),
               "ghost_block")
})

test_that("density_posterior on blockwise validates inputs by name", {
  draws <- make_blockwise_draws(n = 400L)
  cm    <- posterior_correlation(draws)
  part  <- partition_parameters_clusters(cm, threshold = 0.5, min_size = 2L)

  comp <- compress_posterior(
    draws, method = "mclust", n_components = 2L, partition = part
  )

  bad <- draws[1:3, c("g1_a", "g1_b"), drop = FALSE]
  expect_error(
    density_posterior(comp, bad),
    "missing"
  )
})

test_that("blockwise integrates with reconstruct_brmsfit() (full round-trip)", {
  expect_true(requireNamespace("brms", quietly = TRUE))
  expect_true(requireNamespace("cmdstanr", quietly = TRUE))
  expect_true(requireNamespace("posterior", quietly = TRUE))
  expect_no_error(cmdstanr::cmdstan_path())

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

  draws_mat <- as.matrix(posterior::as_draws_matrix(fit))
  cm        <- posterior_correlation(draws_mat)
  # min_size = 2 so the small example fit produces clusters at all.
  part      <- partition_parameters_clusters(cm, threshold = 0.3, min_size = 2L)

  res <- compress_brmsfit(
    fit,
    method               = "mclust",
    n_components         = 2L,
    partition            = part,
    remainder_model_name = "VVI"
  )

  expect_s3_class(res$compressed, "posterior_compressed_blockwise")
  expect_s3_class(res$structure,  "brmsfit")
  expect_s3_class(res$structure,  "brmsfit_stripped")

  recon <- reconstruct_brmsfit(res, n_draws = 400)
  expect_s3_class(recon, "brmsfit")
  expect_true(methods::is(recon$fit, "stanfit"))
  draws_mat2 <- posterior::as_draws_matrix(recon)
  expect_true(nrow(draws_mat2) > 0L)
  expect_equal(attr(recon, "compression_method"), "blockwise")
  expect_equal(attr(recon, "compression_base_method"), "mclust")
  expect_true(is.integer(attr(recon, "compression_blocks")))
  expect_equal(
    sum(attr(recon, "compression_blocks")),
    res$compressed$n_params
  )

  # Both per-chain AND stanfit-level `sampler_params` attributes must
  # be restored so auto-printing the reconstructed brmsfit doesn't trip
  # rstan's diagnostic summary (regression for "do.call(cbind, attr(x,
  # 'sampler_params')) : second argument must be a list").
  for (sm in recon$fit@sim$samples) {
    sp <- attr(sm, "sampler_params")
    expect_true(is.list(sp))
    sp_mat <- do.call(cbind, sp)
    expect_true(is.matrix(sp_mat))
    expect_true("divergent__" %in% colnames(sp_mat))
  }
  sf_sp <- attr(recon$fit, "sampler_params")
  expect_true(is.list(sf_sp))
  expect_equal(length(sf_sp), length(recon$fit@sim$samples))
  for (mat in sf_sp) {
    expect_true(is.matrix(mat))
    expect_true("divergent__" %in% colnames(mat))
  }
  expect_no_error(capture.output(print(recon)))
})

test_that("plain list of character vectors coerces to blockwise partition", {
  draws <- make_blockwise_draws()
  part <- list(
    c("g1_a", "g1_b", "g1_c"),
    c("g2_a", "g2_b", "g2_c")
  )
  comp <- compress_posterior(
    draws,
    method               = "mclust",
    n_components         = 2L,
    partition            = part,
    remainder_model_name = "VVI"
  )
  expect_s3_class(comp, "posterior_compressed_blockwise")
  rem <- comp$blocks$remainder
  expect_setequal(rem$param_names, c("r1", "r2", "r3", "r4"))
})

test_that("partition_blocks() resolves tidyselect against draws colnames", {
  draws <- make_blockwise_draws()
  pb <- partition_blocks(
    c("g1_a", "g1_b", "g1_c"),
    tidyselect::starts_with("g2_"),
    c("r1", "r2")
  )
  comp <- compress_posterior(
    draws,
    method               = "mclust",
    n_components         = 2L,
    partition            = pb,
    remainder_model_name = "VVI"
  )
  expect_s3_class(comp, "posterior_compressed_blockwise")
  b1 <- comp$blocks[[1]]
  expect_setequal(b1$param_names, c("g1_a", "g1_b", "g1_c"))
  b2 <- comp$blocks[[2]]
  expect_setequal(b2$param_names, c("g2_a", "g2_b", "g2_c"))
  rem <- comp$blocks$remainder
  expect_true(all(c("r1", "r2") %in% rem$param_names))
  expect_true("r3" %in% rem$param_names && "r4" %in% rem$param_names)
})

test_that("plain list may include a formula with tidyselect on RHS", {
  draws <- make_blockwise_draws()
  part <- list(
    c("g1_a", "g1_b", "g1_c"),
    ~ tidyselect::starts_with("g2_")
  )
  comp <- compress_posterior(
    draws,
    method       = "mclust",
    n_components = 2L,
    partition    = part
  )
  expect_s3_class(comp, "posterior_compressed_blockwise")
  expect_equal(length(comp$cluster_block_names), 2L)
})

test_that("overlap between manual blocks respects declaration order", {
  draws <- make_blockwise_draws()
  pb <- partition_blocks(
    c("g1_a", "g2_a"),
    c("g2_a", "g2_b", "g2_c")
  )
  comp <- compress_posterior(
    draws,
    method       = "mclust",
    n_components = 2L,
    partition    = pb
  )
  b2 <- comp$blocks[[2]]
  expect_false("g2_a" %in% b2$param_names)
  expect_true(all(c("g2_b", "g2_c") %in% b2$param_names))
})

test_that("partition names absent from draws are ignored with warning", {
  draws <- make_blockwise_draws()
  part <- list(
    c("g1_a", "g1_b", "g1_c", "ghost_a", "ghost_b"),
    c("g2_a", "g2_b", "g2_c")
  )
  expect_warning(
    comp <- compress_posterior(
      draws,
      method       = "mclust",
      n_components = 2L,
      partition    = part
    ),
    "not found in draws"
  )
  expect_s3_class(comp, "posterior_compressed_blockwise")
  expect_false(any(c("ghost_a", "ghost_b") %in% comp$param_names))
})

test_that("draw columns missing from partition are added to remainder with warning", {
  draws <- make_blockwise_draws()
  part <- structure(
    list(
      blocks    = list(
        cluster_1 = c("g1_a", "g1_b", "g1_c"),
        cluster_2 = c("g2_a", "g2_b", "g2_c")
      ),
      remainder = character(0)
    ),
    class = c("poco_partition_clusters", "list")
  )
  expect_warning(
    comp <- compress_posterior(
      draws,
      method       = "mclust",
      n_components = 2L,
      partition    = part
    ),
    "not listed in any cluster or remainder"
  )
  rem <- comp$blocks$remainder
  expect_true(all(c("r1", "r2", "r3", "r4") %in% rem$param_names))
})

test_that("poco_partition_clusters: names only in partition are dropped with warning", {
  draws <- make_blockwise_draws()
  part <- structure(
    list(
      blocks = list(
        cluster_1 = c("g1_a", "g1_b", "g1_c", "phantom_x"),
        cluster_2 = c("g2_a", "g2_b", "g2_c")
      ),
      remainder = c("r1", "r2", "r3", "r4", "phantom_y")
    ),
    class = c("poco_partition_clusters", "list")
  )
  expect_warning(
    comp <- compress_posterior(
      draws,
      method       = "mclust",
      n_components = 2L,
      partition    = part
    ),
    "not present in draws"
  )
  expect_false(any(c("phantom_x", "phantom_y") %in% unlist(
    lapply(comp$blocks, function(b) b$param_names),
    use.names = FALSE
  )))
})

test_that("blockwise compress_posterior runs with explicit cluster_BPPARAM", {
  skip_if_not_installed("BiocParallel")
  draws <- make_blockwise_draws()
  cm    <- posterior_correlation(draws)
  part  <- partition_parameters_clusters(cm, threshold = 0.5, min_size = 2L)
  comp <- compress_posterior(
    draws,
    method          = "mclust",
    n_components    = 2L,
    partition       = part,
    cluster_BPPARAM = BiocParallel::SerialParam()
  )
  expect_s3_class(comp, "posterior_compressed_blockwise")
  expect_equal(comp$n_params, ncol(draws))
})

test_that("cluster_BPPARAM = FALSE forces sequential cluster compression", {
  draws <- make_blockwise_draws()
  cm    <- posterior_correlation(draws)
  part  <- partition_parameters_clusters(cm, threshold = 0.5, min_size = 2L)
  comp <- compress_posterior(
    draws,
    method          = "mclust",
    n_components    = 2L,
    partition       = part,
    cluster_BPPARAM = FALSE
  )
  expect_s3_class(comp, "posterior_compressed_blockwise")
  expect_equal(comp$n_params, ncol(draws))
})
