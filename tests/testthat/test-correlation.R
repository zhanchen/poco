make_block_draws <- function(n = 800L, seed = 1L) {
  set.seed(seed)
  # Block A: 3 strongly correlated parameters.
  a1 <- rnorm(n)
  a2 <- a1 + 0.05 * rnorm(n)
  a3 <- a1 - 0.05 * rnorm(n)
  # Block B: 4 (nearly) independent parameters.
  b <- matrix(rnorm(n * 4L), ncol = 4L)
  draws <- cbind(a1, a2, a3, b)
  colnames(draws) <- c("a1", "a2", "a3", "b1", "b2", "b3", "b4")
  draws
}

test_that("posterior_correlation returns a labelled symmetric matrix", {
  draws <- make_block_draws()
  cm <- posterior_correlation(draws)

  expect_true(is.matrix(cm))
  expect_equal(nrow(cm), ncol(cm))
  expect_equal(rownames(cm), colnames(cm))
  expect_equal(rownames(cm), colnames(draws))
  expect_equal(diag(cm), setNames(rep(1, ncol(cm)), colnames(cm)))
  expect_equal(cm, t(cm), tolerance = 1e-12)
  expect_equal(attr(cm, "type"), "correlation")
  expect_equal(attr(cm, "method"), "pearson")
  expect_equal(attr(cm, "n_draws"), nrow(draws))
  expect_equal(attr(cm, "n_params"), ncol(draws))
})

test_that("posterior_correlation supports covariance and method choice", {
  draws <- make_block_draws()
  cv <- posterior_correlation(draws, type = "covariance")
  expect_equal(attr(cv, "type"), "covariance")
  expect_true(is.na(attr(cv, "method")))

  cs <- posterior_correlation(draws, method = "spearman")
  expect_equal(attr(cs, "method"), "spearman")
})

test_that("posterior_correlation drops zero-variance parameters", {
  draws <- make_block_draws()
  draws <- cbind(draws, dead = 1)
  cm <- posterior_correlation(draws)
  expect_false("dead" %in% rownames(cm))
  expect_equal(attr(cm, "dropped_zero_var"), "dead")
})

test_that("partition_parameters threshold rule isolates block A", {
  draws <- make_block_draws()
  cm <- posterior_correlation(draws)

  part <- partition_parameters(cm, rule = "threshold", threshold = 0.5)
  expect_s3_class(part, "poco_partition")
  expect_setequal(part$block_a, c("a1", "a2", "a3"))
  expect_setequal(part$block_b, c("b1", "b2", "b3", "b4"))
  expect_equal(length(part$abs_max_offdiag), ncol(cm))
})

test_that("partition_parameters top_n / top_prop / min_degree work", {
  draws <- make_block_draws()
  cm <- posterior_correlation(draws)

  p_topn <- partition_parameters(cm, rule = "top_n", n = 3L)
  expect_setequal(p_topn$block_a, c("a1", "a2", "a3"))

  p_prop <- partition_parameters(cm, rule = "top_prop", prop = 3 / 7)
  expect_equal(length(p_prop$block_a), 3L)

  p_deg <- partition_parameters(
    cm, rule = "min_degree", threshold = 0.5, min_degree = 2L
  )
  expect_setequal(p_deg$block_a, c("a1", "a2", "a3"))
})

test_that("partition_parameters validates required args per rule", {
  draws <- make_block_draws()
  cm <- posterior_correlation(draws)

  expect_error(partition_parameters(cm, rule = "threshold"), "threshold")
  expect_error(partition_parameters(cm, rule = "top_n"), "n")
  expect_error(partition_parameters(cm, rule = "top_prop"), "prop")
  expect_error(
    partition_parameters(cm, rule = "min_degree", threshold = 0.5),
    "min_degree"
  )
})

test_that("plot_posterior_correlation returns ggplot objects when ggplot2 available", {
  skip_if_not_installed("ggplot2")
  draws <- make_block_draws()
  cm <- posterior_correlation(draws)

  ps <- plot_posterior_correlation(cm, threshold = 0.3)
  expect_named(ps, c("heatmap", "hist", "max_abs"))
  for (p in ps) expect_s3_class(p, "ggplot")

  one <- plot_posterior_correlation(cm, which = "max_abs", threshold = 0.3)
  expect_s3_class(one, "ggplot")
})

test_that("partition_parameters_clusters returns remainder and rest alias", {
  draws <- make_block_draws()
  cm <- posterior_correlation(draws)
  pc <- partition_parameters_clusters(cm, threshold = 0.5, min_size = 2L)

  expect_s3_class(pc, "poco_partition_clusters")
  expect_true("remainder" %in% names(pc))
  expect_true("rest" %in% names(pc))
  expect_identical(pc$rest, pc$remainder)
  expect_equal(
    length(pc$remainder) + sum(lengths(pc$blocks)),
    length(pc$cluster_id)
  )
  expect_equal(sum(is.na(pc$cluster_id)), length(pc$remainder))
})

test_that("partition_parameters_clusters interprets min_size in (0,1) as proportion", {
  draws <- make_block_draws()
  cm <- posterior_correlation(draws)
  n <- ncol(cm)
  prop <- 0.25
  pc <- suppressMessages(
    partition_parameters_clusters(cm, threshold = 0.5, min_size = prop)
  )
  expect_equal(pc$params$min_size, as.integer(ceiling(prop * n)))
  expect_equal(pc$params$min_size_proportion, prop)
})
