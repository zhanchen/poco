test_that("compress_fit rejects non-fit / non-character inputs", {
  expect_error(compress_fit(NULL), "CmdStanMCMC")
  expect_error(compress_fit(123), "CmdStanMCMC")
})

test_that("compress_fit reads cmdstan CSVs (smoke test)", {
  skip_if_not_installed("cmdstanr")
  skip_on_cran()

  tryCatch(cmdstanr::cmdstan_path(), error = function(e) {
    skip("cmdstan not installed")
  })

  stan_code <- "
data { int<lower=1> N; vector[N] y; }
parameters { real mu; real<lower=0> sigma; }
model { y ~ normal(mu, sigma); }
"
  stan_file <- cmdstanr::write_stan_file(stan_code)
  mod <- cmdstanr::cmdstan_model(stan_file)
  set.seed(1)
  fit <- mod$sample(
    data = list(N = 50L, y = rnorm(50)),
    chains = 2, parallel_chains = 1,
    iter_warmup = 200, iter_sampling = 200,
    refresh = 0
  )

  comp <- compress_fit(fit, method = "mclust", n_components = 2)
  expect_s3_class(comp, "posterior_compressed_mclust")
  expect_true("mu" %in% comp$param_names)

  s <- sample_posterior(comp, n_draws = 100)
  expect_equal(nrow(s), 100L)

  comp_csv <- compress_fit(
    fit$output_files(),
    method = "mclust",
    n_components = 2
  )
  expect_s3_class(comp_csv, "posterior_compressed_mclust")
})
