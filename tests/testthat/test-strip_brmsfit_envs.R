test_that("strip_brmsfit_envs validates input", {
  expect_error(strip_brmsfit_envs(list()), "single brmsfit")
  expect_error(strip_brmsfit_envs("x"), "single brmsfit")
})

test_that("strip_brmsfit_envs clears formula environments and keeps class", {
  env <- new.env(parent = emptyenv())
  env$blob <- rep(1, 1e5)
  f <- stats::as.formula(y ~ x)
  environment(f) <- env

  mock <- structure(
    list(
      formula = list(formula = f),
      data = data.frame(y = 1:3, x = 1:3),
      fit = NULL
    ),
    class = "brmsfit"
  )

  before <- length(serialize(mock, NULL, xdr = FALSE))
  out <- strip_brmsfit_envs(mock)
  after <- length(serialize(out, NULL, xdr = FALSE))

  expect_s3_class(out, "brmsfit")
  expect_true(is.null(environment(out$formula$formula)))
  expect_lt(after, before / 10)
})

test_that("lapply on brmsfit slots is not the supported API", {
  mock <- structure(list(formula = list()), class = "brmsfit")
  expect_error(
    lapply(mock, strip_brmsfit_envs),
    "single brmsfit"
  )
})
