test_that("compress_sccomp rejects objects without an attached fit", {
  obj <- structure(list(), class = "sccomp")
  expect_error(compress_sccomp(obj), "attr")

  attr(obj, "fit") <- "not_a_cmdstan_fit"
  expect_error(compress_sccomp(obj), "cmdstanr backend")
})

test_that("reconstruct_sccomp rejects malformed input", {
  expect_error(reconstruct_sccomp("not a list"), "compressed.*structure")
})
