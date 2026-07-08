
# ---- Layer 1: simulator unit tests (no Stan) --------------------------------

test_that("simulate_incidence_data returns a binary matrix with correct dims", {
  td <- make_detection_fixture()
  expect_true(all(td$W %in% c(0L, 1L)))
  expect_equal(dim(td$W), c(5L, 20L))
})

test_that("simulate_incidence_data is deterministic given the same seed", {
  td1 <- make_detection_fixture(seed = 7L)
  td2 <- make_detection_fixture(seed = 7L)
  expect_identical(td1$W, td2$W)
})

test_that("empirical detection rates match logit-scale beta (large-N check)", {
  X    <- rbind(rep(1, 500))
  beta <- matrix(c(1.0, -1.0), 2, 1)
  W    <- simulate_incidence_data(X, beta, seed = 1)
  expect_equal(rowMeans(W)[1], plogis(1),  tolerance = 0.05)
  expect_equal(rowMeans(W)[2], plogis(-1), tolerance = 0.05)
})


# ---- Layer 2: input validation tests (no Stan) ------------------------------

test_that("iNEXT_stan rejects non-binary W", {
  td <- make_detection_fixture()
  td$W[1, 1] <- 2L
  expect_error(
    iNEXT_stan(td$W, ~x_cov, td$data, chains = 1, iter = 10),
    regexp = "0 or 1"
  )
})

test_that("iNEXT_stan rejects W / data dimension mismatch", {
  td      <- make_detection_fixture(N = 20)
  td$data <- td$data[1:15, , drop = FALSE]
  expect_error(
    iNEXT_stan(td$W, ~x_cov, td$data, chains = 1, iter = 10),
    regexp = "Dimension mismatch"
  )
})


# ---- Layer 3: pipeline smoke test (needs rstan) -----------------------------

test_that("full pipeline W -> Stan -> iNEXT.Sam -> ggiNEXT runs without error", {
  skip_if_not_installed("rstan")
  skip_on_cran()

  td  <- make_detection_fixture(N = 30, S = 6, seed = 42)
  res <- iNEXT_stan(
    W       = td$W,
    formula = ~x_cov,
    data    = td$data,
    chains  = 1, iter = 300, warmup = 150, seed = 1
  )

  expect_type(res, "list")
  expect_named(res, c("incfreq", "stanfit"), ignore.order = TRUE)
  expect_length(res$incfreq, 1L)   # one design point when datnew = NULL

  freq <- res$incfreq[[1]]
  expect_true(freq[1] > 0)
  expect_true(all(freq >= 0))
  expect_equal(freq, round(freq))

  sam_out <- iNEXT.Sam(freq, q = 0, se = FALSE)
  expect_s3_class(sam_out, "data.frame")
  expect_true(nrow(sam_out) > 0)

  inext_out <- iNEXT(freq, q = c(0, 1, 2), datatype = "incidence_freq")
  expect_s3_class(inext_out, "iNEXT")

  p <- ggiNEXT(inext_out)
  expect_s3_class(p, "ggplot")
})
