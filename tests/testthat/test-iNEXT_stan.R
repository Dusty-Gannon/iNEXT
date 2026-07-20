
# ---- Layer 1: simulator unit tests (no Stan) --------------------------------

test_that("simulate_incidence_data returns a binary matrix with correct dims", {
  td <- make_incidence_data()
  expect_true(all(td$W %in% c(0L, 1L)))
  expect_equal(dim(td$W), c(5L, 20L))
})

test_that("simulate_incidence_data is deterministic given the same seed", {
  td1 <- make_incidence_data(seed = 7L)
  td2 <- make_incidence_data(seed = 7L)
  expect_identical(td1$W, td2$W)
})

test_that("empirical detection rates match logit-scale beta (large-N check)", {
  X <- rbind(rep(1, 500))
  beta <- matrix(c(1.0, -1.0), 2, 1)
  W <- simulate_incidence_data(X, beta, seed = 1)
  expect_equal(rowMeans(W)[1], plogis(1),  tolerance = 0.05)
  expect_equal(rowMeans(W)[2], plogis(-1), tolerance = 0.05)
})


# ---- Layer 2: input validation tests (no Stan) ------------------------------

test_that("iNEXT_stan rejects non-binary W", {
  td <- make_incidence_data()
  td$W[1, 1] <- 2L
  expect_error(
    iNEXT_stan(td$W, ~x_cov, td$data, chains = 1, iter = 10),
    regexp = "0 or 1"
  )
})

test_that("iNEXT_stan rejects W / data dimension mismatch", {
  td <- make_incidence_data(N = 20)
  td$data <- as.data.frame(td$X[, 1:15])
  
  expect_error(
    iNEXT_stan(td$W, ~ 1, td$data, chains = 1, iter = 10),
    regexp = "Dimension mismatch"
  )
})


# ---- Layer 3: pipeline smoke test (needs rstan) -----------------------------

test_that("full pipeline W -> Stan -> per-draw curves runs without error", {
  skip_if_not_installed("rstan")
  skip_on_cran()

  td  <- make_incidence_data(N = 30, S = 6, seed = 42)
  td$data <- as.data.frame(t(td$X))
  res <- iNEXT_stan(
    W = td$W,
    formula = ~ 1,
    data = td$data,
    chains = 1, iter = 300, warmup = 150, seed = 1,
    use_ndraws = 50L
  )

  expect_type(res, "list")
  expect_named(res, c("stanfit", "curves"), ignore.order = TRUE)
  expect_length(res$curves, 1L)   # one design point when datnew = NULL

  crv <- res$curves[[1]]
  expect_s3_class(crv, "data.frame")
  expect_named(crv, c("draw", "t", "order_q", "qD"))
  expect_true(all(crv$t > 0))
  expect_setequal(unique(crv$draw), sort(unique(crv$draw)))
})
