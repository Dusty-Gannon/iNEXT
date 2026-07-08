simulate_incidence_data <- function(X, beta, seed = 42) {
  set.seed(seed)
  P <- plogis(beta %*% X)
  matrix(
    rbinom(length(P), size = 1L, prob = as.numeric(P)),
    nrow = nrow(beta), ncol = ncol(P)
  )
}

make_incidence_data <- function(N = 20, S = 5, x_cov = FALSE, seed = 1L) {
  set.seed(seed)
  X <- matrix(1, nrow = 1, ncol = N)
  if(x_cov){
    X <- rbind(X, rnorm(N))
  }
  beta  <- matrix(rnorm(S * nrow(X), sd = 0.5), S, nrow(X))
  W <- simulate_incidence_data(X, beta, seed = seed)
  list(W = W, X = X, beta = beta)
}
