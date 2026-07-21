
# sim_landscape() generates a full spatially-explicit community across an
# n_grid_x by n_grid_y grid. An environmental gradient runs along the vertical
# (y) axis, scaled to [-1, 1]. Each species has a logit-scale intercept and
# gradient slope drawn from normal distributions, so occurrence probability in
# cell (x, j) is plogis(beta[s,1] + beta[s,2] * env[j]).
#
# sample_landscape() draws samples from the full landscape in two modes:
#   - Gradient transect: y_levels as a vector, n_reps x-positions per level.
#   - Single-gradient-point replication: y_levels as a scalar, n_reps > 1.

sim_landscape <- function(
  n_species = 50,
  n_grid_x = 100,
  n_grid_y = 100,
  beta_sd = c(2, 1),
  seed = NULL
) {

  if (!is.null(seed)) set.seed(seed)

  env <- seq(-1, 1, length.out = n_grid_y)

  beta <- cbind(
    rnorm(n_species, -2, beta_sd[1]),
    rnorm(n_species, 0, beta_sd[2])
  )

  # Full landscape: W_full[s, x, y] = presence of species s in cell (x, y)
  W_full <- array(0L, dim = c(n_species, n_grid_x, n_grid_y))
  for (j in seq_len(n_grid_y)) {
    p <- plogis(beta[, 1] + beta[, 2] * env[j])
    W_full[, , j] <- rbinom(n_species * n_grid_x, 1L, p)
  }

  list(
    W_full = W_full,
    beta = beta,
    S_true = n_species,
    env = env
  )
}


sample_landscape <- function(
  landscape,
  y_levels = NULL,
  n_reps = 1
) {

  W_full <- landscape$W_full
  env <- landscape$env
  n_grid_x <- landscape$n_grid_x
  n_grid_y <- landscape$n_grid_y
  n_species <- dim(W_full)[1]

  if (is.null(y_levels)) y_levels <- seq_len(n_grid_y)
  if (n_reps > n_grid_x) {
    stop(paste("n_reps cannot exceed n_grid_x (", n_grid_x, ")"))
  }

  sampled_cells <- do.call(rbind, lapply(y_levels, function(j) {
    x_idx <- sample(seq_len(n_grid_x), n_reps, replace = FALSE)
    data.frame(x_idx = x_idx, y_idx = j)
  }))

  n_plots <- nrow(sampled_cells)

  W <- matrix(0L, nrow = n_species, ncol = n_plots)
  for (plot in seq_len(n_plots)) {
    W[, plot] <- W_full[, sampled_cells$x_idx[plot], sampled_cells$y_idx[plot]]
  }

  detected <- rowSums(W) > 0

  list(
    W = W,
    data = data.frame(env = env[sampled_cells$y_idx]),
    S_obs = sum(detected),
    sampled_cells = sampled_cells
  )
}
