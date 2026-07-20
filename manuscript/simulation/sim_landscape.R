
# Simulate a spatially-explicit community with an environmental gradient and
# sample it using subplots.
#
# The landscape is divided into an n_grid_x by n_grid_y grid. An environmental
# gradient runs along the vertical (y) axis, with one covariate value per
# y-level. Each species has a logit-scale intercept and gradient slope drawn
# from normal distributions, so occurrence probability in each cell follows
# P(s present | y-level j) = plogis(beta[s, 1] + beta[s, 2] * env[j]).
#
# Sampling is specified by y_levels (which rows of the grid to sample) and
# n_reps (how many distinct x-positions to draw within each y-level). This
# lets you replicate samples at the same point along the gradient, matching
# the design passed to iNEXT_stan via its formula and data arguments.

sim_landscape <- function(
  n_species = 50,
  n_grid_x = 10,
  n_grid_y = 10,
  y_levels = NULL,
  n_reps = 1,
  beta_sd = c(1, 1),
  seed = NULL
) {

  if (!is.null(seed)) set.seed(seed)
  if (n_reps > n_grid_x)
    stop(paste("n_reps cannot exceed n_grid_x (", n_grid_x, ")"))

  # ---- Environmental gradient ----
  # Scaled to [-1, 1] across the y-grid rows
  env <- seq(-1, 1, length.out = n_grid_y)

  # ---- Species-level coefficients: intercept and gradient slope ----
  beta <- cbind(
    rnorm(n_species, 0, beta_sd[1]),
    rnorm(n_species, 0, beta_sd[2])
  )

  # ---- Sampling design ----
  if (is.null(y_levels)) y_levels <- 1:n_grid_y

  sampled_cells <- do.call(rbind, lapply(y_levels, function(j) {
    x_idx <- sample(1:n_grid_x, n_reps, replace = FALSE)
    data.frame(x_idx = x_idx, y_idx = j)
  }))

  n_plots <- nrow(sampled_cells)

  # ---- Simulate presence/absence in each sampled cell ----
  W_full <- matrix(0L, nrow = n_species, ncol = n_plots)
  for (plot in 1:n_plots) {
    j <- sampled_cells$y_idx[plot]
    p <- plogis(beta[, 1] + beta[, 2] * env[j])
    W_full[, plot] <- rbinom(n_species, 1L, p)
  }

  # ---- Drop species never detected ----
  detected <- rowSums(W_full) > 0
  W <- W_full[detected, , drop = FALSE]

  list(
    W = W,
    data = data.frame(env = env[sampled_cells$y_idx]),
    beta = beta,
    S_true = n_species,
    S_obs = sum(detected),
    sampled_cells = sampled_cells,
    env = env
  )
}
