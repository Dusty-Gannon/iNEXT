
# ---- Libraries and source files ----

source(here::here("manuscript/simulations/sim_landscape.R"))
devtools::load_all(here::here())

# ---- Simulation ----

# Simulating a large landscape with an environmental gradient along the vertical axis
# The landscape is gridded, then occurrence patterns are determined based on 
# species-specific logistic model of occurrence in the cell. 

# setting the stage
grid_size <- 200
n_samps <- 30

data_full <- sim_landscape(
  n_grid_x = grid_size,
  n_grid_y = grid_size
)

## ---- Sampling from the landscape ----

y_idx <- sample(floor(grid_size / 3):(floor(grid_size / 3) * 2), size = 1)
# sample along the gradient
W_standard <- sample_landscape(
  data_full,
  y_levels = y_idx,
  n_reps = n_samps
)

W_m0 <- W_standard$W[!(rowSums(W_standard$W) == 0), ]

dat_incfreq <- rowSums(W_m0)

dat_incfreq <- c(n_samps, dat_incfreq)

std_estims <- iNEXT(dat_incfreq, datatype = "incidence_freq", q = c(0, 1, 2))

reg_check <- iNEXT_stan(W_m0, formula = ~ 1, data = W_standard$data, use_ndraws = 200, cores = 4)


# ---- Posterior Predictive Checks ----
# Goal: verify that the posterior predictive incidence frequencies are consistent
# with the observed data before comparing extrapolation trajectories.
# The key diagnostics for extrapolation are:
#   1. t_obs (effective sample size) vs actual N -- shifts the baseline t
#   2. S_obs per draw vs observed -- model calibration for richness
#   3. Qk frequency distribution -- Q1/Q2 drive Chao2 extrapolation; misfit here
#      explains divergent extrapolation curves

library(ggplot2)
library(patchwork)

# Observed quantities
obs_counts <- rowSums(W_m0)
obs_qk <- vapply(1:10, function(k) sum(obs_counts == k), integer(1L))
S_obs_data <- nrow(W_m0)

# Raw posterior predictive draws: [n_draws, S_obs + 1, N_new]
# Row 1 of dim-2 = t_obs_int; rows 2+ = per-species detection counts
incfreq_pp <- rstan::extract(reg_check$stanfit, pars = "incfreq_new")$incfreq_new
pp_t <- incfreq_pp[, 1L, 1L]
pp_sp <- incfreq_pp[, -1L, 1L]  # [n_draws x S_obs]

pp_S_obs <- rowSums(pp_sp > 0L)

pp_qk_mat <- t(apply(pp_sp, 1L, function(x)
  vapply(1:10, function(k) sum(x == k), integer(1L))
))

pp_qk_summary <- data.frame(
  k = 1:10,
  med = apply(pp_qk_mat, 2L, median),
  lo  = apply(pp_qk_mat, 2L, quantile, 0.025),
  hi  = apply(pp_qk_mat, 2L, quantile, 0.975),
  obs = obs_qk
)

p_tobs <- ggplot(data.frame(t = pp_t), aes(x = t)) +
  geom_histogram(bins = 40, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = n_samps, color = "tomato", linewidth = 1.2) +
  labs(
    title = "PPC: effective sample size vs observed N",
    x = "t_obs", y = "count",
    caption = paste0("red line = observed N = ", n_samps)
  )

p_sobs <- ggplot(data.frame(S = pp_S_obs), aes(x = S)) +
  geom_histogram(bins = 20, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = S_obs_data, color = "tomato", linewidth = 1.2) +
  labs(
    title = "PPC: species detected per draw vs observed",
    x = "S_obs per draw", y = "count",
    caption = paste0("red line = observed S_obs = ", S_obs_data)
  )

p_qk <- ggplot(pp_qk_summary, aes(x = k)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "steelblue", alpha = 0.3) +
  geom_line(aes(y = med), color = "steelblue", linewidth = 1) +
  geom_point(aes(y = obs), color = "tomato", size = 3) +
  geom_line(aes(y = obs), color = "tomato", linewidth = 1) +
  scale_x_continuous(breaks = 1:10) +
  labs(
    title = "PPC: incidence frequency distribution (Qk)",
    subtitle = "blue = PP median + 95% CI, red = observed",
    x = "k (number of sites species detected in)", y = "number of species"
  )

(p_tobs + p_sobs) / p_qk



