
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

dat_incfreq <- rowSums(W_standard$W)

dat_incfreq <- c(n_samps, dat_incfreq)

std_estims <- iNEXT(dat_incfreq, datatype = "incidence_freq")






