
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
  n_grid_y = grid_size,
  y_levels = seq(1, grid_size, by = floor(grid_size / n_samps)),
  
)





