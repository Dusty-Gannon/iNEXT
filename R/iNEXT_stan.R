

.get_stan_model <- function() {
  if (!requireNamespace("rstan", quietly = TRUE))
    stop("Package 'rstan' is required. Install it with install.packages('rstan').")
  stan_file <- system.file("stan", "detection_glm.stan", package = "iNEXT")
  rstan::stan_model(stan_file, auto_write = TRUE)
}


iNEXT_stan <- function(W, formula, data, q = c(0, 1, 2), datnew = NULL, ...){
  
  # ---- Create objects and check them ----
  if(any(!(W == 0L | W == 1L))){
    stop("x must be a raw incidence matrix with S_obs rows and N columns.\n
         Found some elements not equal to 0 or 1")
  }
  
  ## ---- Model matrix construction ----
  X <- stats::model.matrix(formula, data) |> t()
  if(ncol(X) != ncol(W)) {
    stop(paste(
      "Dimension mismatch. Found",
      ncol(X),
      "sites/replicates in data, but",
      ncol(W), "in W."
    ))
  }
  
  mod <- .get_stan_model()
  
  ## ---- Construct datnew model matrix ----
  if(is.null(datnew)){
    X_new = matrix(rowMeans(X), ncol = 1)
  } else {
    X_new = stats::model.matrix(formula, datnew) |> t()
  }
  
  ## ---- Construct data list and arguments to pass to Stan ----
  datlist <- list(
    N = ncol(W),
    S_obs = nrow(W),
    K = nrow(X),
    W = W,
    X = X,
    N_new = ncol(X_new),
    X_new = X_new
  )
  
  # copied the defaults for stan
  stan_args <- list(
    chains = 4, iter = 2000, warmup = floor(iter/2), thin = 1,
    seed = sample.int(.Machine$integer.max, 1), 
    init = 'random', check_data = TRUE, 
    sample_file = NULL, diagnostic_file = NULL, verbose = FALSE, 
    algorithm = c("NUTS", "HMC", "Fixed_param"),
    control = NULL, include = TRUE, 
    cores = getOption("mc.cores", 1L),
    open_progress = interactive() && !isatty(stdout()) &&
      !identical(Sys.getenv("RSTUDIO"), "1"),
    show_messages = TRUE
  )
  
  # change args based on user input
  stan_args <- modifyList(stan_args, as.list(...))
  
  # ---- Fitting the model ----
  fit <- with(
    c(list(object = mod, data = datlist),
      stan_args), 
    rstan::sampling()
  )
  
  ## ---- Extract generated quantities and estimate Qk^(s) ----
  P_draws <- as.matrix(rstan::extract(fit, pars = "P"))
  
  
  
}




