
.get_stan_model <- function() {
  if (!requireNamespace("rstan", quietly = TRUE))
    stop("Package 'rstan' is required. Install it with install.packages('rstan').")
  stan_file <- system.file("stan", "detection_model.stan", package = "iNEXT")
  rstan::stan_model(stan_file, auto_write = TRUE)
}


iNEXT_stan <- function(W, formula, data, q = c(0, 1, 2), datnew = NULL, ...){

  # ---- Validate inputs ----
  if (any(!(W == 0L | W == 1L)))
    stop("x must be a raw incidence matrix with S_obs rows and N columns.\n
         Found some elements not equal to 0 or 1")

  # ---- Model matrix construction ----
  X <- t(stats::model.matrix(formula, data))
  if (ncol(X) != ncol(W))
    stop(paste(
      "Dimension mismatch. Found",
      ncol(X),
      "sites/replicates in data, but",
      ncol(W), "in W."
    ))

  # ---- Prediction design ----
  if (is.null(datnew)) {
    X_new <- matrix(rowMeans(X), ncol = 1)
  } else {
    X_new <- t(stats::model.matrix(formula, datnew))
  }

  # ---- Resolve sampling args before building the list ----
  dots <- list(...)
  iter <- if (!is.null(dots$iter)){ dots$iter } else 2000L
  warmup  <- if (!is.null(dots$warmup)){ dots$warmup } else floor(iter / 2)
  chains  <- if (!is.null(dots$chains)){ dots$chains } else 4L
  seed    <- if (!is.null(dots$seed)){ dots$seed } else sample.int(.Machine$integer.max, 1)

  stan_args <- list(
    chains = chains,
    iter = iter,
    warmup = warmup,
    thin = if (!is.null(dots$thin)){ dots$thin } else 1L,
    seed = seed,
    cores = if (!is.null(dots$cores)){ dots$cores } else getOption("mc.cores", 1L),
    verbose = if (!is.null(dots$verbose)){ dots$verbose } else FALSE
  )

  # ---- Data list for Stan ----
  datlist <- list(
    N = ncol(W),
    S_obs = nrow(W),
    K = nrow(X),
    W = W,
    X = X,
    N_new = ncol(X_new),
    X_new = X_new,
    t_obs_max = 3L * ncol(W)
  )

  # ---- Fit ----
  mod <- .get_stan_model()
  fit <- do.call(
    rstan::sampling,
    c(list(object = mod, data = datlist), stan_args)
  )

  # ---- Extract Qt_new and pipe into iNEXT pipeline ----
  # Qt_new draws: array [n_draws, t_obs_max+1, N_new]
  Qt_draws <- rstan::extract(fit, pars = "Qt_new")$Qt_new
  
  # convert this to a list of truncated Qt_new per design point
  Qt_draws_list <- lapply(
    X = 1:dim(Qt_draws)[1],
    function(i, A) {
      lapply(
        1:dim(Qt_draws)[3],
        function(k){
          return( A[i, 1:(A[i, 1, k] + 1), k] )
        }
      )
    },
    A = Qt_draws
  )
  
}
