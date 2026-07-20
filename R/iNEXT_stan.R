
.get_stan_model <- function() {
  if (!requireNamespace("rstan", quietly = TRUE))
    stop("Package 'rstan' is required. Install it with install.packages('rstan').")
  stan_file <- system.file("stan", "detection_model.stan", package = "iNEXT")
  rstan::stan_model(stan_file, auto_write = TRUE)
}


iNEXT_stan <- function(W, formula, data, q = c(0, 1, 2), datnew = NULL,
                       use_ndraws = 200L, ...){

  # ---- Validate inputs ----
  if (any(!(W == 0L | W == 1L))) {
    stop("x must be a raw incidence matrix with S_obs rows and N columns.\n
         Found some elements not equal to 0 or 1")
  }
  
  # ---- Model matrix construction ----
  X <- t(stats::model.matrix(formula, data))
  if (ncol(X) != ncol(W)) {
    stop(paste(
      "Dimension mismatch. Found",
      ncol(X),
      "sites/replicates in data, but",
      ncol(W), "in W."
    ))
  }
  
  ## ---- X_new model matrix ----
  if (is.null(datnew)) {
    X_new <- matrix(rowMeans(X), ncol = 1)
  } else {
    X_new <- t(stats::model.matrix(formula, datnew))
  }
  N_new <- ncol(X_new)

  # ---- Resolve sampling args before building the list to pass to Stan ----
  dots <- list(...)
  iter <- if (!is.null(dots$iter)){ dots$iter } else 2000L
  warmup <- if (!is.null(dots$warmup)){ dots$warmup } else floor(iter / 2)
  chains <- if (!is.null(dots$chains)){ dots$chains } else 4L
  seed <- if (!is.null(dots$seed)){ dots$seed } else sample.int(.Machine$integer.max, 1)

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
  t_obs_max <- 3L * ncol(W)
  datlist <- list(
    N = ncol(W),
    S_obs = nrow(W),
    K = nrow(X),
    W = W,
    X = X,
    N_new = N_new,
    X_new = X_new,
    t_obs_max = t_obs_max
  )

  # ---- Fit ----
  mod <- .get_stan_model()
  fit <- do.call(
    rstan::sampling,
    c(list(object = mod, data = datlist), stan_args)
  )

  # ---- Extract Qt_new and subsample posterior draws ----
  # Qt_draws: array [n_draws_total, t_obs_max+1, N_new]
  Qt_draws <- rstan::extract(fit, pars = "Qt_new")$Qt_new

  n_draws_total <- dim(Qt_draws)[1]
  use_ndraws <- min(use_ndraws, n_draws_total)
  draw_idx <- sort(sample.int(n_draws_total, use_ndraws))
  Qt_draws <- Qt_draws[draw_idx, , , drop = FALSE]

  # Truncate each draw x design-point to its effective sample size.
  # Qt_draws[i, 1, k] = t_obs_int, so the valid incidence-frequency vector
  # spans indices 1:(t_obs_int + 1).
  Qt_draws_list <- lapply(
    X = 1:use_ndraws,
    function(i, A) {
      lapply(
        1:N_new,
        function(k){
          return( A[i, 1:(A[i, 1, k] + 1), k] )
        }
      )
    },
    A = Qt_draws
  )

  # ---- Build per-draw accumulation curves for each design point ----
  # Common t grid spanning rarefaction through 2x suggested max
  t_pred_median <- sapply(
    Qt_draws_list,
    function(x){ length(x[[1]]) - 1 }
  ) |> median()
  t_grid <- unique(round(seq(1, 2 * t_pred_median, length.out = 40)))

  pt_names <- if (!is.null(datnew) && !is.null(rownames(datnew))) {
    rownames(datnew)
  } else {
    paste0("design_pt_", 1:N_new)
  }

  curves <- lapply(1:N_new, function(k) {
    draw_frames <- lapply(1:use_ndraws, function(i) {
      qt_ik <- Qt_draws_list[[i]][[k]]
      qD <- tryCatch(
        TD.m.est_inc(y = qt_ik, t_ = t_grid, qs = q),
        error = function(e) rep(NA_real_, length(t_grid) * length(q))
      )
      data.frame(
        draw = draw_idx[i],
        t = rep(t_grid, length(q)),
        order_q = as.factor(rep(q, each = length(t_grid))),
        qD = qD
      )
    })
    do.call(rbind, draw_frames)
  })
  names(curves) <- pt_names

  return( list(stanfit = fit, curves = curves) )
}
