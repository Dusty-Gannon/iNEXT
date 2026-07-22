
.create_data_info <- function(incfreq_draws_list, use_ndraws, N_new,
                              pt_names, formula, data, datnew, conf_probs) {
  used_vars <- all.vars(formula)
  if (length(used_vars) == 0L) {
    cov_df <- data.frame(matrix(nrow = N_new, ncol = 0L))
  } else if (is.null(datnew)) {
    cov_df <- as.data.frame(lapply(data[used_vars], mean))
    cov_df <- cov_df[rep(1L, N_new), , drop = FALSE]
  } else {
    cov_df <- as.data.frame(datnew[, used_vars, drop = FALSE])
  }
  rownames(cov_df) <- NULL

  .summ <- function(x, digits) {
    round(c(med = median(x),
            lo = unname(quantile(x, conf_probs[1L])),
            hi = unname(quantile(x, conf_probs[2L]))),
          digits)
  }

  out <- do.call(rbind, lapply(seq_len(N_new), function(k) {

    t_vals <- vapply(seq_len(use_ndraws),
                     function(i) incfreq_draws_list[[i]][[k]][1L],
                     numeric(1L))
    s_vals <- vapply(seq_len(use_ndraws),
                     function(i) sum(incfreq_draws_list[[i]][[k]][-1L] > 0L),
                     numeric(1L))
    sc_vals <- vapply(seq_len(use_ndraws), function(i) {
      y <- incfreq_draws_list[[i]][[k]]
      Chat.Sam(y, y[1L])
    }, numeric(1L))
    # Qk_mat[j, i] = number of species with exactly j detections in draw i
    Qk_mat <- vapply(seq_len(use_ndraws), function(i) {
      counts <- incfreq_draws_list[[i]][[k]][-1L]
      vapply(1:10, function(j) sum(counts == j), integer(1L))
    }, integer(10L))

    t_s <- .summ(t_vals, digits = 0L)
    s_s <- .summ(s_vals, digits = 0L)
    sc_s <- .summ(sc_vals, digits = 4L)
    # apply over rows (one per Q1..Q10) -> 3 x 10 matrix (rows: med/lo/hi)
    Qk_s <- apply(Qk_mat, 1L, .summ, digits = 0L)
    # as.vector reads column-major -> Q1_med, Q1_lo, Q1_hi, Q2_med, ...
    Qk_flat <- setNames(
      as.list(as.vector(Qk_s)),
      paste0(rep(paste0("Q", 1:10), each = 3L), c("_med", "_lo", "_hi"))
    )

    row_k <- cbind(
      data.frame(Assemblage = pt_names[k], stringsAsFactors = FALSE),
      cov_df[k, , drop = FALSE],
      data.frame(
        t_obs_med = t_s["med"], t_obs_lo = t_s["lo"], t_obs_hi = t_s["hi"],
        S.obs_med = s_s["med"], S.obs_lo = s_s["lo"], S.obs_hi = s_s["hi"],
        SC_med = sc_s["med"], SC_lo = sc_s["lo"], SC_hi = sc_s["hi"],
        stringsAsFactors = FALSE
      ),
      as.data.frame(Qk_flat, stringsAsFactors = FALSE)
    )
    rownames(row_k) <- NULL
    row_k
  }))
  rownames(out) <- NULL
  out
}


# curves[[k]] rows are stacked draw blocks, each of length n_tq = length(t_grid) * length(q).
# Within each block: t cycles fastest inside q (matching TD.m.est_inc output order).
# SC is stored rep(sc_t, times = length(q)) so it is identical across q blocks.
.create_size_based <- function(curves, incfreq_draws_list, use_ndraws, N_new,
                               pt_names, t_grid, q, conf_probs) {
  n_t <- length(t_grid)
  n_tq <- n_t * length(q)

  do.call(rbind, lapply(seq_len(N_new), function(k) {

    T_k <- round(median(vapply(seq_len(use_ndraws),
                               function(i) incfreq_draws_list[[i]][[k]][1L],
                               numeric(1L))))

    crv_k <- curves[[k]]

    # Each draw's n_tq rows are contiguous; reshape into [use_ndraws x n_tq] matrices
    qD_mat <- matrix(crv_k$qD, nrow = use_ndraws, ncol = n_tq, byrow = TRUE)
    # SC is per-t (same across q); pull the first n_t entries from each draw block
    sc_idx <- rep((seq_len(use_ndraws) - 1L) * n_tq, each = n_t) +
              rep(seq_len(n_t), times = use_ndraws)
    sc_mat <- matrix(crv_k$SC[sc_idx], nrow = use_ndraws, ncol = n_t, byrow = TRUE)

    qD_mean <- colMeans(qD_mat, na.rm = TRUE)
    qD_lcl <- apply(qD_mat, 2L, quantile, probs = conf_probs[1L], na.rm = TRUE)
    qD_ucl <- apply(qD_mat, 2L, quantile, probs = conf_probs[2L], na.rm = TRUE)
    sc_mean <- colMeans(sc_mat, na.rm = TRUE)
    sc_lcl <- apply(sc_mat, 2L, quantile, probs = conf_probs[1L], na.rm = TRUE)
    sc_ucl <- apply(sc_mat, 2L, quantile, probs = conf_probs[2L], na.rm = TRUE)

    t_rep <- rep(t_grid, times = length(q))
    q_rep <- rep(q, each = n_t)

    out <- data.frame(
      Assemblage = pt_names[k],
      t = t_rep,
      Method = ifelse(t_rep < T_k, "Rarefaction",
               ifelse(t_rep == T_k, "Observed", "Extrapolation")),
      Order.q = q_rep,
      qD = qD_mean,
      qD.LCL = pmax(0, qD_lcl),
      qD.UCL = qD_ucl,
      SC = rep(sc_mean, times = length(q)),
      SC.LCL = pmax(0, rep(sc_lcl, times = length(q))),
      SC.UCL = pmin(1, rep(sc_ucl, times = length(q))),
      stringsAsFactors = FALSE
    )
    rownames(out) <- NULL
    out
  }))
}


# Each incfreq vector is a posterior predictive draw that integrates over both
# parameter uncertainty (the MCMC draw of beta) and sampling process uncertainty
# (bernoulli_rng in generated quantities). Calling invChat.Sam per draw therefore
# propagates both sources through to the diversity CIs, making them broader than
# standard iNEXT bootstrap CIs which treat the detection model as fixed.
.create_coverage_based <- function(incfreq_draws_list, use_ndraws, N_new,
                                   pt_names, size_based, q, conf_probs) {
  do.call(rbind, lapply(seq_len(N_new), function(k) {

    T_k <- round(median(vapply(seq_len(use_ndraws),
                               function(i) incfreq_draws_list[[i]][[k]][1L],
                               numeric(1L))))

    # SC grid from posterior-mean coverage already in size_based
    goalSC <- unique(size_based$SC[size_based$Assemblage == pt_names[k] &
                                   size_based$Order.q == q[1]])
    goalSC <- goalSC[goalSC > 0 & goalSC < 1]
    if (length(goalSC) == 0) return(NULL)

    n_sg <- length(goalSC) * length(q)
    qD_mat <- matrix(NA_real_, nrow = use_ndraws, ncol = n_sg)
    t_mat <- matrix(NA_real_, nrow = use_ndraws, ncol = n_sg)

    for (i in seq_len(use_ndraws)) {
      y_ik <- incfreq_draws_list[[i]][[k]]
      res <- tryCatch(
        invChat.Sam(y_ik, q, goalSC),
        error = function(e) NULL
      )
      if (!is.null(res)) {
        qD_mat[i, ] <- res$qD
        t_mat[i, ] <- res$t
      }
    }

    qD_mean <- colMeans(qD_mat, na.rm = TRUE)
    qD_lcl <- apply(qD_mat, 2L, quantile, probs = conf_probs[1L], na.rm = TRUE)
    qD_ucl <- apply(qD_mat, 2L, quantile, probs = conf_probs[2L], na.rm = TRUE)
    t_mean <- round(colMeans(t_mat, na.rm = TRUE))

    # invChat.Sam output order: all goalSC for q[1], then all goalSC for q[2], ...
    sc_rep <- rep(goalSC, times = length(q))
    q_rep <- rep(q, each = length(goalSC))

    out <- data.frame(
      Assemblage = pt_names[k],
      SC = sc_rep,
      t = t_mean,
      Method = ifelse(t_mean < T_k, "Rarefaction",
               ifelse(t_mean == T_k, "Observed", "Extrapolation")),
      Order.q = q_rep,
      qD = pmax(0, qD_mean),
      qD.LCL = pmax(0, qD_lcl),
      qD.UCL = qD_ucl,
      stringsAsFactors = FALSE
    )
    rownames(out) <- NULL
    out
  }))
}


.create_asy_est <- function(curves, incfreq_draws_list, use_ndraws, N_new,
                            pt_names, t_grid, q, conf_probs) {
  t_max <- max(t_grid)
  all_div_names <- c("0" = "Species richness",
                     "1" = "Shannon diversity",
                     "2" = "Simpson diversity")
  div_names <- all_div_names[as.character(q)]
  div_names[is.na(div_names)] <- paste0("Hill q=", q[is.na(div_names)])

  do.call(rbind, lapply(seq_len(N_new), function(k) {

    T_k <- round(median(vapply(seq_len(use_ndraws),
                               function(i) incfreq_draws_list[[i]][[k]][1L],
                               numeric(1L))))
    t_ref <- t_grid[which.min(abs(t_grid - T_k))]
    crv_k <- curves[[k]]

    do.call(rbind, lapply(seq_along(q), function(j) {
      q_sel <- as.numeric(as.character(crv_k$order_q)) == q[j]

      obs_draws <- crv_k$qD[crv_k$t == t_ref & q_sel]
      est_draws <- crv_k$qD[crv_k$t == t_max & q_sel]

      obs_val <- mean(obs_draws, na.rm = TRUE)
      est_val <- median(est_draws, na.rm = TRUE)
      se_val <- sd(est_draws, na.rm = TRUE)
      lcl_val <- unname(quantile(est_draws, probs = conf_probs[1L], na.rm = TRUE))
      ucl_val <- unname(quantile(est_draws, probs = conf_probs[2L], na.rm = TRUE))
      if (q[j] == 0) lcl_val <- max(lcl_val, obs_val)

      data.frame(
        Assemblage = pt_names[k],
        Diversity = div_names[j],
        Observed = round(obs_val, 3),
        Estimator = round(est_val, 3),
        s.e. = round(se_val, 3),
        LCL = round(lcl_val, 3),
        UCL = round(ucl_val, 3),
        stringsAsFactors = FALSE
      )
    }))
  }))
}


.get_stan_model <- function() {
  if (!requireNamespace("rstan", quietly = TRUE))
    stop("Package 'rstan' is required. Install it with install.packages('rstan').")
  stan_file <- system.file("stan", "detection_model.stan", package = "iNEXT")
  rstan::stan_model(stan_file, auto_write = TRUE)
}


iNEXT_stan <- function(W, formula, data, q = c(0, 1, 2), datnew = NULL,
                       use_ndraws = 200L, conf = 0.95, ...){

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

  # ---- Extract incfreq_new and subsample posterior draws ----
  # incfreq_draws: array [n_draws_total, S_obs + 1, N_new]
  # [draw, 1, k]     = t_obs_int[k]
  # [draw, s + 1, k] = detection count for species s at design point k
  # This is exactly the incidence-freq format expected by TD.m.est_inc and Chat.Sam.
  incfreq_draws <- rstan::extract(fit, pars = "incfreq_new")$incfreq_new

  n_draws_total <- dim(incfreq_draws)[1]
  use_ndraws <- min(use_ndraws, n_draws_total)
  draw_idx <- sort(sample.int(n_draws_total, use_ndraws))
  incfreq_draws <- incfreq_draws[draw_idx, , , drop = FALSE]

  # incfreq_draws_list[[i]][[k]] = c(t_obs_int, sp1_count, ..., spS_obs_count)
  incfreq_draws_list <- lapply(
    X = seq_len(use_ndraws),
    function(i, A) lapply(seq_len(N_new), function(k) A[i, , k]),
    A = incfreq_draws
  )

  # ---- Build per-draw accumulation curves for each design point ----
  # t_obs_int is always the first element of each incfreq vector
  t_pred_median <- median(sapply(
    seq_len(use_ndraws),
    function(i) incfreq_draws_list[[i]][[1]][1]
  ))
  t_grid <- unique(round(seq(1, 2 * t_pred_median, length.out = 40)))

  pt_names <- if (!is.null(datnew) && !is.null(rownames(datnew))) {
    rownames(datnew)
  } else {
    paste0("design_pt_", 1:N_new)
  }

  curves <- lapply(1:N_new, function(k) {
    draw_frames <- lapply(1:use_ndraws, function(i) {
      y_ik <- incfreq_draws_list[[i]][[k]]
      qD <- tryCatch(
        TD.m.est_inc(y = y_ik, t_ = t_grid, qs = q),
        error = function(e) rep(NA_real_, length(t_grid) * length(q))
      )
      sc <- tryCatch(
        Chat.Sam(y_ik, t_grid),
        error = function(e) rep(NA_real_, length(t_grid))
      )
      data.frame(
        draw = draw_idx[i],
        t = rep(t_grid, length(q)),
        order_q = as.factor(rep(q, each = length(t_grid))),
        qD = qD,
        SC = rep(sc, times = length(q))
      )
    })
    do.call(rbind, draw_frames)
  })
  names(curves) <- pt_names

  # ---- Build DataInfo, iNextEst, and AsyEst ----
  conf_probs <- c((1 - conf) / 2, 1 - (1 - conf) / 2)
  DataInfo <- .create_data_info(incfreq_draws_list, use_ndraws, N_new,
                                pt_names, formula, data, datnew, conf_probs)
  size_based <- .create_size_based(curves, incfreq_draws_list, use_ndraws, N_new,
                                   pt_names, t_grid, q, conf_probs)
  coverage_based <- .create_coverage_based(incfreq_draws_list, use_ndraws, N_new,
                                           pt_names, size_based, q, conf_probs)
  AsyEst <- .create_asy_est(curves, incfreq_draws_list, use_ndraws, N_new,
                             pt_names, t_grid, q, conf_probs)

  out <- list(DataInfo = DataInfo,
              iNextEst = list(size_based = size_based,
                              coverage_based = coverage_based),
              AsyEst = AsyEst,
              stanfit = fit,
              curves = curves)
  class(out) <- c("iNEXT_stan", "iNEXT")
  out
}
