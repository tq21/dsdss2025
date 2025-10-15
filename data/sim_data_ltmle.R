sim_data <- function(n,
                     K = 4,
                     a_bar = NULL) {

  # baseline covariates
  W1 <- rnorm(n)
  W2 <- rbinom(n, 1, plogis(W1))
  dt <- data.frame(W1 = W1, W2 = W2)

  # initialize histories
  L_hist <- list(); A_hist <- list(); C_hist <- list(); Y_hist <- list()

  L_prev <- rnorm(n, mean = W1, sd = 1)
  if (is.null(a_bar)) {
    A_prev <- rbinom(n, 1, plogis(-0.25 + 0.4*W1 + 0.4*W2 + 0.4*L_prev))
  } else {
    A_prev <- rep(a_bar[1], n)
  }

  # generate over time
  for (t in 1:K) {
    # past summaries up to t-1
    past_L <- if (t == 1) 0 else Reduce(`+`, lapply(1:(t-1), function(s) L_hist[[s]]))
    past_A <- if (t == 1) 0 else Reduce(`+`, lapply(1:(t-1), function(s) A_hist[[s]]))
    past_Y <- if (t == 1) 0 else Reduce(`+`, lapply(1:(t-1), function(s) Y_hist[[s]]))

    # L_t
    mu_Lt <- 0.25*W1 + 0.25*W2 + 0.30*past_L + 0.20*past_A + 0.15*past_Y
    L_t <- rnorm(n, mean = mu_Lt + 0.3*L_prev + 0.2*A_prev, sd = 1)

    # A_t
    if (is.null(a_bar)) {
      A_t <- rbinom(n, 1, plogis(-0.2 + 0.35*W1 + 0.15*W2 + 0.6*L_prev + 0.5*L_t))
    } else {
      A_t <- rep(a_bar[t], n)
    }

    # C_t: 1 = observed, 0 = censored
    if (is.null(a_bar)) {
      C_t <- rbinom(n, 1, plogis(2.6 - 0.4*W2 - 0.25*past_L - 0.25*past_A - 0.20*past_Y - 0.20*L_t - 0.20*A_t))
    } else {
      C_t <- rep(1, n)
    }

    # Y_t
    pY_t <- plogis(-1.0 + 0.7*A_t + 0.6*L_t + 0.25*W1 + 0.10*past_L + 0.10*past_A + 0.10*past_Y)
    Y_t <- rbinom(n, 1, pY_t)

    # write columns
    dt[[paste0("L", t)]] <- L_t
    dt[[paste0("A", t)]] <- A_t
    dt[[paste0("C", t)]] <- C_t
    dt[[paste0("Y", t)]] <- Y_t

    # store histories
    L_hist[[t]] <- L_t; A_hist[[t]] <- A_t; Y_hist[[t]] <- Y_t

    # update lags
    L_prev <- L_t; A_prev <- A_t
  }

  # apply monotone censoring
  if (is.null(a_bar)) {
    Cnames <- paste0("C", 1:K)
    Ynames <- paste0("Y", 1:K)
    Cmat <- as.matrix(dt[, Cnames, drop = FALSE])

    # cumprod row-wise: becomes 0 from the first censoring time onward
    cumobs <- t(apply(Cmat, 1, cumprod))

    for (t in 1:K) {
      drop_idx <- cumobs[, t] == 0
      if (any(drop_idx)) {
        dt[drop_idx, Ynames[t:K]] <- NA
        if (t < K) {
          dt[drop_idx, Cnames[(t+1):K]] <- NA
          dt[drop_idx, paste0("L", (t+1):K)] <- NA
          dt[drop_idx, paste0("A", (t+1):K)] <- NA
        }
      }
    }
  }

  return(dt)
}

get_truth <- function(n_large = 1e7, K = 4) {
  data_A1 <- sim_data(n_large, K = K, a_bar = rep(1, K))
  data_A0 <- sim_data(n_large, K = K, a_bar = rep(0, K))
  mean_A1 <- mean(data_A1[[paste0("Y", K)]])
  mean_A0 <- mean(data_A0[[paste0("Y", K)]])
  risk_diff <- mean_A1 - mean_A0
  risk_ratio <- mean_A1 / mean_A0

  return(list(mean_A1 = mean_A1,
              mean_A0 = mean_A0,
              risk_diff = risk_diff,
              risk_ratio = risk_ratio) )
}
