#' @param n RCT + RWD sample size
#' @param rct_prop Proportion of RCT sample in combined data
sim_data <- function(n,
                     rct_prop = 0.4,
                     A_counter = NULL) {

  # baseline covariates (with mild correlation)
  Sigma <- matrix(c(1, 0.2, 0.1, 0.2,
                    0.2, 1, 0.15, 0.25,
                    0.1, 0.15, 1, 0.1,
                    0.2, 0.25, 0.1, 1), 4, 4)
  W <- MASS::mvrnorm(n, mu = rep(0, 4), Sigma = Sigma)
  W1 <- W[,1]
  W2 <- W[,2]
  W3 <- rbinom(n, 1, plogis(0.5*W[,3]))
  W4 <- 0.7*W[,4] + 0.3*sin(W[,1])

  # RCT participation
  lp_S <- -0.2 + 0.25*W1 - 0.25*W2 + 0.4*W3 - 0.2*W4
  intercept_adj <- uniroot(function(b0) mean(plogis(b0 + lp_S)) - rct_prop,
                           interval = c(-2, 2))$root
  S <- rbinom(n, 1, plogis(intercept_adj + lp_S))

  # treatment assignment:
  # - In RCT (S=1): randomized A ~ Bern(0.5)
  # - In RWD (S=0): depends on W1 and W2
  if (is.null(A_counter)) {
    g_rct <- rep(0.5, n)
    g_rwd <- plogis(-0.3 + 0.8*W1 - 0.9*W2 + 0.4*W1*W2)
    g <- ifelse(S == 1, g_rct, g_rwd)
    A <- rbinom(n, 1, g)
  } else {
    A <- rep(A_counter, n)
  }

  # outcome
  mu_base <- -0.3 + 0.9*W1 - 0.6*W2 + 0.3*W3 + 0.2*W4 + 0.2*W1*W4
  site_shift <- (S == 0)*(0.9+0.2*W1)*A
  tau <- 0.15 + 0.6*W1 - 0.4*W2 + 0.3*W3 - 0.2*W4 + 0.3*sin(W1)
  eps <- rnorm(n, 0, 1)
  Y <- mu_base + site_shift + A * tau + eps

  if (is.null(A_counter)) {
    return(data.frame(W1, W2, W3, W4, S, A, Y))
  } else {
    return(data.frame(W1, W2, W3, W4, S, A, Y, tau))
  }
}

get_truth <- function(n_large=1e8) {
  data_large <- sim_data(n_large, A_counter=1)
  return(mean(data_large$tau))
}
