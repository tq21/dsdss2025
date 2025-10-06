hinge <- function(x, u) pmax(x - u, 0)

sim_data <- function(n,
                     A_counter = NULL) {

  # baseline covariates
  W1 <- rnorm(n, 0, 1)
  W2 <- rbinom(n, 1, 0.5)
  W3 <- rnorm(n, 0, 1)

  # treatment mechanism (only depends on W1 and W2)
  g1W <- plogis(0.1-0.6*W1+0.8*W2+0.3*W1*W2)
  if (is.null(A_counter)) {
    A <- rbinom(n, 1, g1W)
  } else {
    A <- rep(A_counter, n)
  }

  mu0 <- 1-0.25*W1^2+0.35*W1+0.4*W2+0.6*sin(pi*W3/2)
  tau <- 0.5+0.2*(W1 > 0)+0.4*hinge(W3, 0.2)
  Y <- mu0+tau*A+rnorm(n, 0, 1)

  return(data.frame(W1 = W1,
                    W2 = W2,
                    W3 = W3,
                    A = A,
                    Y = Y))
}

get_truth <- function(n_large=1e7) {
  data_A1 <- sim_data(n_large, A_counter=1)
  data_A0 <- sim_data(n_large, A_counter=0)
  return(mean(data_A1$Y-data_A0$Y))
}
