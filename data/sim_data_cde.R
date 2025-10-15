sim_data_cde <- function(n, grid=NULL) {
  gen_A <- function(W1, W2, dens=FALSE) {
    mu_A <- 1+0.3*W1-1.5*W2
    sigma_A <- 2+0.3*abs(W1)
    if (dens) {
      return(dnorm(grid$A, mean = mu_A, sd = sigma_A))
    } else {
      return(rnorm(length(W1), mean = mu_A, sd = sigma_A))
    }
  }

  if (is.null(grid)) {
    W1 <- rnorm(n, mean = 0, sd = 1)
    W2 <- rnorm(n, mean = 0, sd = 1)
    A <- gen_A(W1, W2)
  } else {
    W1 <- grid$W1
    W2 <- grid$W2
    A <- gen_A(W1, W2, dens=TRUE)
  }

  return(data.frame(W1 = W1,
                    W2 = W2,
                    A = A))
}
