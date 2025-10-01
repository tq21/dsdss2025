sim_data_cate <- function(n, grid=NULL) {
  UY <- rnorm(n, 0, 1)

  W1 <- round(runif(n, -1, 1), 2)
  W2 <- round(runif(n, -1, 1), 2)
  W3 <- round(runif(n, -1, 1), 2)

  A <- rbinom(n, 1, plogis(-0.25*W1+W2))

  # CATE
  tau <- function(W1, W2) return(0.5*W1+sin(2*pi*W2))
  #tau <- function(W1, W2) return(0.5*W1+sin(0.8*pi*W2))
  if (is.null(grid)) {
    Y <- -0.5+W1+0.5*W2+0.3*W3+A*tau(W1, W2)+UY
    return(data.frame(W1 = W1,
                      W2 = W2,
                      W3 = W3,
                      A = A,
                      Y = Y))
  } else {
    return(tau(grid$W1, grid$W2))
  }
}
