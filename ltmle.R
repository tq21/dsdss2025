library(ltmle)

rexpit <- function(x) rbinom(n=length(x), size=1, prob=plogis(x))
n <- 10000
W <- rnorm(n)
A <- rexpit(W)
Y <- rexpit(W + A)
data <- data.frame(W, A, Y)
head(data)

result <- ltmle(data, Anodes = "A", Ynodes = "Y", abar = list(1, 0))
result

n.large <- 1e6
W <- rnorm(n.large)
A <- 1
Y1 <- rexpit(W + A)
A <- 0
Y0 <- rexpit(W + A)
mean(Y1 - Y0)



n <- 1000
W1 <- rnorm(n)
W2 <- rbinom(n, size=1, prob=0.3)
W3 <- rnorm(n)
A <- rexpit(-1 + 2 * W1 + W3)
Y <- rexpit(-0.5 + 2 * W1^2 + 0.5 * W2 - 0.5 * A + 0.2 * W3 * A - 1.1 * W3)
data <- data.frame(W1, W2, W3, A, Y)

result <- ltmle(data, Anodes="A", Lnodes=NULL, Ynodes="Y", abar=1,
                SL.library="default")
summary(result)
summary(result, estimator="iptw")


result <- ltmle(data, Anodes="A", Lnodes=NULL, Ynodes="Y",
                Qform=c(Y="Q.kplus1 ~ I(W1^2) + W2 + W3*A"),
                gform="A ~ W1 + W3", abar=1, SL.library="default")
summary(result)

result <- ltmle(data, Anodes="A", Lnodes=NULL, Ynodes="Y",
                Qform=c(Y="Q.kplus1 ~ I(W1^2) + W2 + W3*A"), gform="A ~ W1 + W3",
                abar=1, SL.library=NULL)
summary(result)

result <- ltmle(data, Anodes="A", Lnodes=NULL, Ynodes="Y",
                abar=list(1, 0), SL.library="default")
summary(result)


# censoring
n <- 100000
W <- rnorm(n)
C <- BinaryToCensoring(is.censored = rexpit(W))
summary(C)
Y <- rep(NA, n)
Y[C == "uncensored"] <- rexpit(W[C == "uncensored"])
data <- data.frame(W, C, Y)
head(data, 20)
result <- ltmle(data, Anodes = NULL, Cnodes = "C", Ynodes = "Y", abar = NULL)
summary(result)



# longitudinal data
n <- 1000
W <- rnorm(n)
A1 <- rexpit(W)
L <- 0.3 * W + 0.2 * A1 + rnorm(n)
A2 <- rexpit(W + A1 + L)
Y <- rexpit(W - 0.6 * A1 + L - 0.8 * A2)
data <- data.frame(W, A1, L, A2, Y)
head(data)

ltmle(data, Anodes=c("A1", "A2"), Lnodes="L", Ynodes="Y", abar=c(0, 0))


# longitudinal data w censoring
n <- 1000
W <- rnorm(n)
A1 <- rexpit(W)
C <- BinaryToCensoring(is.censored = rexpit(0.6 * W - 0.5 * A1))
uncensored <- C == "uncensored"
L <- A2 <- Y <- rep(NA, n)
L[uncensored] <- (0.3 * W[uncensored] + 0.2 * A1[uncensored] + rnorm(sum(uncensored)))
A2[uncensored] <- rexpit(W[uncensored] + A1[uncensored] + L[uncensored])
Y[uncensored] <- rexpit(W[uncensored] - 0.6 * A1[uncensored] + L[uncensored] - 0.8 * A2[uncensored])
data <- data.frame(W, A1, C, L, A2, Y)
head(data)


# dynamic regime
abar <- matrix(nrow=n, ncol=2)
abar[, 1] <- 1
abar[, 2] <- L > 0

result.abar <- ltmle(data, Anodes=c("A1", "A2"), Cnodes = "C", Lnodes="L", Ynodes="Y", abar=abar)
summary(result.abar)


rule <- function(row) c(1, row["L"] > 0)
result.rule <- ltmle(data, Anodes=c("A1", "A2"), Cnodes = "C", Lnodes="L", Ynodes="Y", rule=rule)
summary(result.rule)
# rule is a function applied to each row of data which returns a numeric vector
# of the same length as Anodes





ltmle(data, Anodes=c("A1", "A2"), Cnodes = "C", Lnodes="L", Ynodes="Y", abar=c(1, 0))
