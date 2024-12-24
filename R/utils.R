#' @param Y input time series data matrix
#' @param lag.max a integer indicating max lag to test
#' @param criterion "AIC" or "BIC"
#' @return optimal lag based on AIC or BIC
which.lag <- function(data, lag.max = 12, criterion = "AIC") {
  T <- dim(data)[1]
  n.var <- dim(data)[2]
  AIC <- vector("numeric", lag.max)
  BIC <- vector("numeric", lag.max)
  for (p in 1:lag.max) {
    T.est <- T - p
    Y <- data[(p + 1):T, ]
    Y.lag <- data[p:(T - 1), ]
    if (p >= 2) {
      for (j in 2:p) {
        Y.lag <- rbind(Y.lag, data[(p - j + 1):(T - j), ])
      }
    }
    X <- rbind(matrix(1, T.est, 1), Y.lag)
    beta <- solve((t(X) %*% X), t(X) %*% Y)
    U <- Y - X %*% beta
    Sigma.tilde <- (t(U) %*% U) / T.est
    AIC[p] <- log(det(Sigma.tilde)) + 2 * n.var * (n.var * p + 1) / T.est
    BIC[p] <- log(det(Sigma.tilde)) + log(T.est) * n.var * (n.var * p + 1) / T.est
  }
  if (criterion == "AIC") {
    p.out <- which(AIC == min(AIC))
  } else {
    p.out <- which(BIC == min(BIC))
  }
  return(p.out)
}

mat2cube <- function(mat) {
  n.var <- sqrt(dim(mat)[1])
  hor <- dim(mat)[2]
  cube <- array(data = mat, dim = c(n.var, n.var, hor))
  cube.t <- aperm(cube, c(2, 1, 3))
  return(cube.t)
}
