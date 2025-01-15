#' Find optimal lag of VAR according to certain criterion
#'
#' @param Y input time series data matrix
#' @param lag.max a integer indicating max lag to test
#' @param criterion "AIC" or "BIC"
#'
#' @return optimal lag based on AIC or BIC
#' @export
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

#' Transforming IRF matrix to 3D array
#'
#' @param mat The IRF matrix, row variable shock combination, col horizon
#'
#' @returns a 3D array, with each slice the IRF matrix of that horizon
#' @seealso [IRF.compute.batch()]
#' @export
mat2cube <- function(mat) {
  n.var <- sqrt(dim(mat)[1])
  hor <- dim(mat)[2]
  cube <- array(data = mat, dim = c(n.var, n.var, hor))
  cube.t <- aperm(cube, c(2, 1, 3))
  return(cube.t)
}


#' Plot IRF
#'
#' @param obj a cVAR or bVAR object, provide full variable names
#' @param IRFs IRF computed from cVAR or bVAR
#' @param hor horizons you want to plot
#' @param shk.names which shocks you want to include
#' @param var.names which variables you want to include
#'
#' @returns a figure with sub figues, each containing an IRF of 1 variable, 1 shock
#' @seealso [bVAR$tool()], [cVAR$tool()]
#' @export
IRF.plot <- function(obj, IRFs, hor, shk.names, var.names) {
  num.var <- length(var.names)
  num.shk <- length(shk.names)
  mynames <- obj$var.names
  idx.var <- match(var.names, mynames)
  idx.shk <- match(shk.names, mynames)
  num.IRF <- length(IRFs)
  colors <- c("blue", "red")
  shadows <- c("#C0C0C0AE", "#FFC0CBAB")

  par(mfrow = c(num.var, num.shk))
  # par(family="serif", cex.main = 1.5, cex.lab = 1.2, cex.axis=1.2)
  for (i in idx.var) {
    for (j in idx.shk) {
      shk.name <- mynames[j]
      var.name <- mynames[i]
      for (k in 1:num.IRF) {
        myIRF <- IRFs[[k]]
        base <- myIRF$base[i, j, ]
        ub <- myIRF$ub[i, j, ]
        lb <- myIRF$lb[i, j, ]
        xaxis <- 1:(hor+1)
        y.range <- range(c(base, lb, ub))
        plot(xaxis, base,
             type = "l", col = colors[k], ylim = y.range, lwd=3,
             xlab = "Horizons", ylab = "Percent", main = paste(var.name,"to",shk.name,"shock")
        )
        polygon(c(xaxis, rev(xaxis)), c(ub, rev(lb)), col = shadows[k], border = NA)
        lines(xaxis, base, col = colors[k], lty = 1, lwd = 3)
        abline(h = 0, lwd = 3, lty = 2)
      }
    }
  }
}
