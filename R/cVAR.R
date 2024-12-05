cVAR <- R6::R6Class(
  "Classic VAR",
  inherit = VAR,
  public = list(
    # leave place for other variables
    rank = NULL,
    est.method = NULL,
    identify.method = NULL,
    IV = NULL,
    IV.provide = NULL,
    initialize = function(data = NA, p.lag = NA, data.IV = NA, cons = 1) {
      super$initialize(data, p.lag, cons)
      self$rank <- (det(tcrossprod(self$X)) != 0)
      self$IV <- matrix(0, self$T.est, self$n.var)
      # write a function to extract IV
      if (any(!is.na(data.IV))){
        IV.names <- colnames(data.IV)
        self$IV[which(self$var.names == IV.names)] = as.matrix(data.IV[(self$p.lag + 1):self$T, ])
        self$IV.provide <- sort(which(self$var.names == IV.names))
        cat("* IV for variable", self$IV.provide, "provided.\n")
      }
    },
    est = function(method = "OLS") {
      self$est.method == method
      if (!rank) stop("the matrix X'X not invertible")
      if (method == "OLS") {
        self$beta <- solve(t(self$X) %*% self$X, t(self$X) %*% self$Y)
        self$U <- self$Y - self$X %*% self$beta
        self$Sigma <- (t(self$U) %*% self$U) / (self$T.est)
        cat("* reduced form VAR parameters estimated using OLS.\n")
      }
    },
    identify = function(method = "recursive", boot.method = "block", repeats = 1000) {
      self$identify.method <- method
      ## identify section
      if (method == "recursive") {
        self$B <- t(chol(self$Sigma)) # to get B as lower triangular factor
        B.Tinv <- solve(t(self$B)) # Note here eps is actually eps' and U is U', and eps'=u'B'^{-1}
        self$eps <- self$U %*% B.Tinv # use B matrix and reduced-form shock to compute structural shock
        cat("* model identified using ", method, " approach.\n")
      } else if (method == "IV") {
        for (i in self$IV.provide) {
          yy <- self$U[, i]
          xx <- IV[, i]
          test <- summary(lm(yy ~ xx, na.action = na.exclude()))
          t.stat <- test$coefficients[2, 3]
          if (t.stat < 1.96) stop("Weak IV !")
        }
        self$B <- (t(self$U) %*% self$IV) / self$T.est
        B.modify <- diag(self$B)
        B.modify[which(B.modify == 0)] = 1
        self$B <- self$B / B.modify
        B.Tinv <- solve(t(self$B)) # BUG how could eps be computed when only first col of B is identified ?
        self$eps <- self$U %*% B.Tinv
        cat("* model identified using ", method, " approach.\n")
        cat("* Note: only column ", self$IV.provide, " of B matrix and IRF are credible.\n")
      } else {
        stop("current classic VAR only support recursive and IV identification!")
      }
      ## bootstrap section
      cat("* start bootstrapping.\n")
      msg <- paste("* inference using", boot.method, "bootstrap method,", repeats, "repeats time usage")
      tic(msg)
      boot.out <- bootstrap_c(self$Y.start, self$beta, self$U, self$IV, method, boot.method, repeats, self$T)
      beta.set <- boot.out$beta.set
      U.set <- boot.out$U.set
      Sigma.set <- boot.out$Sigma.set
      B.set <- boot.out$B.set
      eps.set <- boot.out$eps.set
      toc()
    }
  )
)

