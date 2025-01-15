#' R6 class of classic VAR
#' @description
#' cVAR: classic VAR model
#' @details
#' The `cVAR` class implements classic VAR estimation an identification.
#' * estimation: OLS
#' * identification: recursive and IV
#' * inference: bootstrap
#'
#' @useDynLib VARS
#' @import Rcpp
#' @import RcppArmadillo
#' @import R6
#' @import stats
#' @import tictoc
#' @export
cVAR <- R6::R6Class(
  "Classic VAR",
  inherit = VAR,
  public = list(
    #' @field rank store rank of OLS `XX'` matrix
    rank = NULL,
    #' @field est.method store name of estimation method, should be "OLS"
    est.method = NULL,
    #' @field identify.method store name of identify method, should be "recursive" or "IV"
    identify.method = NULL,
    #' @field IV store IV data series, should be named matrix
    IV = NULL,
    #' @field IV.provide bool, indicate whether IV is provided
    IV.provide = NULL,
    #' @description
    #' `cVAR` class constructor
    #'
    #' @param data time series data frame or names matrix
    #' @param p.lag int, number of lags in model
    #' @param data.IV named matrix, the IV series; the IV name should **be set as the same as the related variable name**
    #' @param cons bool, contain constant or not in model
    #'
    #' @return  an R6 class, `cVAR`
    #' @export
    initialize = function(data = NA, p.lag = NA, data.IV = NA, cons = 1) {
      super$initialize(data, p.lag, cons)
      self$rank <- (det(t(self$X) %*% self$X) != 0)
      self$IV <- matrix(0, self$T.est, self$n.var)
      # write a function to extract IV
      if (any(!is.na(data.IV))){
        IV.names <- colnames(data.IV)
        colnames(self$IV) <- self$var.names
        self$IV[, IV.names] = as.matrix(data.IV[(self$p.lag + 1):self$T, ])
        self$IV.provide <- sort(which(self$var.names == IV.names))
        cat("* IV for variable", self$IV.provide, "provided.\n")
      }
    },
    #' @description
    #' estimate reduced form VAR using OLS
    #'
    #' @param method name of estimation method, currently only support "OLS"
    #' @return nothing, estimation result saved in `self$beta`, `self$U` and `self$Sigma`
    #' @export
    est = function(method = "OLS") {
      self$est.method == method
      if (!self$rank) stop("the matrix X'X not invertible")
      if (method == "OLS") {
        self$beta <- solve(t(self$X) %*% self$X, t(self$X) %*% self$Y)
        self$U <- self$Y - self$X %*% self$beta
        self$Sigma <- (t(self$U) %*% self$U) / (self$T.est)
        cat("* reduced form VAR parameters estimated using OLS.\n")
      }
    },
    #' @description
    #' idenitfy classic VAR model using "recursive" or "IV" method, then do bootstrap
    #' @param method identification method, should be "recursive" or "IV"
    #' @param boot.method bootstrap method, should be "wild" or "block"
    #' @param repeats bootstrap repeat times
    #'
    #' @return nothing, results saved in `self$B`, `self$eps` and `self$*.set`
    #' @export
    identify = function(method = NA_character_, boot.method = NA_character_, repeats = 1000) {
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
          xx <- self$IV[, i]
          test <- summary(lm(yy ~ xx, na.action = na.exclude))
          t.stat <- test$coefficients[2, 3]
          if (t.stat < 1.96) stop("Weak IV !")
        }
        self$B <- (t(self$U) %*% self$IV) / self$T.est
        B.modify <- diag(self$B)
        B.modify[which(B.modify == 0)] = 1
        B.div <- matrix(rep(B.modify, each = self$n.var), self$n.var, self$n.var)
        self$B <- self$B / B.div
        # B.Tinv <- solve(t(self$B)) # BUG how could eps be computed when only first col of B is identified ?
        # self$eps <- self$U %*% B.Tinv
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
      self$beta.set <- boot.out$beta.set
      self$U.set <- boot.out$U.set
      self$Sigma.set <- boot.out$Sigma.set
      self$B.set <- boot.out$B.set
      self$eps.set <- boot.out$eps.set
      toc()
    },
    #' @description
    #' compute VAR tools, currently support IRF, FEVD and HDC
    #' @param tool name of VAR tools to be computed
    #' @param hor horizons
    #' @param prob probability of CI
    #' @param start start of interested period, only required when `tool=="HDC"'
    #' @param end end of interested period, only required when `tool=="HDC"'
    #' @param which.var the variable name to do HDC, only required when `tool=="HDC"'
    #'
    #' @return related VAR tools, returned as a list with it's average, median and upper lower bound values
    #' @export
    tool = function(tool = NA_character_, hor = NA_integer_, prob = NA_real_, start = NA_integer_, end = NA_integer_, which.var = NA_character_) {
      if (tool == "IRF") {
        self$IRF.compute(hor = hor)
        IRF.out <- self$IRF.compute.batch(hor = hor, prob = prob)
        IRF.out$base <- self$IRF
        return(IRF.out)
      } else if (tool == "FECD") {
        self$FEVD.compute(hor = hor)
        FEVD.out <- self$FEVD.compute.batch(hor = hor, prob = prob)
        FEVD.out$base <- self$FEVD
        return(FEVD.out)
      } else if (tool == "HDC") {
        HDC.base <- self$HDC.compute(start, end, which.var)
        HDC.boot <- self$HDC.compute.batch(start, end, which.var, prob)
        HDC.boot$base <- HDC.base
        return(HDC.out)
      } else {
        stop("* Please specify correct VAR tool name.")
      }
    }
  )
)

