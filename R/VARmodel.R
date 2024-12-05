#' VAR: basic VAR model
#'
#' The `VAR` class implements basic data manipulation, info extraction and VAR tools computation, which are inherited by other specific VAR classes.
#'
#' @field T int, total time periods of data
#' @field cons bool, indicating whether include constant in model
#' @field T.est int, total time periods used for estimation
#' @field p.lag int, lags of the VAR model
#' @field hor int, number of horizons interested in inference
#' @field data **named** matrix of dataframe, each column indicating one time series, the first column should be time stramp
#' @field var.names character vector, variable names
#' @field time character of Date vector, time stramp
#' @field Y numeric matrix, dependent variables
#' @field Y.start numeric matrix, first `p.lag` periods of data
#' @field X numeric matrix, dependent variables
#' @field beta numeric matrix, VAR reduced form coefficients
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(name, value)}}{Creates a new instance of the class.}
#'   \item{\code{greet()}}{Prints a greeting message.}
#'   \item{\code{add(x, y)}}{Adds two numbers and returns the result.}
#' }
#' @useDynLib VARS
#' @import Rcpp
#' @import RcppArmadillo
#' @import R6
#' @import stats
#' @import tictoc
#' @export
VAR <- R6::R6Class(
  "basic VAR model",
  public = list(
    T = NULL,
    cons = NULL,
    T.est = NULL,
    n.var = NULL,
    p.lag = NULL,
    hor = NULL,
    data = NULL,
    var.names = NULL,
    time = NULL,
    Y = NULL,
    Y.start = NULL,
    X = NULL,
    beta = NULL,
    beta.set = NULL,
    U = NULL,
    U.set = NULL,
    Sigma = NULL,
    Sigma.set = NULL,
    B = NULL,
    B.set = NULL,
    eps = NULL,
    eps.set = NULL,
    Psi = NULL,
    Psi.set = NULL,
    IRF = NULL,
    IRF.set = NULL,
    FEVD = NULL,
    FEVD.set = NULL,
    initialize = function(data = NA, p.lag = NA, cons = 1) {
      # TODO check data is a matrix, p.lag is a positive integer, add a `validate()` function
      self$T <- dim(data)[1]
      self$cons <- cons
      self$n.var <- dim(data)[2] - 1
      self$time <- data[, 1]
      self$data <- as.matrix(data[, 2:(self$n.var + 1)])
      self$p.lag <- p.lag
      self$T.est <- self$T - p.lag
      self$var.names <- colnames(self$data)
      self$Y <- self$data[(1 + p.lag):self$T, ]
      Y.lag <- self$data[p.lag:(self$T - 1), ]
      for (i in 2:p.lag) {
        Y.lag <- cbind(Y.lag, self$data[(1 + p.lag - i):(self$T - i), ])
      }
      if (self$cons == 1) {
        self$X <- cbind(rep(1, self$T.est), Y.lag)
      } else {
        self$X <- Y.lag
      }
      self$Y.start <- self$data[1:(self$T - p.lag), ]
      cat("* VAR class initialized.\n")
      cat("-> variables: ", paste0(self$var.names, collapse = ", "), "\n")
      cat("-> time period: ", as.character(self$time[1]), " to ", as.character(self$time[self$T]), "\n")
      cat("-> model: ", as.character(self$p.lag), " lags of all variables\n")
    },
    IRF.compute = function(hor = NA_integer_) {
      # compute Psi and IRF using basic OLS result
      IRF.out <- IRF_compute(self$beta, self$B, hor, self$n.var, self$p.lag, self$cons)
      self$hor <- hor
      self$Psi <- IRF.out$Psi
      self$IRF <- IRF.out$IRF
    },
    FEVD.compute = function(hor = self$hor) {
      # compute FEVD using basic OLS result
      self$FEVD <- FEVD_compute(self$Sigma, self$B, self$Psi, self$n.var, self$hor)
    },
    HDC.compute = function(start = NA, end = NA, which.var = NA_character_) {
      # compute HDC using basic OLS result, start and end should be as.Date() variables
      var.id <- which(self$var.names == which.var)
      start.id <- which(self$time == start) - self$p.lag
      end.id <- which(self$time == end) - self$p.lag
      HDC <- get_HDs_ts(start.id, end.id, var.id, self$IRF, self$eps)
      return(HDC)
    },
    IRF.compute.batch = function(hor = NA_integer_, prob = NA_real_) {
      msg <- paste("* IRF computed and", prob, "HDP get, time usage:")
      HDP <- c((1 - prob) / 2, 1 - (1 - prob) / 2)
      tic(msg)
      self$Psi.set <- IRF_compute_batch(self$beta.set, self$B.set, hor, self$cons)$Psi
      self$IRF.set <- IRF_compute_batch(self$beta.set, self$B.set, hor, self$cons)$IRF
      IRF.avg <- apply(self$IRF.set, c(1, 2), median)
      IRF.ub <- apply(self$IRF.set, c(1, 2), function(x) quantile(x, probs = HDP[2]))
      IRF.lb <- apply(self$IRF.set, c(1, 2), function(x) quantile(x, probs = HDP[1]))
      toc()
      return(list("avg" = IRF.avg, "ub" = IRF.ub, "lb" = IRF.lb))
    },
    FEVD.compute.batch = function(hor = NA_integer_, prob = NA_real_) {
      msg <- paste("* FEVD computed and", prob, "HDP get, time usage:")
      HDP <- c((1 - prob) / 2, 1 - (1 - prob) / 2)
      tic(msg)
      self$FEVD.set <- FEVD_compute_batch(self$beta.set, self$B.set, self$Sigma.set, self$Psi.set, hor)
      FEVD.avg <- apply(self$FEVD.set, c(1, 2), median)
      FEVD.ub <- apply(self$FEVD.set, c(1, 2), function(x) quantile(x, probs = HDP[2]))
      FEVD.lb <- apply(self$FEVD.set, c(1, 2), function(x) quantile(x, probs = HDP[1]))
      toc()
      return(list("avg" = FEVD.avg, "ub" = FEVD.ub, "lb" = FEVD.lb))
    },
    HDC.compute.batch = function(start = NA, end = NA, which.var = NA_character_) {
      # note: IRF and FEVD computes all, while HDC only compute the historical decomposition of "which.var"
      var.id <- which(self$var.names == which.var)
      start.id <- which(self$time == start) - self$p.lag
      end.id <- which(self$time == end) - self$p.lag
      msg <- paste("* HDC of", which.var, "computed and", prob, "HDP get, time usage:")
      HDP <- c((1 - prob) / 2, 1 - (1 - prob) / 2)
      tic(msg)
      HDC <- HDs_compute_batch(start.id, end.id, var.id, self$IRF.set, self$eps.set)
      HDC.avg <- apply(HDC, c(1, 2), median)
      HDC.ub <- apply(HDC, c(1, 2), function(x) quantile(x, probs = HDP[2]))
      HDC.lb <- apply(HDC, c(1, 2), function(x) quantile(x, probs = HDP[1]))
      toc()
      return(list("avg" = HDC.avg, "ub" = HDC.ub, "lb" = HDC.lb))
    }
  )
)
