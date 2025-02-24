#' R6 class to construct basic VAR model
#' @description
#' VAR: basic VAR model
#' @details
#' The `VAR` class implements basic data manipulation, info extraction and VAR tools computation, which are inherited by other specific VAR classes.
#'
#' @section Fields:
#' @field T int, total time periods of data
#' @field cons bool, indicating whether include constant in model
#' @field T.est int, total time periods used for estimation
#' @field n.var int number of variables
#' @field p.lag int, lags of the VAR model
#' @field hor int, number of horizons interested in inference
#' @field data **named** matrix of dataframe, each column indicating one time series, the first column should be time stramp
#' @field var.names character vector, variable names
#' @field time character of Date vector, time stramp
#' @field Y numeric matrix, dependent variables
#' @field Y.start numeric matrix, first `p.lag` periods of data
#' @field X numeric matrix, dependent variables
#' @field beta numeric matrix, VAR reduced form coefficients
#' @field U reduced form VAR residuals
#' @field Sigma reduced form VAR
#' @field B the structural impact matrix of SVAR
#' @field eps the structural shocks
#' @field Psi the un-orthogonal IRF
#' @field IRF the orthogonal IRF
#' @field FEVD the forecast error-variance decomposition matrix
#' @field beta.set store related bootstrap values
#' @field U.set store related bootstrap values
#' @field Sigma.set store related bootstrap values
#' @field B.set store related bootstrap values
#' @field eps.set store related bootstrap values
#' @field Psi.set store related bootstrap values
#' @field IRF.set store related bootstrap values
#' @field FEVD.set store related bootstrap values
#'
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

    #' @description
    #' `VAR` class constructor
    #' @param data time series data frame or names matrix
    #' @param p.lag int, number of lags in model
    #' @param cons bool, contain constant or not in model
    #' @return  an R6 class, `VAR`
    #' @export
    initialize = function(data = NA, p.lag = NA, cons = 1) {
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
      self$Y.start <- self$data[1:(self$p.lag), ]
      cat("* VAR class initialized.\n")
      cat("-> variables: ", paste0(self$var.names, collapse = ", "), "\n")
      cat("-> time period: ", as.character(self$time[1]), " to ", as.character(self$time[self$T]), "\n")
      cat("-> model: ", as.character(self$p.lag), " lags of all variables\n")
    },
    #' @description
    #' Compute IRF
    #' @param hor number of horizons to compute
    #' @return nothing, store IRF in `self$IRF`
    #' @export
    IRF.compute = function(hor = NA_integer_) {
      # compute Psi and IRF using basic OLS result
      IRF.out <- IRF_compute(self$beta, self$B, hor, self$n.var, self$p.lag)
      self$hor <- hor
      self$Psi <- IRF.out$Psi
      self$IRF <- IRF.out$IRF
    },
    #' @description
    #' Compute forecast error variance decomposition
    #' @param hor number of horizons to compute
    #' @return nothing, store FEVD in `self$FEVD`
    #' @export
    FEVD.compute = function(hor = self$hor) {
      self$FEVD <- FEVD_compute(self$Sigma, self$B, self$Psi, self$n.var, self$hor)
    },
    #' @description
    #' compute historical decomposition
    #' @param start `as.Date()` format, start of interested period
    #' @param end `as.Date()` format, end of interested period
    #' @param which.var the variable you want to do historical decomposition
    #' @return historical decomposition matrix, each row denotes the contribution of one shock
    #' @export
    HDC.compute = function(start = NA, end = NA, which.var = NA_character_) {
      # compute HDC using basic OLS result, start and end should be as.Date() variables
      var.id <- which(self$var.names == which.var)
      start.id <- which(self$time == start) - self$p.lag
      end.id <- which(self$time == end) - self$p.lag
      HDC <- get_HDs_ts(start.id, end.id, var.id, self$IRF, self$eps)
      return(HDC)
    },
    #' @description
    #' compute IRFs of parameter sets
    #' @param hor horizon
    #' @param prob probability of CI
    #' @return a list, containing median, average, upper bound and lower bound of the CI
    #' @seealso [VAR$IRF.compute()]
    #' @export
    IRF.compute.batch = function(hor = NA_integer_, prob = NA_real_) {
      msg <- paste("* IRF computed and", prob, "HDP get, time usage")
      HDP <- c((1 - prob) / 2, 1 - (1 - prob) / 2)
      tic(msg)
      IRF.batch <- IRF_compute_batch(self$beta.set, self$B.set, hor)
      self$Psi.set <- IRF.batch$Psi
      self$IRF.set <- IRF.batch$IRF
      IRF.med <- apply(self$IRF.set, c(1, 2), median)
      IRF.avg <- apply(self$IRF.set, c(1, 2), mean)
      IRF.ub <- apply(self$IRF.set, c(1, 2), function(x) quantile(x, probs = HDP[2]))
      IRF.lb <- apply(self$IRF.set, c(1, 2), function(x) quantile(x, probs = HDP[1]))
      toc()
      IRF.list <- list("avg" = IRF.avg, "med" = IRF.med, "ub" = IRF.ub, "lb" = IRF.lb)
      IRF.out <- lapply(IRF.list, mat2cube)
      return(IRF.out)
    },
    #' @description
    #' compute FEVDs of parameter sets
    #' @param hor horizon
    #' @param prob probability of CI
    #' @return a list, containing median, average, upper bound and lower bound of the CI
    #' @seealso [VAR$FEVD.compute()]
    #' @export
    FEVD.compute.batch = function(hor = NA_integer_, prob = NA_real_) {
      msg <- paste("* FEVD computed and", prob, "HDP get, time usage:")
      HDP <- c((1 - prob) / 2, 1 - (1 - prob) / 2)
      tic(msg)
      self$FEVD.set <- FEVD_compute_batch(self$beta.set, self$B.set, self$Sigma.set, self$Psi.set, hor)
      FEVD.avg <- apply(self$FEVD.set, c(1, 2), mean)
      FEVD.med <- apply(self$FEVD.set, c(1, 2), median)
      FEVD.ub <- apply(self$FEVD.set, c(1, 2), function(x) quantile(x, probs = HDP[2]))
      FEVD.lb <- apply(self$FEVD.set, c(1, 2), function(x) quantile(x, probs = HDP[1]))
      toc()
      FEVD.list <- list("avg" = FEVD.avg, "med" = FEVD.med, "ub" = FEVD.ub, "lb" = FEVD.lb)
      FEVD.out <- lapply(FEVD.list, mat2cube)
      return(FEVD.out)
    },
    #' @description
    #' compute HDCs of parameter sets
    #' @param prob probability of CI
    #' @param start start of the period of interest
    #' @param end end of the period of interest
    #' @param which.var the variable you want to do historical decomposition
    #' @return a list, containing median, average, upper bound and lower bound of the CI
    #' @seealso [VAR$HDC.compute()]
    #' @export
    HDC.compute.batch = function(start = NA, end = NA, which.var = NA_character_, prob = NA_real_) {
      # note: IRF and FEVD computes all, while HDC only compute the historical decomposition of "which.var"
      var.id <- which(self$var.names == which.var)
      start.id <- which(self$time == start) - self$p.lag
      end.id <- which(self$time == end) - self$p.lag
      msg <- paste("* HDC of", which.var, "computed and", prob, "HDP get, time usage:")
      HDP <- c((1 - prob) / 2, 1 - (1 - prob) / 2)
      tic(msg)
      HDC <- HDs_compute_batch(start.id, end.id, var.id, self$IRF.set, self$eps.set)
      HDC.avg <- apply(HDC, c(1, 2), mean)
      HDC.med <- apply(HDC, c(1, 2), median)
      HDC.ub <- apply(HDC, c(1, 2), function(x) quantile(x, probs = HDP[2]))
      HDC.lb <- apply(HDC, c(1, 2), function(x) quantile(x, probs = HDP[1]))
      toc()
      HDC.list <- list("avg" = HDC.avg, "med" = HDC.med, "ub" = HDC.ub, "lb" = HDC.lb)
      HDC.out <- lapply(HDC.list, mat2cube)
      return(HDC.out)
    }
  )
)
